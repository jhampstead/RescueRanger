#include "annotation.h"

#include <string.h>
#include <htslib/faidx.h>
#include <htslib/regidx.h>

#include "bcftools/gff.h"
#include "logging/log.h"

int find_gf_tscripts_by_iregion(IRegion iregion, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts) {
    regidx_t *gff_tscripts = (regidx_t *) gff_get(gff, idx_tscript);
    regitr_t *itr = regitr_init(gff_tscripts);

    const char *chr = bcf_hdr_id2name(hdr, iregion.rid);
    if (!regidx_overlap(gff_tscripts, chr, iregion.start, iregion_last_base(iregion), itr)) { regitr_destroy(itr); return -1; } // The index takes an inclusive end

    gf_tscripts *tscripts = malloc(sizeof(gf_tscripts));
    kv_init(*tscripts);

    while (regitr_overlap(itr)) {
        kv_push(gf_tscript_t *, *tscripts, regitr_payload(itr, gf_tscript_t *));
    }

    regitr_destroy(itr);
    *out_gf_tscripts = tscripts;

    return 0;
}

void match_gf_tscripts_to_vep_tscripts(gf_tscripts *tscripts, Record *r, gff_t *gff, gf_tscripts **out_gf_scripts) {
    gf_tscripts *record_tscripts = malloc(sizeof(gf_tscripts));
    kv_init(*record_tscripts);

    for (int i = 0; i < r->vep->n; i++) {
        int has_matched = 0;

        if (tscripts == NULL) { kv_push(gf_tscript_t *, *record_tscripts, NULL); continue; }
        for (int j = 0; j < tscripts->n; j++) {
            gf_tscript_t *t = kv_A(*tscripts, j);
            const char *tid = gff_id2string(gff, 0, t->id);

            if (strcmp(tid, kv_A(*r->vep, i)->tid) != 0) continue;
            if (t->type != GF_PROTEIN_CODING) break; // Transcript is not coding, so should not be matched and set to NULL instead

            has_matched = 1;
            kv_push(gf_tscript_t *, *record_tscripts, t);
            break;
        }

        if (!has_matched) kv_push(gf_tscript_t *, *record_tscripts, NULL);
    }

    *out_gf_scripts = record_tscripts;
}

int find_gf_tscripts_for_record(Record *r, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts) {
    IRegion region = variant_to_iregion(r->v);
    gf_tscripts *tscripts = NULL;
    find_gf_tscripts_by_iregion(region, hdr, gff, &tscripts);

    gf_tscripts *matched_tscripts; match_gf_tscripts_to_vep_tscripts(tscripts, r, gff, &matched_tscripts);
    if (tscripts != NULL) { kv_destroy(*tscripts); free(tscripts); }

    *out_gf_tscripts = matched_tscripts;
    return 0;
}

// The exon ranges of a GFF transcript, ascending and inclusive, as the GFF index returns them
static void get_gff_exons(gff_t *gff, const char *chr, const gf_tscript_t *gf_t, Ranges *out) {
    regidx_t *exon_index = (regidx_t *) gff_get(gff, idx_exon);
    regitr_t *itr = regitr_init(exon_index);

    if (regidx_overlap(exon_index, chr, gf_t->beg, gf_t->end, itr) == 0) {
        log_warn("No overlapping exons found for %s:%i-%i", chr, gf_t->beg, gf_t->end);
        regitr_destroy(itr);
        return;
    }

    while (regitr_overlap(itr)) {
        gf_exon_t *gf_e = regitr_payload(itr, gf_exon_t *);

        // Exons for other transcripts and genes can overlap, but we don't want them!
        if (gf_e->tr->id != gf_t->id) continue;

        kv_push(Range, *out, ((Range) { gf_e->beg, gf_e->end }));
    }

    regitr_destroy(itr);
}

static Transcript *build_transcript(gff_t *gff, const char *chr, const gf_tscript_t *gf_t) {
    Range *cds = calloc(gf_t->ncds, sizeof(Range));
    int *constitutive = calloc(gf_t->ncds, sizeof(int));
    for (size_t k = 0; k < gf_t->ncds; k++) {
        cds[k] = (Range) { gf_t->cds[k]->beg, gf_t->cds[k]->beg + gf_t->cds[k]->len - 1 }; // Inclusive, like the GFF
        constitutive[k] = gf_t->cds[k]->constitutive;
    }

    Ranges exons;
    kv_init(exons);
    get_gff_exons(gff, chr, gf_t, &exons);

    TranscriptData data = {
        .id = gf_t->id,
        .start = gf_t->beg,
        .end = gf_t->end,
        .strand = gf_t->strand,
        .cds = cds,
        .constitutive = constitutive,
        .ncds = gf_t->ncds,
        .exons = exons.a,
        .nexons = kv_size(exons)
    };
    Transcript *transcript = Transcript_build(&data);

    free(cds);
    free(constitutive);
    kv_destroy(exons);

    return transcript;
}

// Returns NULL, after logging, if the sequence cannot be fetched
static Gene *build_gene(faidx_t *fai, const char *chr, const gf_gene_t *gf_gene) {
    int slen = 0;
    char *seq = faidx_fetch_seq(fai, chr, gf_gene->beg, gf_gene->end, &slen); // faidx takes an inclusive end
    if (seq == NULL) {
        log_error("Failed to fetch the sequence of gene %s at %s:%u-%u", gf_gene->name, chr, gf_gene->beg, gf_gene->end);
        return NULL;
    }

    Gene *gene = malloc(sizeof(Gene));
    gene->gff_id = gf_gene->id;
    gene->name = gf_gene->name; // Borrowed from the GFF
    gene->coords = (Range) { gf_gene->beg, gf_gene->end + 1 }; // Half-open, the GFF end is inclusive
    gene->seq = (kstring_t) { slen, slen + 1, seq };
    gene->strand = gf_gene->strand;

    return gene;
}

static int annotate_record(Annotation *annotation, const Session *s, Record *r, RecordAnnotation *out) {
    if (r->vep == NULL) return 0;

    gf_tscripts *gf_tscripts_of_record = NULL;
    find_gf_tscripts_for_record(r, s->hdr, s->gff, &gf_tscripts_of_record);
    const char *chr = bcf_hdr_id2name(s->hdr, r->v->rid);

    int status = 0;
    for (size_t j = 0; j < kv_size(*gf_tscripts_of_record); j++) {
        gf_tscript_t *gf_t = kv_A(*gf_tscripts_of_record, j);
        if (gf_t == NULL) {
            kv_push(Transcript *, out->tscripts, NULL);
            kv_push(Gene *, out->transcript_genes, NULL);
            continue;
        }

        khint_t k = kh_get_TranscriptHash(annotation->transcripts, gf_t->id);
        if (k == kh_end(annotation->transcripts)) {
            int is_new;
            k = kh_put_TranscriptHash(annotation->transcripts, gf_t->id, &is_new);
            kh_val(annotation->transcripts, k) = build_transcript(s->gff, chr, gf_t);
        }
        kv_push(Transcript *, out->tscripts, kh_val(annotation->transcripts, k));

        k = kh_get_GeneHash(annotation->genes, gf_t->gene->id);
        if (k == kh_end(annotation->genes)) {
            Gene *gene = build_gene(s->fai, chr, gf_t->gene);
            if (gene == NULL) { status = -1; break; }

            int is_new;
            k = kh_put_GeneHash(annotation->genes, gf_t->gene->id, &is_new);
            kh_val(annotation->genes, k) = gene;
        }

        Gene *gene = kh_val(annotation->genes, k);
        kv_push(Gene *, out->transcript_genes, gene);

        int has_gene = 0;
        for (size_t g = 0; g < kv_size(out->genes); g++) {
            if (kv_A(out->genes, g)->gff_id == gene->gff_id) { has_gene = 1; break; }
        }
        if (!has_gene) kv_push(Gene *, out->genes, gene);
    }

    kv_destroy(*gf_tscripts_of_record);
    free(gf_tscripts_of_record);

    return status;
}

Annotation *Annotation_build(const Records *records, const Session *s) {
    Annotation *annotation = malloc(sizeof(Annotation));
    annotation->transcripts = kh_init_TranscriptHash();
    annotation->genes = kh_init_GeneHash();
    annotation->n_records = kv_size(*records);
    annotation->records = calloc(annotation->n_records, sizeof(RecordAnnotation)); // All zero is an empty kvec

    for (size_t i = 0; i < annotation->n_records; i++) {
        if (annotate_record(annotation, s, kv_A(*records, i), &annotation->records[i]) == 0) continue;

        Annotation_destroy(annotation);
        return NULL;
    }

    return annotation;
}

void Annotation_destroy(Annotation *annotation) {
    if (annotation == NULL) return;

    for (khint_t k = kh_begin(annotation->transcripts); k != kh_end(annotation->transcripts); ++k) {
        if (!kh_exist(annotation->transcripts, k)) continue;
        Transcript_destroy(kh_val(annotation->transcripts, k));
    }
    kh_destroy_TranscriptHash(annotation->transcripts);

    for (khint_t k = kh_begin(annotation->genes); k != kh_end(annotation->genes); ++k) {
        if (!kh_exist(annotation->genes, k)) continue;

        Gene *gene = kh_val(annotation->genes, k);
        free(gene->seq.s); // The name is borrowed from the GFF, so it is not freed
        free(gene);
    }
    kh_destroy_GeneHash(annotation->genes);

    for (size_t i = 0; i < annotation->n_records; i++) {
        kv_destroy(annotation->records[i].tscripts);
        kv_destroy(annotation->records[i].transcript_genes);
        kv_destroy(annotation->records[i].genes);
    }
    free(annotation->records);

    free(annotation);
}
