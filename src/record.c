#include "record.h"
#include "klib/kvec.h"
#include "utils.h"
#include "triggers.h"
#include "logging/log.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

const int n_csq_fields = 29;
const char vcf_delim = ',';
const char *csq_field = "CSQ";
const char *bcftools_csq_field = "BCSQ";
const char *bcftools_transcript = "transcript"; // The column of bcftools csq with the ID of the transcript
const char bcftools_pointer_prefix = '@';       // An entry '@position' means the consequence is given at another record
const char bcftools_marker_prefix = '*';        // Marks a consequence that bcftools csq flags; the name after it is the consequence's
const char *hdr_desc = "Description";
const char vep_delim = '|';
const char *vep_desc_prefix = "Format: ";
const int vep_desc_prefix_len = 8;
const char consequence_delim = '&';

const char *vep_consequence = "Consequence";
const char *vep_feature_type = "Feature_type";
const char *vep_feature = "Feature";
const char *vep_feature_type_transcript = "Transcript";


void init_splice_site(SpliceSite **site, size_t *n, size_t *m) {
    *m = 1;
    *n = 0;
    *site = malloc(sizeof(SpliceSite));
}

void resize_splice_site(SpliceSite **site, size_t *m) {
    *m *= 2;
    *site = realloc(*site, *m * sizeof(SpliceSite));
}


/*
 *
 * WARN: Does not add last field if it is empty (e.g., 28 fields instead of 29 for VEP)
 */
int *str_split(kstring_t *s, char delim, int *num_splits) {
    int i = 0;
    kvec_t(int) field_idxs; kv_init(field_idxs);
    while (i < s->l) {
        kv_push(int, field_idxs, i);

        while(*(s->s + i) && *(s->s + i) != delim ) i++;

        s->s[i++] = '\0';
    }

    *num_splits = field_idxs.n;
    return field_idxs.a;
}

/*
*
* Usage: 
*   int id; khash_str2int_get(idHash, "Consequence", &id);
*   fields[id];
*/
const char *Annotation_default_tag(const AnnotationKind kind) {
    return kind == ANNOTATION_BCFTOOLS ? bcftools_csq_field : csq_field;
}

static const char *annotation_option_name(const AnnotationKind kind) {
    return kind == ANNOTATION_BCFTOOLS ? "bcftools" : "vep";
}

// The column of that name, or -1
static int column_of(kh_str2int_t *str2id, const char *name) { // Not const: khash_str2int_get takes a void *
    int id;
    return khash_str2int_get(str2id, name, &id) == 0 ? id : -1;
}

HdrFormat *parse_csq_hdr_fmt(bcf_hdr_t *hdr, const AnnotationKind kind, const char *tag) {
    if (tag == NULL || tag[0] == '\0') tag = Annotation_default_tag(kind);

    if (bcf_hdr_id2int(hdr, BCF_DT_ID, tag) < 0) {
        log_error("No INFO/%s field in the header of the variants; it is where the consequences (annotation %s) are looked for. "
                  "Use --csq-tag if the annotation is in another field", tag, annotation_option_name(kind));
        return NULL;
    }

    // Get the Description of the annotation, which ends with its columns
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, NULL, tag, NULL);
    const int key = hrec == NULL ? -1 : bcf_hrec_find_key(hrec, hdr_desc);
    const char *format = key < 0 ? NULL : strstr(hrec->vals[key], vep_desc_prefix);
    if (format == NULL) {
        log_error("The header of INFO/%s has no \"%s...\" in its Description, so its columns are unknown", tag, vep_desc_prefix);
        return NULL;
    }
    const char *description = format + vep_desc_prefix_len; // Skip "Format: "
    int len = strlen(description);
    if (len > 0 && description[len - 1] == '"') len--; // The closing quote of the Description

    // Split the columns
    kstring_t s = {0}; kputsn(description, len, &s);
    int *field_idxs, n; field_idxs = ksplit(&s, vep_delim, &n);

    kh_str2int_t *idHash = kh_init_str2int();
    char **fields = malloc(n * sizeof(char *));
    int ret;
    for (int i = 0; i < n; i++) {
        kh_put_str2int(idHash, s.s + field_idxs[i], &ret);
        khash_str2int_set(idHash, s.s + field_idxs[i], i);
        fields[i] = s.s + field_idxs[i];
    }
    free(field_idxs);

    HdrFormat *hdr_fmt = malloc(sizeof(HdrFormat));
    *hdr_fmt = (HdrFormat) {
        .data_ = s.s, .id2str = fields, .str2id = idHash, .kind = kind, .tag = strdup(tag),
        .consequence_id = column_of(idHash, vep_consequence),
        .feature_type_id = kind == ANNOTATION_VEP ? column_of(idHash, vep_feature_type) : -1,
        .feature_id = column_of(idHash, kind == ANNOTATION_VEP ? vep_feature : bcftools_transcript),
    };

    // The columns of the other kind are what a file of the other kind looks like
    const int missing = hdr_fmt->consequence_id < 0 || hdr_fmt->feature_id < 0 || (kind == ANNOTATION_VEP && hdr_fmt->feature_type_id < 0);
    if (missing) {
        log_error("INFO/%s does not have the columns of %s annotation (%s), its columns are: %s. Check --annotation and --csq-tag",
                  tag, annotation_option_name(kind),
                  kind == ANNOTATION_VEP ? "Consequence, Feature_type, Feature" : "Consequence, transcript", description);
        destroy_hdr_fmt(hdr_fmt);
        return NULL;
    }

    return hdr_fmt;
}

void destroy_hdr_fmt(HdrFormat *hdr_fmt) {
    kh_destroy_str2int(hdr_fmt->str2id);
    free(hdr_fmt->id2str);
    free(hdr_fmt->data_);
    free(hdr_fmt->tag);
    free(hdr_fmt);
}

// The consequences of bcftools csq that have another name in VEP, which is the one the triggers look for
static const struct { const char *bcftools, *vep; } bcftools_consequence_names[] = {
    { "frameshift", FRAMESHIFT_VARIANT },
    { "splice_donor", SPLICE_DONOR_VARIANT },
    { "splice_acceptor", SPLICE_ACCEPTOR_VARIANT },
    { "splice_region", SPLICE_REGION_VARIANT },
    { "non_coding", "non_coding_transcript_variant" },
};

// A consequence as VEP names it, in memory of its own
static char *consequence_as_vep(const AnnotationKind kind, const char *consequence) {
    if (kind != ANNOTATION_BCFTOOLS) return strdup(consequence);

    if (consequence[0] == bcftools_marker_prefix) consequence++;
    for (size_t i = 0; i < sizeof(bcftools_consequence_names) / sizeof(bcftools_consequence_names[0]); i++) {
        if (strcmp(consequence, bcftools_consequence_names[i].bcftools) == 0) return strdup(bcftools_consequence_names[i].vep);
    }

    return strdup(consequence);
}

// The field with that index of a split entry. An entry can have fewer fields than the format: an empty last field is not split off
// (see str_split), and an entry of bcftools that points to another record has only one. NULL then.
static const char *field_at(const kstring_t *entry, const int *field_idxs, const int n_fields, const int id) {
    if (id < 0 || id >= n_fields) return NULL;
    return entry->s + field_idxs[id];
}

// Adds an entry for the transcript, or with merge the consequences to the entry that the transcript already has. Takes cons.
static void push_transcript_csqs(VepCsqs *transcripts, const char *tid, char **cons, const int n_cons, const bool merge) {
    if (merge) {
        for (int i = 0; i < transcripts->n; i++) {
            TscriptCsqs *existing = kv_A(*transcripts, i);
            if (strcmp(existing->tid, tid) != 0) continue;

            existing->consequences = realloc(existing->consequences, (existing->num_consequences + n_cons + 1) * sizeof(char *));
            for (int j = 0; j < n_cons; j++) existing->consequences[existing->num_consequences++] = cons[j];
            free(cons);
            return;
        }
    }

    TscriptCsqs *t = calloc(1, sizeof(TscriptCsqs)); // No trigger yet
    t->tid = strdup(tid);
    t->num_consequences = n_cons;
    t->consequences = cons;
    kv_push(TscriptCsqs *, *transcripts, t);
}

VepCsqs *parse_vep_transcript_consequences(bcf1_t *v, const HdrFormat *vep_fmt, bcf_hdr_t *hdr) {
    char *vep_str = NULL; int n_vep_values = 0;
    const int n_chars = bcf_get_info_string(hdr, v, vep_fmt->tag, &vep_str, &n_vep_values);

    VepCsqs *transcripts = malloc(sizeof(VepCsqs)); kv_init(*transcripts);

    // A record without INFO/CSQ has no consequences, so no transcripts: the same as a record whose features are no transcripts
    if (n_chars <= 0 || vep_str == NULL) {
        free(vep_str);
        return transcripts;
    }

    kstring_t s = {0}; kputs(vep_str, &s);
    int *field_idxs, n; field_idxs = str_split(&s, vcf_delim, &n);
    for (int i = 0; i < n; i++) {
        int flen = strlen(s.s + field_idxs[i]);
        kstring_t transcript_csq = { flen, flen+1, s.s + field_idxs[i] }; // WARN: NOT MALLOC'ED, SHOULD NOT BE FREED
        int *vep_field_idxs, vep_n; vep_field_idxs = str_split(&transcript_csq, vep_delim, &vep_n);

        const char *consequence_field = field_at(&transcript_csq, vep_field_idxs, vep_n, vep_fmt->consequence_id);
        const char *tr_id = field_at(&transcript_csq, vep_field_idxs, vep_n, vep_fmt->feature_id);

        // VEP has features that are no transcript. bcftools csq has no such column: it has entries without a transcript (intergenic)
        // and entries that point to another record, which have no consequence of their own
        int skip = consequence_field == NULL || tr_id == NULL;
        if (!skip && vep_fmt->kind == ANNOTATION_VEP) {
            const char *feature_type = field_at(&transcript_csq, vep_field_idxs, vep_n, vep_fmt->feature_type_id);
            skip = feature_type == NULL || strcmp(feature_type, vep_feature_type_transcript) != 0;
        } else if (!skip) {
            skip = tr_id[0] == '\0' || consequence_field[0] == bcftools_pointer_prefix;
        }
        if (skip) { free(vep_field_idxs); continue; }

        // Parse consequences for transcript, vep allows multiple with '&' separators
        int clen = strlen(consequence_field);
        kstring_t consequences = { clen, clen+1, (char *) consequence_field };
        int *consequence_field_idxs, consequence_n; consequence_field_idxs = str_split(&consequences, consequence_delim, &consequence_n);
        char **cons = malloc((consequence_n + 1) * sizeof(char *));
        for (int j = 0; j < consequence_n; j++) {
            cons[j] = consequence_as_vep(vep_fmt->kind, consequences.s + consequence_field_idxs[j]);
        }

        // Entries of one transcript are joined for bcftools csq; VEP has one entry per transcript
        push_transcript_csqs(transcripts, tr_id, cons, consequence_n, vep_fmt->kind == ANNOTATION_BCFTOOLS);

        free(consequence_field_idxs);
        free(vep_field_idxs);
    }
    free(s.s);
    free(vep_str);
    free(field_idxs);

    return transcripts;
}

void destroy_vep_csqs(VepCsqs *vep_csqs) {
    for (int i = 0; i < vep_csqs->n; i++) {
        TscriptCsqs *vep_csq_field = kv_A(*vep_csqs, i);

        for (int j = 0; j < vep_csq_field->num_consequences; j++) {
            free(vep_csq_field->consequences[j]);
        }
        free(vep_csq_field->consequences);
        free(vep_csq_field->tid);

        free(vep_csq_field);
    }

    kv_destroy(*vep_csqs);
    free(vep_csqs);
}

Alleles Alleles_from_variant(const bcf1_t *v) {
    const char *ref = v->d.allele[0];
    const char *alt = v->n_allele > 1 ? v->d.allele[1] : ref + strlen(ref); // No ALT reads as an empty one

    int prefix_len = 0;
    while (ref[prefix_len] != '\0' && ref[prefix_len] == alt[prefix_len]) prefix_len++;

    return (Alleles) {
        v->pos + prefix_len,
        prefix_len,
        (int) strlen(ref) - prefix_len,
        (int) strlen(alt) - prefix_len,
        alt + prefix_len
    };
}

bool bcf1_variant_in_region(const bcf1_t *v, const Range region) {
    const Alleles a = Alleles_from_variant(v);
    return variant_in_region(&a, region);
}

Record *Record_init(bcf1_t *v, VepCsqs *csqs) {
    Record *r = malloc(sizeof(Record));

    bcf1_t *v_dup = bcf_dup(v);
    bcf_unpack(v_dup, BCF_UN_ALL);
    r->v = v_dup;
    r->alleles = Alleles_from_variant(v_dup);

    r->genes = NULL;
    r->tscripts = NULL;
    r->transcript_genes = NULL;
    r->vep = NULL;
    r->isoforms = NULL;
    r->flags = NULL;
    r->native_scores = NULL;
    r->diff_scores = NULL;
    r->haplotype_scores = NULL;

    r->skipped = 0;
    r->unphased = 0;

    r->vep = csqs;

    return r;
}

void set_record_flag(Record *r, int transcript_index, int phase_index, uint32_t flag_value) {
    r->flags[transcript_index * NUM_PHASES + phase_index] |= flag_value;
}

void Isoforms_destroy(Isoforms *isoforms) {
    kv_destroy(isoforms->isoforms);

    RescueProbabilities *probs = isoforms->rescue_type_probabilities;
    if (probs != NULL) {
        free(probs->n_isoforms);
        free(probs->rescues);
        free(probs->rescue_types);
        free(probs->rescue_probabilities);
        free(probs->affected_bases);
        free(probs);
    }
}

static void free_scores(Score **scores, const size_t n_slots) {
    for (size_t i = 0; i < n_slots; i++) free(scores[i]);
    free(scores);
}

void Record_destroy(Record *r) {
    bcf_destroy(r->v);

    if (r->vep != NULL) destroy_vep_csqs(r->vep);

    // The isoform array has a slot per transcript and phase, so it goes before the transcripts. Only built slots own memory.
    if (r->isoforms != NULL) {
        for (size_t i = 0; i < kv_size(*r->tscripts) * NUM_PHASES; i++) {
            if (r->isoforms[i].valid) Isoforms_destroy(&r->isoforms[i]);
        }
        free(r->isoforms);
    }

    free(r->flags);

    // The score arrays have a slot per gene and phase, so they go before the genes. If native scores are set, all three are
    // set. A slot is NULL until it is scored; the alt and gene a Score points at are not its own.
    if (r->native_scores != NULL) {
        const size_t n_slots = kv_size(*r->genes) * NUM_PHASES;
        free_scores(r->native_scores, n_slots);
        free_scores(r->diff_scores, n_slots);
        free_scores(r->haplotype_scores, n_slots);
    }

    if (r->tscripts != NULL) {
        kv_destroy(*r->tscripts);
        free(r->tscripts);
    }

    if (r->transcript_genes != NULL) {
        kv_destroy(*r->transcript_genes);
        free(r->transcript_genes);
    }

    if (r->genes != NULL) {
        kv_destroy(*r->genes);
        free(r->genes);
    }

    free(r);
}

Introns Introns_from_exons(const Range *exons, size_t n_exons) {
    Introns introns;
    kv_init(introns);

    if (n_exons < 2) return introns; // Also keeps n_exons - 1 from underflowing
    kv_resize(Intron, introns, n_exons - 1);

    int64_t coding_offset = 0;
    for (size_t i = 1; i < n_exons; i++) {
        int64_t upstream_end = exons[i-1].end;
        int64_t downstream_start = exons[i].start;
        assert(downstream_start >= upstream_end && "exons must be sorted ascending and must not overlap");

        Intron intron;
        intron.start = upstream_end;
        intron.end = downstream_start;
        intron.len = intron.end - intron.start;

        coding_offset += intron.len;
        intron.coding_offset = coding_offset;

        kv_push(Intron, introns, intron);
    }

    return introns;
}

int64_t Exon_donor_base(const Exon *exon, uint32_t strand) {
    return strand == STRAND_REV ? exon->ref_bounds->start : exon->ref_bounds->end - 1;
}

int64_t Exon_acceptor_base(const Exon *exon, uint32_t strand) {
    return strand == STRAND_REV ? exon->ref_bounds->end - 1 : exon->ref_bounds->start;
}

Range Intron_between_sites(uint32_t strand, int64_t donor_base, int64_t acceptor_base) {
    if (strand == STRAND_REV) return (Range) { acceptor_base + 1, donor_base };
    return (Range) { donor_base + 1, acceptor_base };
}

int64_t Intron_donor_base(Range intron, uint32_t strand) {
    return strand == STRAND_REV ? intron.end : intron.start - 1;
}

int64_t Intron_acceptor_base(Range intron, uint32_t strand) {
    return strand == STRAND_REV ? intron.start - 1 : intron.end;
}

Range Utr_range(enum UTRType type, uint32_t strand, Range gene, Range boundary_cds) {
    // The UTR on the low side of the gene is the 5' UTR on the forward strand and the 3' UTR on the reverse strand
    const bool low_side = (type == UTR5) == (strand != STRAND_REV);
    if (low_side) return (Range) { gene.start, boundary_cds.start };
    return (Range) { boundary_cds.end, gene.end };
}

bool Exon_splice_window(const Transcript *transcript, size_t exon_index, bool donor, Range *window) {
    if (transcript->strand == STRAND_UNK) return false;

    const size_t n_exons = kv_size(transcript->exons);
    const Range bounds = *kv_A(transcript->exons, exon_index)->ref_bounds;

    const bool upper_border = donor == (transcript->strand == STRAND_FWD);
    if (upper_border) {
        if (exon_index + 1 >= n_exons) return false;
        *window = (Range) { bounds.start + 1, kv_A(transcript->exons, exon_index + 1)->ref_bounds->start };
    } else {
        if (exon_index == 0) return false;
        *window = (Range) { kv_A(transcript->exons, exon_index - 1)->ref_bounds->end - 1, bounds.end - 2 };
    }

    return true;
}

Transcript *Transcript_build(const TranscriptData *data) {
    Transcript *t = calloc(1, sizeof(Transcript));
    t->id = data->id;
    t->start = data->start;
    t->end = data->end + 1;
    t->strand = data->strand;

    t->ncds = data->ncds;
    t->mcds = data->ncds;
    t->cds = data->ncds ? malloc(data->ncds * sizeof(Range)) : NULL;
    for (size_t k = 0; k < data->ncds; k++) {
        t->cds[k] = (Range) { data->cds[k].start, data->cds[k].end + 1 };
    }

    kv_init(t->exons);
    kv_init(t->introns);
    if (data->nexons == 0) return t;

    kv_resize(Exon *, t->exons, data->nexons);
    Range bounds[data->nexons];

    size_t ci = 0;
    for (size_t i = 0; i < data->nexons; i++) {
        Range *ref_bounds = malloc(sizeof(Range));
        *ref_bounds = (Range) { data->exons[i].start, data->exons[i].end + 1 };
        bounds[i] = *ref_bounds;

        Exon *exon = malloc(sizeof(Exon));
        kv_push(Exon *, t->exons, exon);

        // Does this exon have a CDS? If not, the coding region stays available for the next exon
        Range *cds = ci < data->ncds ? &t->cds[ci] : NULL;
        Range cds_intersection = {0};
        if (cds == NULL || !intersect_ranges(*ref_bounds, *cds, &cds_intersection)) {
            *exon = (Exon) { .cds = NULL, .ref_bounds = ref_bounds, .splice_sites = NULL, .constitutive = false, .in_frame = false };
            continue;
        }

        *exon = (Exon) { .cds = cds, .ref_bounds = ref_bounds, .splice_sites = NULL, .constitutive = data->constitutive[ci] != 0, .in_frame = range_len(*cds) % 3 == 0 };
        ci++;
    }

    t->introns = Introns_from_exons(bounds, data->nexons);

    return t;
}

void SpliceSites_destroy(SpliceSites *sites) {
    SpliceSite *arrays[4] = { sites->canonical_donors, sites->canonical_acceptors, sites->cryptic_donors, sites->cryptic_acceptors };
    size_t counts[4] = { sites->n_canonical_donors, sites->n_canonical_acceptors, sites->n_cryptic_donors, sites->n_cryptic_acceptors };

    for (int a = 0; a < 4; a++) {
        for (size_t i = 0; i < counts[a]; i++) {
            Records *records = (Records *) arrays[a][i].records_ptr;
            if (records == NULL) continue;

            kv_destroy(*records);
            free(records);
        }
        free(arrays[a]);
    }

    free(sites);
}

Transcript *Transcript_phase_copy(const Transcript *shared) {
    Transcript *copy = malloc(sizeof(Transcript));
    *copy = *shared; // cds and introns keep pointing at the shared data

    const size_t n_exons = kv_size(shared->exons);
    kv_init(copy->exons);
    if (n_exons > 0) kv_resize(Exon *, copy->exons, n_exons);

    for (size_t i = 0; i < n_exons; i++) {
        Exon *exon = malloc(sizeof(Exon));
        *exon = *kv_A(shared->exons, i); // ref_bounds and cds keep pointing at the shared data
        exon->splice_sites = NULL;
        kv_push(Exon *, copy->exons, exon);
    }

    return copy;
}

void Transcript_phase_destroy(Transcript *copy) {
    for (size_t i = 0; i < kv_size(copy->exons); i++) {
        Exon *exon = kv_A(copy->exons, i);

        if (exon->splice_sites != NULL) SpliceSites_destroy(exon->splice_sites);
        free(exon);
    }

    kv_destroy(copy->exons);
    free(copy);
}

NmdRegions Transcript_nmd_regions(const Transcript *transcript) {
    NmdRegions regions;
    kv_init(regions.start_proximal);
    kv_init(regions.long_exon);
    kv_init(regions.distal);

    const bool reverse = transcript->strand == STRAND_REV;
    const size_t n_exons = kv_size(transcript->exons);

    // The first 150 coding bases in transcript order, then listed in ascending order
    int64_t *taken = calloc(transcript->ncds ? transcript->ncds : 1, sizeof(int64_t));
    int64_t remaining = NMD_START_PROXIMAL_LENGTH;
    for (size_t k = 0; k < transcript->ncds && remaining > 0; k++) {
        const size_t ci = reverse ? transcript->ncds - 1 - k : k;
        taken[ci] = range_len(transcript->cds[ci]) < remaining ? range_len(transcript->cds[ci]) : remaining;
        remaining -= taken[ci];
    }
    for (size_t ci = 0; ci < transcript->ncds; ci++) {
        if (taken[ci] == 0) continue;

        const Range cds = transcript->cds[ci];
        kv_push(Range, regions.start_proximal, (reverse ? (Range) { cds.end - taken[ci], cds.end } : (Range) { cds.start, cds.start + taken[ci] }));
    }
    free(taken);

    // The part of every exon that is more than 400 bases before its last base. The last exon has no donor
    for (size_t ei = 0; ei < n_exons; ei++) {
        const bool is_last = reverse ? ei == 0 : ei == n_exons - 1;
        if (is_last) continue;

        const Range bounds = *kv_A(transcript->exons, ei)->ref_bounds;
        const int64_t length = range_len(bounds) - (NMD_LONG_EXON_LENGTH + 1);
        if (length <= 0) continue;

        kv_push(Range, regions.long_exon, (reverse ? (Range) { bounds.end - length, bounds.end } : (Range) { bounds.start, bounds.start + length }));
    }

    // The last exon, the intron before it and the last 50 bases of the exon before that make one stretch
    if (n_exons > 0) {
        const Range last = *kv_A(transcript->exons, reverse ? 0 : n_exons - 1)->ref_bounds;
        Range distal = last;

        if (n_exons > 1) {
            const Range penultimate = *kv_A(transcript->exons, reverse ? 1 : n_exons - 2)->ref_bounds;
            const int64_t tail = range_len(penultimate) < NMD_PENULTIMATE_EXON_LENGTH ? range_len(penultimate) : NMD_PENULTIMATE_EXON_LENGTH;

            if (reverse) distal.end = penultimate.start + tail;
            else distal.start = penultimate.end - tail;
        }

        kv_push(Range, regions.distal, distal);
    }

    return regions;
}

void NmdRegions_destroy(NmdRegions *regions) {
    kv_destroy(regions->start_proximal);
    kv_destroy(regions->long_exon);
    kv_destroy(regions->distal);
}

void Transcript_destroy(Transcript *transcript) {
    // Exon.cds points into transcript->cds, and splice_sites is phase state owned elsewhere
    for (size_t i = 0; i < kv_size(transcript->exons); i++) {
        free(kv_A(transcript->exons, i)->ref_bounds);
        free(kv_A(transcript->exons, i));
    }

    kv_destroy(transcript->exons);
    kv_destroy(transcript->introns);
    free(transcript->cds);
    free(transcript);
}

