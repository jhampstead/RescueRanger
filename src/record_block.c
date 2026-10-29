#include "record_block.h"

#include <htslib/regidx.h>

#include "bcftools/gff.h"
#include "logging/log.h"
#include "region.h"
#include "utils.h"

static int find_gene_region_for_variant(bcf1_t *v, bcf_hdr_t *hdr, gff_t *gff, IRegion *region) {
    regidx_t *genes = gff_get(gff, idx_gene);
    regitr_t *itr = regitr_init(genes);
    const char *chr = bcf_hdr_id2name(hdr, v->rid);

    // The REF bases of the variant: the index takes an inclusive end, so the last REF base and not the one after it
    const IRegion query = variant_to_iregion(v);
    if (!regidx_overlap(genes, chr, query.start, iregion_last_base(query), itr)) { regitr_destroy(itr); return 1; };

    IRegion r = {0};
    r.rid = v->rid;
    bool found = false;
    bcf_unpack(v, BCF_UN_STR);
    while (regitr_overlap(itr)) {
        gf_gene_t *gene = regitr_payload(itr, gf_gene_t *);
        if (gene->type != GF_PROTEIN_CODING) continue;

        // The variant belongs to the gene when a base it changes lies in it (an insertion: both bases around it), the same test
        // the later stages use, not when its REF range, with the attachment base, touches the gene
        const IRegion gene_region = iregion_from_bases(v->rid, gene->beg, gene->end); // The GFF end is inclusive
        if (!bcf1_variant_in_region(v, (Range) { gene_region.start, gene_region.end })) continue;

        // The region of the genes the variant is in, half-open like every IRegion
        if (!found || gene_region.start < r.start) r.start = gene_region.start;
        if (!found || gene_region.end > r.end) r.end = gene_region.end;
        found = true;
    }

    regitr_destroy(itr);
    if (!found) return 1; // No protein coding gene was found

    *region = r;

    return 0;
}

enum ReadStatus classify_read(const int bcf_read_result, const int errcode) {
    if (bcf_read_result >= 0) return READ_OK;
    if (bcf_read_result == -1 && errcode == 0) return READ_END; // No next line and no error code: the end of the file

    return READ_FAILED;
}

enum ReadStatus RecordBuffer_fill(htsFile *fp, bcf_hdr_t *hdr, RecordBuffer *buffer) {
    if (buffer->failed) return READ_FAILED;
    if (buffer->unused) return READ_OK; // Record hasn't been processed yet, no need to load in a new one

    const int result = bcf_read(fp, hdr, buffer->v);
    const enum ReadStatus status = classify_read(result, buffer->v->errcode);

    if (status == READ_FAILED) {
        log_error("Could not read the input (result %d, error code %d), stopping after the last record that was read. Line: %s", result, buffer->v->errcode, fp->line.s);
        buffer->failed = true;
    } else if (status == READ_OK) {
        // htslib gives 0 for a record; a result above 0 is not defined, the record is used as it is
        if (result > 0) log_warn("Unexpected result %d when reading line: %s", result, fp->line.s);
        buffer->unused = true;
    }

    return status;
}

// Gives a record the state that depends on its annotation only, so it is set up once and not once per phase.
// A record without any coding transcript (or with more than one ALT, see finish_block) keeps NULL for all of it, which is how
// the rest of the program skips it.
static void init_record_state(Record *r, const RecordAnnotation *annotation) {
    if (kv_size(annotation->tscripts) == 0) return; // A variant in a coding gene can still have no overlapping coding transcripts

    // The transcripts themselves are copied per phase, see RecordBlock_begin_phase
    r->tscripts = malloc(sizeof(Transcripts));
    kv_init(*r->tscripts);
    for (size_t i = 0; i < kv_size(annotation->tscripts); i++) kv_push(Transcript *, *r->tscripts, NULL);

    r->transcript_genes = malloc(sizeof(Genes));
    kv_init(*r->transcript_genes);
    for (size_t i = 0; i < kv_size(annotation->transcript_genes); i++) kv_push(Gene *, *r->transcript_genes, kv_A(annotation->transcript_genes, i));

    r->genes = malloc(sizeof(Genes));
    kv_init(*r->genes);
    for (size_t i = 0; i < kv_size(annotation->genes); i++) kv_push(Gene *, *r->genes, kv_A(annotation->genes, i));

    r->native_scores = malloc(r->genes->n * NUM_PHASES * sizeof(Score *));
    for (int k = 0; k < r->genes->n * NUM_PHASES; k++) r->native_scores[k] = NULL;

    r->diff_scores = malloc(r->genes->n * NUM_PHASES * sizeof(Score *));
    for (int k = 0; k < r->genes->n * NUM_PHASES; k++) r->diff_scores[k] = NULL;

    r->haplotype_scores = malloc(r->genes->n * NUM_PHASES * sizeof(Score *));
    for (int k = 0; k < r->genes->n * NUM_PHASES; k++) r->haplotype_scores[k] = NULL;

    r->isoforms = malloc(r->tscripts->n * NUM_PHASES * sizeof(Isoforms));
    for (int k = 0; k < r->tscripts->n * NUM_PHASES; k++) r->isoforms[k].valid = 0;

    r->flags = calloc(r->tscripts->n * NUM_PHASES, sizeof(uint32_t));
}

// Returns 0 on success, and -1 if the annotation could not be built
static int finish_block(const Session *s, RecordBlock *block) {
    if (!block->is_coding) return 0;

    block->annotation = Annotation_build(&block->records, s);
    if (block->annotation == NULL) return -1;

    for (size_t i = 0; i < kv_size(block->records); i++) {
        Record *r = kv_A(block->records, i);

        // Only the first ALT is read, which need not be the one on the haplotype. Without coding state the record takes no
        // part in any phase, like a record without coding transcripts, and is written unchanged.
        if (r->v->n_allele > 2) {
            log_warn("Skipping on %s:%li, record has %d ALT alleles (split them first, for example with bcftools norm -m-)", bcf_hdr_id2name(s->hdr, r->v->rid), r->v->pos + 1, r->v->n_allele - 1);
            r->skipped = 1;
            continue;
        }

        init_record_state(r, &block->annotation->records[i]);
    }

    return 0;
}

static int emit_block(const Session *s, RecordBlock block, RecordBlock *block_out) {
    if (finish_block(s, &block) != 0) {
        RecordBlock_destroy(block);
        return 2;
    }

    *block_out = block;
    return 0;
}

int RecordBlock_next(const Session *s, RecordBuffer *buffer, RecordBlock *block_out) {
    RecordBlock block = {0};
    kv_init(block.records);

    IRegion batch_region = {0};
    enum ReadStatus read_status;
    while((read_status = RecordBuffer_fill(s->fp_in, s->hdr, buffer)) == READ_OK) {

        // Check if record is coding
        IRegion record_region = {0};
        int is_non_coding = find_gene_region_for_variant(buffer->v, s->hdr, s->gff, &record_region);
        if (is_non_coding) {
            if (block.records.n == 0) {
                kv_push(Record *, block.records, Record_init(buffer->v, NULL));
                buffer->unused = false;
            }

            return emit_block(s, block, block_out);
        }


        // Check if the coding region the record falls in overlaps the current batch
        block.is_coding = 1;
        int overlaps_batch = 1;
        if (block.records.n == 0) {
            batch_region = record_region; // First record in batch has nothing to overlap with
        } else {
            overlaps_batch = iregion_intersects(batch_region, record_region);
            if (!overlaps_batch) {
                return emit_block(s, block, block_out);
            }

            // Extends the batch region to allow for overlapping genes
            if (overlaps_batch & (record_region.end > batch_region.end)) batch_region.end = record_region.end;
        }

        VepCsqs *vep_csqs = parse_vep_transcript_consequences(buffer->v, s->vep_fmt, s->hdr);
        kv_push(Record *, block.records, Record_init(buffer->v, vep_csqs));
        buffer->unused = false;
    }

    // The records that were read before a read failure are a block of their own. The failure is kept in the buffer, so the next
    // call reports it.
    if (block.records.n > 0) {
        return emit_block(s, block, block_out);
    }

    return read_status == READ_FAILED ? 2 : 1;
}

void RecordBlock_destroy(RecordBlock block) {
    for (int i = 0; i < block.records.n; i++) {
        Record_destroy(kv_A(block.records, i));
    }
    kv_destroy(block.records);

    Annotation_destroy(block.annotation);
}

void RecordBlock_begin_phase(RecordBlock *block) {
    for (size_t i = 0; i < kv_size(block->records); i++) {
        Record *r = kv_A(block->records, i);
        if (r->tscripts == NULL) continue;

        const RecordAnnotation *annotation = &block->annotation->records[i];
        for (size_t ti = 0; ti < kv_size(*r->tscripts); ti++) {
            const Transcript *shared = kv_A(annotation->tscripts, ti);
            kv_A(*r->tscripts, ti) = shared == NULL ? NULL : Transcript_phase_copy(shared);
        }
    }
}

void RecordBlock_end_phase(RecordBlock *block) {
    for (size_t i = 0; i < kv_size(block->records); i++) {
        Record *r = kv_A(block->records, i);
        if (r->tscripts == NULL) continue;

        for (size_t ti = 0; ti < kv_size(*r->tscripts); ti++) {
            Transcript *copy = kv_A(*r->tscripts, ti);
            if (copy == NULL) continue;

            Transcript_phase_destroy(copy);
            kv_A(*r->tscripts, ti) = NULL;
        }
    }
}
