#include <cpliceai/range.h>
#include <cpliceai/utils.h>
#include <htslib/kstring.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include <klib/kvec.h>
#include <klib/kstring.h>
#include <klib/klist.h>
#include <klib/ksort.h>

#include <time.h>
#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <cpliceai/predict.h>

#include "bcftools/gff.h"
#include "bcftools/regidx.h"

#include "logging/log.h"
#include "region.h"
#include "record.h"
#include "utils.h"
#include "haplotype.h"
#include "scores.h"

#define AF_FIELD "AF"
#define MNV_RANGE 3

#define BUFFER_SIZE 2024
#define ID 0
#define TYPE 2
#define BEGIN 3
#define END 4
#define DIRECTION 6
#define PHASE 7
#define ATTRIBUTE 8

#define SPLICE_DONOR_VARIANT "splice_donor_variant"
#define SPLICE_ACCEPTOR_VARIANT "splice_acceptor_variant"
#define SPLICE_REGION_VARIANT "splice_region_variant"
#define FRAMESHIFT_VARIANT "frameshift_variant"
#define STOP_GAINED_VARIANT "stop_gained"

#define NUM_PHASES 2

#define MAX_FRAMESHIFT_RADIUS 80
#define MAX_SPLICE_RADIUS 80

#define CAI_DISTANCE 50
#define PRED_WIDTH 101
#define WIDTH 10101
#define VARIANT_POS 5050
#define SPLICEAI_THRESHOLD 0.8
#define SPLICEAI_RESCUE_THRESHOLD 0.2
#define SPLICEAI_CRYPTIC_SPLICE_THRESHOLD 0.5

#define SPLICEAI_DS_AG_INDEX 2  // Index of DS_AG in the INFO/SpliceAI field
#define SPLICEAI_DS_AL_INDEX 3  // Index of DS_AL in the INFO/SpliceAI field
#define SPLICEAI_DS_DG_INDEX 4  // Index of DS_DG in the INFO/SpliceAI field
#define SPLICEAI_DS_DL_INDEX 5  // Index of DS_DL in the INFO/SpliceAI field

#define SPLICE_TRIGGER 0
#define FRAMESHIFT_TRIGGER 1
#define STOP_TRIGGER 2

#define CANONICAL_SPLICE_VALUE 1
#define AF_SPLICE_VALUE 2

#define STOP_RESCUE_VALUE 1
#define STOP_GAIN_VALUE 2

#define PROXIMITY_THRESHOLD 15
#define AF_THRESHOLD 0.01

#define START_CODON "ATG"
#define COMPLEMENTARY_START_CODON "CAT"
#define STOP_CODONS {"TAG", "TGA", "TAA"}
#define COMPLEMENTARY_STOP_CODONS {"ATC", "ACT", "ATT"}

#define SPLICEAI_HAP_TAG "SpliceAI_Haplotype"
#define SPLICEAI_HAP_DESC "##INFO=<ID=SpliceAI_Haplotype,Number=.,Type=String,Description=\"CpliceAIv1.0 haplotype variant annotation from reference to complete haplotype. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"

#define CS_TAG "CS"
#define CS_DESC "##INFO=<ID=CS,Number=.,Type=String,Description=\"Change in splicing between haplotype + variant and haplotype at positions DP_AG, DP_AL, DP_DG, DP_DL. Format: CS_AG|CS_AL|CS_DG|CS_DL|CS_SUM\">"

#define RESCUE_TAG "RESCUE"
#define RESCUE_DESC "##INFO=<ID=RESCUE,Number=.,Type=String,Description=\"Boolean rescue prediction\">"

#define RESCUE_PROB_TAG "RESCUE_PROB"
#define RESCUE_PROB_DESC "##INFO=<ID=RESCUE_PROB,Number=.,Type=String,Description=\"Probability of rescue\">"

#define RESCUE_TYPE_TAG "RESCUE_TYPE"
#define RESCUE_TYPE_DESC "##INFO=<ID=RESCUE_TYPE,Number=.,Type=String,Description=\"Information on type of rescue. Format: RESCUE_TYPE|EXON_NUMBER|PROBABILITY|AFFECTED_BASES\">"

#define PTV_TAG "PTV"
#define PTV_DESC "##INFO=<ID=PTV,Number=.,Type=String,Description=\"Newly predicted PTV consequence\">"

#define PTV_PROB_TAG "PTV_PROB"
#define PTV_PROB_DESC "##INFO=<ID=PTV_PROB,Number=.,Type=String,Description=\"Probability of induced PTV\">"

#define PTV_TYPE_TAG "PTV_TYPE"
#define PTV_TYPE_DESC "##INFO=<ID=PTV_TYPE,Number=.,Type=String,Description=\"Type of PTV\">"

#define PREDICTION_FLAGS_TAG "PREDICTION_FLAGS"
#define PREDICTION_FLAGS_DESC "##INFO=<ID=PREDICTION_FLAGS,Number=.,Type=String,Description=\"Additional flags for the context in which the prediction happened.\">"


#define INFO_TAGS { SPLICEAI_TAG, SPLICEAI_HAP_TAG, CS_TAG, RESCUE_TAG, RESCUE_PROB_TAG, RESCUE_TYPE_TAG, PTV_TAG, PTV_PROB_TAG, PTV_TYPE_TAG, PREDICTION_FLAGS_TAG }
#define INFO_TAG_DESCS { SPLICEAI_DESC, SPLICEAI_HAP_DESC, CS_DESC, RESCUE_DESC, RESCUE_PROB_DESC, RESCUE_TYPE_DESC, PTV_DESC, PTV_PROB_DESC, PTV_TYPE_DESC, PREDICTION_FLAGS_DESC }

typedef struct {
    char *models;
    char *gff;
    char *reference;
    char *in;
    char *out;
    char *scores;
    char *use_af;
    char *use_gt;
    float af_threshold;
    int include_unphased;
    // int max_unphased;
} Arguments;


// enum RescueType {
//     COMPLETE,
//     POSSIBLE_REINITIATION,
//     SKIPPED_EXON,
//     PARTIAL_EXON_SKIPPING,
//     PARTIAL_INTRON_RETENTION,
//     EXON_SKIPPING,
//     INTRON_RETENTION,
//     NONE
// };

enum SpliceType {
    ACCEPTOR = 1,
    DONOR = 2
};

typedef struct SpliceProbability {
    Record *record;
    enum SpliceType type;
    float ref_prior_probability;
    float prior_probability;
    float posterior_probability;
    uint32_t position;
    uint32_t canonical_position;
} SpliceProbability;

typedef kvec_t(SpliceProbability) SpliceProbabilities;

typedef struct RecordBlock {
    int is_coding;
    Records *records;
} RecordBlock;

typedef struct RecordBlock2 {
    int is_coding;
    Records records;
} RecordBlock2;

typedef kvec_t(Records *) PhasedRecords;

typedef struct RecordBuffer {
    bcf1_t *v;
    bool unused;
} RecordBuffer;

enum ShiftType {
    INDEL = 1,
    SPLICE = 2
};

typedef struct IsoformFrameshift {
    int coding;
    int partial;
    enum ShiftType shift_type;
    uint32_t seq_pos;
    int shift;
    int modified_bases;
    Record *r;
} IsoformVariant;

#define empty_free(x) 
KLIST_INIT(RecordContext, IsoformVariant, empty_free);


#define frameshift_lt(a, b) ((a).seq_pos < (b).seq_pos)
KSORT_INIT(frameshift_sort, IsoformVariant, frameshift_lt)

#define frameshift_gt(a, b) ((a).seq_pos > (b).seq_pos)
KSORT_INIT(frameshift_reverse_sort, IsoformVariant, frameshift_gt)

int find_gf_tscripts_by_iregion(IRegion iregion, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts) {
    regidx_t *gff_tscripts = (regidx_t *) gff_get(gff, idx_tscript);
    regitr_t *itr = regitr_init(gff_tscripts);

    const char *chr = bcf_hdr_id2name(hdr, iregion.rid);
    if (!regidx_overlap(gff_tscripts, chr, iregion.start, iregion.end, itr)) { regitr_destroy(itr); return -1; }

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

void predict_for_seq(Model *models, const char *seq, int seq_len, int strand, float *predictions[], int *num_predictions) {
    float *ohe; 
    int enc_len = one_hot_encode(seq, seq_len, &ohe);

    if (strand == STRAND_REV) reverse_encoding(ohe, enc_len);

    // log_info("Enc len: %i", enc_len);
    if (predict(models, enc_len, 1, &ohe, num_predictions, predictions)) {
        log_info("Predict failed.");
        return;
    }
    free(ohe);

    if (strand == STRAND_REV) reverse_prediction(*predictions, *num_predictions, 3);
}

// void get_best_insertion_score(const float *insert_start, const int length, float *best) {
//     for (int i = 0; i < length; i+=3) {
//         for (int j = 0; j < 3; j++) {
//             if (insert_start[i + j] > best[j]) best[j] = insert_start[i+j];
//         }
//     }
// }

void unpack_variants(Records *vs) {
    for (int i = 0; i < vs->n; i++) {
        bcf_unpack(kv_A(*vs, i)->v, BCF_UN_ALL);
    }
}

PhasedRecords *split_records_by_phase(Records *vs, bcf_hdr_t *hdr) {
    PhasedRecords *pvs = malloc(sizeof(PhasedRecords)); kv_init(*pvs); 
    for (int i = 0; i < NUM_PHASES; i++) {
        Records *tmp = malloc(sizeof(Records)); kv_init(*tmp);
        kv_push(Records *, *pvs, tmp);
    }

    for (int i = 0; i < vs->n; i++) {
        Record *r = kv_A(*vs, i);

        int *gt = NULL; int ngts = 0;
        bcf_get_genotypes(hdr, r->v, &gt, &ngts);

        if (gt[0] == gt[1]) { // Homozygous variant, equal opportunity phaser
            kv_push(Record *, *(pvs->a[0]), r);
            kv_push(Record *, *(pvs->a[1]), r);
            free(gt);
            continue;
        }

        // TODO: Add a CLI option to set this behaviour
        if (!bcf_gt_is_phased(gt[0]) && !bcf_gt_is_phased(gt[1])) { // No phasing for heterozygous variants, including on both for now
            r->unphased = 1;
            kv_push(Record *, *(pvs->a[0]), r);
            kv_push(Record *, *(pvs->a[1]), r);
            log_debug("Adding %s:%li (heterozygous, no phase)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos); 
            free(gt);
            continue; 
        }

        //Something is phased, check what. Enforced split of multiallelic sites makes this easy
        if (gt[0] > gt[1]) kv_push(Record *, *(pvs->a[0]), r); // gt[0]=4, gt[1]=3, 1|0
        else kv_push(Record *, *(pvs->a[1]), r); // gt[0]=2, gt[1]=5, 0|1
        free(gt);
    }

    return pvs;
}

void filter_overlapping_variants(Records *vs, const bcf_hdr_t *hdr) {
    Range last_variant_range = {0};
    int last_variant_type = 0;

    Range current_variant_range = {0};
    int current_variant_type = {0};

    // int *should_delete = calloc(vs->n, sizeof(int));
    int should_delete[vs->n]; for (int i = 0; i < vs->n; i++) should_delete[i] = 0;

    for (int i = 0; i < vs->n; i++) {
        bcf1_t *v = kv_A(*vs, i)->v;

        if (current_variant_range.end > last_variant_range.end) {
            last_variant_range = current_variant_range;
            last_variant_type = current_variant_type;
        }

        current_variant_range = variant_to_range(v, 0);
        Range intersection;
        int overlaps = intersect_ranges(current_variant_range, last_variant_range, &intersection);

        if (range_len(variant_to_range(v, 1)) > 1) { // Insertion
            current_variant_type = 1;
            current_variant_range.start++;

            if (last_variant_type == 0) { // SNV
                continue; // Cannot overlap a SNV
            } else if (last_variant_type == -1) { // Deletion
                // Insertions are not skipped as long as they don't fall inside of another deletion
                if (current_variant_range.start <= last_variant_range.start | current_variant_range.end >= last_variant_range.end) continue;
            } else { // Insertion
                // Insertions are not skipped as long as their start does not overlap another insertions
                if (last_variant_range.start != current_variant_range.start | last_variant_range.end != current_variant_range.end) continue;
            }
        } else if (range_len(current_variant_range) > 1) { // Deletion
            current_variant_type = -1;
            current_variant_range.start++;
            if (current_variant_range.start >= last_variant_range.end) continue;
        } else { // SNV
            current_variant_type = 0;
            if (!overlaps || range_len(intersection) == 0) continue;
        }

        log_warn("Skipping on %s:%li, overlaps previous variant", bcf_hdr_id2name(hdr, v->rid), v->pos + 1);
        should_delete[i] = 1;
    }

    int num_variants = vs->n;
    for (int i = 0, vs_i = 0; i < num_variants; i++, vs_i++) {
        if (!should_delete[i]) continue;

        int current_num_variants = vs->n;
        if (vs_i == current_num_variants-1) {
            kv_pop(*vs);
            vs_i--;
            continue;
        }

        for (int j = vs_i; j < vs->n-1; j++) {
            vs->a[j] = vs->a[j+1];
        }
        vs->n--; vs_i--;
    }
}

// int prediction_offset(bcf1_t *v, Records *rs) {
//     Range pred_range = {v->pos - (CAI_DISTANCE), v->pos + CAI_DISTANCE+1};
//
//     int offset = 0;
//     for (int i = 0; i < rs->n; i++) {
//         Record *r = kv_A(*rs, i);
//
//         Range ref;
//         if (intersect_ranges((Range) {r->v->pos, r->v->pos + r->v->rlen}, pred_range, &ref) == 0) continue;
//         int ref_len = range_len(ref);
//         if (ref_len == 0) continue;
//
//         Range alt = variant_to_range(r->v, 1);
//         Range intersect_alt;
//         int alt_len = 0;
//         if (intersect_ranges(alt, pred_range, &intersect_alt) > 0) {
//             alt_len = range_len(alt);
//         }
//         offset += ref_len - alt_len;
//     }
//
//     return offset;
// }

void add_phasing_for_prediction_records(Records *rs, Range prediction_range) {
    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        Range ref;
        if (!intersect_ranges((Range) {r->v->pos, r->v->pos + r->v->rlen}, prediction_range, &ref) || range_len(ref) == 0) continue;

        Range alt;
        if (!intersect_ranges((Range) {r->v->pos, r->v->pos + strlen(r->v->d.allele[1])}, prediction_range, &alt) || range_len(alt) == 0) continue;

        for (int j = 0; j < r->transcript_genes->n; j++) {
            if (r->skipped) r->flags[j] |= SKIPPED_VARIANT_IN_PREDICTION;
            if (r->unphased) r->flags[j] |= UNPHASED_VARIANT_IN_PREDICTION;
        }
    }

}

int find_gene_region_for_variant(bcf1_t *v, bcf_hdr_t *hdr, gff_t *gff, IRegion *region) {
    regidx_t *genes = gff_get(gff, idx_gene);
    regitr_t *itr = regitr_init(genes);
    const char *chr = bcf_hdr_id2name(hdr, v->rid);

    if (!regidx_overlap(genes, chr, v->pos, v->pos + v->rlen, itr)) { regitr_destroy(itr); return 1; };

    IRegion r = {0};
    r.rid = v->rid;
    bcf_unpack(v, BCF_UN_STR);
    while (regitr_overlap(itr)) {
        gf_gene_t *gene = regitr_payload(itr, gf_gene_t *);
        if (gene->type != GF_PROTEIN_CODING) continue;

        Range gene_range = (Range) { gene->beg, gene->end };
        Range variant_range = variant_to_range(v, 0);
        Range intersection;
        if (intersect_ranges(gene_range, variant_range, &intersection) == 0) continue;
        if (range_len(intersection) == 0) continue;

        if (r.start == 0) r.start = gene->beg;
        if (itr->end > r.end) r.end = gene->end;
    }

    regitr_destroy(itr);
    if (r.start == 0 && r.end == 0) return 1; // No protein coding gene was found

    *region = r;

    return 0;
}

// float *match_predictions_to_reference(const float *predictions, const bcf1_t *v, const Records *rs, const size_t num_reference_predictions) {
//     float *corrected_predictions = malloc(num_reference_predictions * sizeof(float));
//     Range pred_range = {v->pos - (CAI_DISTANCE), v->pos + CAI_DISTANCE+1}; // +1 for open ended coordinates
//
//     int ref_idx = 0, hap_idx = 0;
//     for (int l = 0; l < rs->n; l++) {
//         Record *r = kv_A(*rs, l);
//
//         Range ref;
//         if (intersect_ranges(variant_to_range(r->v, 0), pred_range, &ref) == 0) continue;
//         shift_range(&ref, -pred_range.start);
//         int ref_len = range_len(ref);
//         if (ref_len == 0) continue;
//
//         Range alt = variant_to_range(r->v, 1);
//         Range intersected_alt;
//         int alt_segment = 0;
//         if (intersect_ranges(alt, pred_range, &intersected_alt) > 0) {
//             shift_range(&alt, -pred_range.start);
//             alt_segment = range_len(alt);
//         }
//
//         int skipped_alt_bases = (intersected_alt.start - alt.start);
//         if (alt_segment > 0) alt_segment -= skipped_alt_bases;
//         Range ref_segment = scale_range(ref, NUM_SCORES);
//
//         // Copy segments between variants/boundaries
//         int segment_len = (ref_segment.start - ref_idx);
//         if (segment_len > 0) memcpy(corrected_predictions + ref_idx, predictions + hap_idx, segment_len * sizeof(float));
//
//         int variant_offset = haplotype_offset(rs, (Range) { pred_range.start, pred_range.start + ref.start }) * NUM_SCORES;
//         int offset = ref_len - alt_segment; // Naive offset in end position vs ref and alt
//         float *corrected_start = corrected_predictions + ref_segment.start;
//         if (offset == 0) {
//             memcpy(corrected_start, predictions + ref_segment.start - variant_offset, range_len(ref_segment) * sizeof(float));
//             // memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, range_len(ref_segment) * sizeof(float));
//         } else if (offset < 0) {
//             float best_scores[NUM_SCORES] = {0};
//             get_best_insertion_score(predictions + ref_segment.start - variant_offset, alt_segment, best_scores);
//             memcpy(corrected_start, best_scores, NUM_SCORES * sizeof(float));
//             // memcpy(corrected_predictions + ref_segment.start, best_scores, NUM_SCORES * sizeof(float));
//         } else { // offset > 0
//             // Only set if it would not already have been set
//             // if (ref_segment.start > ref_idx) 
//             memcpy(corrected_start, predictions + ref_segment.start - variant_offset, NUM_SCORES * sizeof(float));
//             memset(corrected_start + NUM_SCORES, 0, offset * NUM_SCORES + sizeof(float));
//             // memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, NUM_SCORES * sizeof(float));
//             // memset(corrected_predictions + (ref_segment.start + NUM_SCORES), 0, offset * NUM_SCORES * sizeof(float));
//         }
//
//
//         // hap_idx += segment_len + alt_len*NUM_SCORES;
//         if (ref_segment.end > ref_idx) ref_idx = ref_segment.end;
//     }
//
//     int variant_offset = haplotype_offset(rs, pred_range) * NUM_SCORES;
//     memcpy(corrected_predictions + ref_idx, predictions + ref_idx - variant_offset, (num_reference_predictions - ref_idx) * sizeof(float));
//
//     return corrected_predictions;
// }

static int get_range_in_isoform(const bcf1_t *v, const NonCodingRegions *ncrs, const NonCodingRegions *hap_ncrs, const Records *rs, const int strand, const Range bounds, int64_t *out_position, int allele_lengths[2]) {
    const Range ref = variant_to_range(v, 0);
    const Range alt = variant_to_range(v, 1);
    const int is_insertion = range_len(alt) - range_len(ref) > 0;
    const int is_deletion = range_len(ref) - range_len(alt) > 0;
    const int has_proceeding_reference = range_len(alt) != range_len(ref);

    const int64_t variant_start = ref.start;
    const int64_t variant_hap_start = variant_start - haplotype_offset(rs, (Range) { bounds.start, ref.start });
    const int64_t variant_end = variant_start + range_len(ref);
    int skipped_bases = 0;
    if (strand == STRAND_FWD) {
        for (int i = 0; i < hap_ncrs->n; i++) {
            NonCodingRegion ncr = kv_A(*ncrs, i);
            NonCodingRegion hap_ncr = kv_A(*hap_ncrs, i);

            // When a noncoding region is ablated (only happens in case of intron retention due to splice site loss)
            if (!ncr.is_used) continue;

            // Convert to actual start of intron, instead of closed end of exon (not applicable if it's a UTR)
            if (!ncr.is_utr) {
                ncr.coords.start++;
                hap_ncr.coords.start++;
            }

            // Range is open-ended
            // If ther variant is an insertion, this is generally included as it is before the splice site.
            // But in case of a UTRs it's not bounded by a splice site so this should not happen.

            if (variant_start + (is_insertion * !ncr.is_utr) >= ncr.coords.end) {
                if (ncr.is_utr) continue;
                skipped_bases += range_len(hap_ncr.coords);
                continue;
            }

            if (variant_end + (is_insertion * !!ncr.is_utr) <= ncr.coords.start) break;

            if (is_deletion && variant_start + 1 < ncr.coords.start) {
                Range noncoding_variant;

                // If this if hits, it's a UTR bounding variant
                if (intersect_ranges(ref, ncr.coords, &noncoding_variant) == 0) {
                    log_error("Attempted an overlap of a variant to an intron boundary without actually overlapping. Big bad.");
                    exit(1);
                }

                *out_position = variant_hap_start - skipped_bases;

                allele_lengths[0] = range_len(ref) - range_len(noncoding_variant);
                allele_lengths[1] = 1;
                return 0;

            }

            if (is_deletion && variant_end > ncr.coords.end) {
                Range noncoding_variant;

                // If this if hits, it's a UTR bounding variant
                if (intersect_ranges(ref, ncr.coords, &noncoding_variant) == 0) {
                    log_error("Attempted an overlap of a variant to an intron boundary without actually overlapping. Big bad.");
                    exit(1);
                }

                *out_position = variant_hap_start - (skipped_bases + range_len(hap_ncr.coords));

                allele_lengths[0] = range_len(ref) - range_len(noncoding_variant);
                allele_lengths[1] = 0;
                return 0;
            }

            *out_position = variant_hap_start - skipped_bases;
            return 1; // Variant falls within a non-coding region
        }
    } else if (strand == STRAND_REV) {
        for (int i = 0; i < ncrs->n; i++) {
            NonCodingRegion ncr = kv_A(*ncrs, i);
            NonCodingRegion hap_ncr = kv_A(*hap_ncrs, i);

            if (!ncr.is_used) continue;

            // Convert to actual start of intron, instead of closed end of exon
            if (!ncr.is_utr) { 
                ncr.coords.start--;
                hap_ncr.coords.start--;
            }

            if (variant_end + (is_insertion * !!ncr.is_utr) <= ncr.coords.end + 1) {
                if (ncr.is_utr) continue;
                skipped_bases += range_len(hap_ncr.coords);
                continue;
            }

            if (variant_start + has_proceeding_reference > ncr.coords.start) break;

            if (is_deletion && variant_end > ncr.coords.start + 1) {
                Range noncoding_variant;
                if (intersect_ranges(ref, ncr.coords, &noncoding_variant) == 0) {
                    log_error("Attempted an overlap of a variant to an intron boundary without actually overlapping. Big bad.");
                    exit(1);
                }

                *out_position = variant_hap_start + skipped_bases;
                allele_lengths[0] = range_len(ref) - range_len(noncoding_variant);
                allele_lengths[1] = 0;
            }

            if (is_deletion && variant_start < ncr.coords.end) {
                Range noncoding_variant;
                if (intersect_ranges(ref, ncr.coords, &noncoding_variant) == 0) {
                    log_error("Attempted an overlap of a variant to an intron boundary without actually overlapping. Big bad.");
                    exit(1);
                }

                *out_position = variant_hap_start + skipped_bases + range_len(hap_ncr.coords);
                allele_lengths[0] = range_len(ref) - range_len(noncoding_variant);
                allele_lengths[1] = 1;
            }

            *out_position = variant_start;
            return 1; // Variant falls within a non-coding region
        }
    }

    *out_position = strand == STRAND_FWD ? variant_hap_start - skipped_bases : variant_hap_start + skipped_bases;
    allele_lengths[0] = range_len(ref);
    allele_lengths[1] = range_len(alt);

    return 0; // Coding variant
}

static int get_position_in_isoform(const int64_t position, const NonCodingRegions *ncrs, const Records *rs, const int strand, const Range bounds, int64_t *out_position) {
    int64_t isoform_position = position; //  - haplotype_offset(rs, bounds);
    int skipped_bases = 0;
    if (strand == STRAND_FWD) {
        // isoform_position -= haplotype_offset(rs, (Range) { bounds.start, position-1 });
        for (int i = 0; i < ncrs->n; i++) {
            NonCodingRegion ncr = kv_A(*ncrs, i);

            if (position < ncr.coords.start) break;

            if (position >= ncr.coords.end) {
                if (ncr.is_utr) continue;
                skipped_bases += range_len(ncr.coords) - 1;
                continue;
            }

            isoform_position -= skipped_bases;
            *out_position = isoform_position;
            return 1; // Variant falls within a non-coding region
        }
        isoform_position -= skipped_bases;
    } else if (strand == STRAND_REV) {
        // isoform_position -= haplotype_offset(rs, (Range) { bounds.start, position-1 });
        for (int i = 0; i < ncrs->n; i++) {
            NonCodingRegion ncr = kv_A(*ncrs, i);

            if (position > ncr.coords.start) break; // Approaches from other end (3' <-- 5')

            if (position <= ncr.coords.end) {
                if (ncr.is_utr) continue;
                skipped_bases += range_len(ncr.coords) - 1;
                continue;
            }

            isoform_position += skipped_bases;
            *out_position = isoform_position;
            return 1; // Variant falls within a non-coding region
        }
        isoform_position += skipped_bases;
    }

    *out_position = isoform_position;

    return 0;
}

static enum UTRType is_utr_variant(const Record *r, const NonCodingRegions *ncrs, const Records *rs, const int strand, const Range bounds) {
    const int64_t position = r->v->pos;
    const int64_t isoform_position = position - haplotype_offset(rs, (Range) { bounds.start, position-1 });
    const Range variant = { isoform_position, isoform_position + r->v->rlen };

    for (int i = 0; i < ncrs->n; i++) {
        NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_utr) continue;

        if (ranges_overlap(ncr.coords, variant)) return ncr.is_utr;
    }

    return 0;
}

void calculate_rescue_probabilities(Isoforms *isoforms_struct) {
    isoforms_struct->rescue_probability = 0.0f;
    RescueProbabilities *rescue_probs = malloc(sizeof(RescueProbabilities));

    rescue_probs->n_isoforms = malloc(isoforms_struct->m_rescue_types * sizeof(int32_t));
    rescue_probs->rescues = malloc(isoforms_struct->m_rescue_types * sizeof(uint8_t));
    rescue_probs->rescue_types = malloc(isoforms_struct->m_rescue_types * sizeof(enum RescueType));
    rescue_probs->rescue_probabilities = malloc(isoforms_struct->m_rescue_types * sizeof(float));
    rescue_probs->affected_bases = malloc(isoforms_struct->m_rescue_types * sizeof(int32_t));

    for (int i = 0; i < isoforms_struct->m_rescue_types; i++) {
        rescue_probs->n_isoforms[i] = 0;
        rescue_probs->rescues[i] = 0;
        rescue_probs->rescue_probabilities[i] = 0.0f;
        rescue_probs->rescue_types[i] = NO_RESCUE; // initialise as NO_RESCUE
        rescue_probs->affected_bases[i] = 0;
    }

    for (int i = 0; i < isoforms_struct->n_isoforms; i++) {
        Isoform *isoform = isoforms_struct->isoforms[i];

        // If rescue is predicted on the isoform, add isoform probability to the rescue probability [INFO/RESCUE_PROB]
        if (isoform->rescue == 1) {
            isoforms_struct->rescue_probability += isoform->posterior_probability;
        }

        // Get probabilities by rescue type [INFO/RESCUE_TYPE]. If the rescue type is already recorded, add the isoform and sum the probability
        int found = 0;
        for (int j = 0; j < isoforms_struct->n_rescue_types; j++) {
            if ((rescue_probs->rescue_types[j] != isoform->rescue_type) || (rescue_probs->affected_bases[j] != isoform->affected_bases)) continue;

            rescue_probs->n_isoforms[j] += 1;
            rescue_probs->rescue_probabilities[j] += isoform->posterior_probability;
            if (isoform->affected_bases > rescue_probs->affected_bases[j]) rescue_probs->affected_bases[j] = isoform->affected_bases;
            found = 1;
            break;
        }

        // If rescue type isn't found, add a new rescue type, dynamically reallocating memory if necessary
        if (found) continue;

        if (isoforms_struct->n_rescue_types >= isoforms_struct->m_rescue_types) {
            isoforms_struct->m_rescue_types *= 2; // double the size
            rescue_probs->n_isoforms = realloc(rescue_probs->n_isoforms, isoforms_struct->m_rescue_types * sizeof(int32_t));
            rescue_probs->rescues = realloc(rescue_probs->rescues, isoforms_struct->m_rescue_types * sizeof(uint8_t));
            rescue_probs->rescue_probabilities = realloc(rescue_probs->rescue_probabilities, isoforms_struct->m_rescue_types * sizeof(float));
            rescue_probs->rescue_types = realloc(rescue_probs->rescue_types, isoforms_struct->m_rescue_types * sizeof(enum RescueType));
            rescue_probs->affected_bases = realloc(rescue_probs->affected_bases, isoforms_struct->m_rescue_types * sizeof(int32_t));

            for (int k = isoforms_struct->n_rescue_types; k < isoforms_struct->m_rescue_types; k++) {
                rescue_probs->n_isoforms[k] = 0;
                rescue_probs->rescues[k] = 0;
                rescue_probs->rescue_probabilities[k] = 0.0f;
                rescue_probs->rescue_types[k] = NO_RESCUE; // initialise with NULL
                rescue_probs->affected_bases[k] = 0;
            }
        }

        rescue_probs->n_isoforms[isoforms_struct->n_rescue_types] = 1;
        rescue_probs->rescues[isoforms_struct->n_rescue_types] = isoform->rescue;
        rescue_probs->rescue_types[isoforms_struct->n_rescue_types] = isoform->rescue_type;
        rescue_probs->rescue_probabilities[isoforms_struct->n_rescue_types] = isoform->posterior_probability;
        rescue_probs->affected_bases[isoforms_struct->n_rescue_types] = isoform->affected_bases;
        isoforms_struct->n_rescue_types++;
    }

    isoforms_struct->rescue = (isoforms_struct->rescue_probability > 0.5) ? 1 : 0; // set rescue value on the basis of rescue_probability
    isoforms_struct->rescue_type_probabilities = rescue_probs; // assign probabilities per rescue type back to struct

}

void calculate_ptv_probabilities(Isoforms *isoforms_struct) {
    float ptv_probability = 0.0;
    float stop_gained_probability = 0.0;
    float frameshift_variant_probability = 0.0;
    float long_inframe_indel_probability = 0.0;

    for (int i = 0; i < isoforms_struct->n_isoforms; i++) {
        Isoform *isoform = isoforms_struct->isoforms[i];
        if (isoform->ptv == 0) continue;

        ptv_probability += isoform->posterior_probability;

        if (isoform->ptv_type == STOP_GAINED) {
            stop_gained_probability += isoform->posterior_probability;
        } else if (isoform->ptv_type == FRAMESHIFT) {
            frameshift_variant_probability += isoform->posterior_probability;
        } else if (isoform->ptv_type == LONG_INFRAME_INDEL) {
            long_inframe_indel_probability += isoform->posterior_probability;
        }
    }

    isoforms_struct->ptv_probability = ptv_probability;
    isoforms_struct->ptv = (isoforms_struct->ptv_probability > 0.5) ? 1 : 0;

    if (stop_gained_probability > frameshift_variant_probability) {
        isoforms_struct->ptv_type = STOP_GAINED;
    } else if (frameshift_variant_probability > stop_gained_probability) {
        isoforms_struct->ptv_type = FRAMESHIFT;
    } else if (long_inframe_indel_probability > stop_gained_probability && long_inframe_indel_probability > stop_gained_probability) {
        isoforms_struct->ptv_type = LONG_INFRAME_INDEL;
    } else {
        isoforms_struct->ptv_type = NO_PTV; // tie or no PTV
    }
}

kstring_t generate_rescue_type_tag(Isoforms *isoforms_struct) { // INFO/RESCUE_TYPE rescue_type|probability|affected_bases
    // Calculate total length for allocation
    kstring_t tag = {0};

    for (int i = 0; i < isoforms_struct->n_rescue_types; i++) {
        if (tag.l > 0) kputc(':', &tag);

        kputw(isoforms_struct->rescue_type_probabilities->rescues[i], &tag);
        kputc('|', &tag);

        enum RescueType rescue_type = isoforms_struct->rescue_type_probabilities->rescue_types[i];
        switch(rescue_type) {
            case COMPLETE_RESCUE:
                kputs("complete_rescue", &tag);
                break;
            case POSSIBLE_REINITIATION_SITE:
                kputs("possible_reinitiation", &tag);
                break;
            case FRAME_RESTORING_INDEL:
                kputs("frame_restoring_indel", &tag);
                break;
            case PARTIAL_EXON_SKIPPING:
                kputs("partial_exon_skipping", &tag);
                break;
            case PARTIAL_INTRON_RETENTION:
                kputs("partial_intron_retention", &tag);
                break;
            case EXON_SKIPPING:
                kputs("exon_skipping", &tag);
                break;
            case INTRON_RETENTION:
                kputs("intron_retention", &tag);
                break;
            case NO_RESCUE:
                kputs("no_rescue", &tag);
                break;
        }

        kputc('|', &tag);
        char prob_s[7];
        sprintf(prob_s, "%.4f", isoforms_struct->rescue_type_probabilities->rescue_probabilities[i]);
        kputsn(prob_s, 6, &tag);
        kputc('|', &tag);
        kputw(isoforms_struct->rescue_type_probabilities->affected_bases[i], &tag);

    }

    return tag;
}

void free_isoforms(Isoforms *isoforms_struct) {
    for (int i = 0; i < isoforms_struct->n_isoforms; i++) {
        free(isoforms_struct->isoforms[i]);
    }
    free(isoforms_struct->isoforms);

    free(isoforms_struct->rescue_type_probabilities->rescue_types);
    free(isoforms_struct->rescue_type_probabilities->n_isoforms);
    free(isoforms_struct->rescue_type_probabilities->rescue_probabilities);
    free(isoforms_struct->rescue_type_probabilities->affected_bases);
    free(isoforms_struct->rescue_type_probabilities);
}

#define TAG_SIZE 10
static int add_rescueranger_tags(bcf_hdr_t *hdr) {
    const char *tags[TAG_SIZE] = INFO_TAGS;
    const char *descriptions[TAG_SIZE] = INFO_TAG_DESCS;

    for (int i = 0; i < TAG_SIZE; i++) {
        if (bcf_hdr_append(hdr, descriptions[i]) != 0) {
            log_error("Failed to append description for tag %s to vcf header.", tags[i]);
            return -1;
        }
    }

    return 0;
}
#undef TAG_SIZE

static void filter_to_most_likely_prediction(size_t num_probabilities, float *probabilities) {
    for (int i = 0; i < num_probabilities; i += NUM_SCORES) {
        float max = 0;
        int max_idx = 0;

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if (probabilities[pred_idx] <= max) continue;

            max = probabilities[pred_idx];
            max_idx = pred_idx;
        }

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if (pred_idx == max_idx) continue;

            probabilities[pred_idx] = 0;
        }
    }
}

static int transform_to_posterior_probabilities(const size_t num_probabilities, const float *probabilties, float *posterior_probabilities) {
    float sum = 0;
    int num_values = 0;

    for (int i = 0; i < num_probabilities; i++) {
        if (probabilties[i] == 0) continue;

        num_values++;
        sum += probabilties[i];
    }
    memset(posterior_probabilities, 0.0, num_probabilities * sizeof(float));

    for (int i = 0; i < num_probabilities; i++) {
        if (probabilties[i] == 0) continue;

        posterior_probabilities[i] = probabilties[i] / sum;
    }

    return num_values;
}

static int transform_to_posterior_site_probabilities(const int target, const size_t num_probabilities, const float *probabilities, float *posterior_probabilities) {
    int i = 0;
    float sum = 0;
    int num_values = 0;

    for (i = target; i < num_probabilities; i+=NUM_SCORES) {
        if (probabilities[i] == 0) continue;

        num_values++;
        sum += probabilities[i];
    }

    memset(posterior_probabilities, 0.0, num_probabilities * sizeof(float));

    if (num_values == 0) return 0;

    for (i = target; i < num_probabilities; i+=NUM_SCORES) {
        if (probabilities[i] == 0) continue;

        posterior_probabilities[i] = probabilities[i] / sum;
    }

    return num_values;
}

static int set_transcript_exons(Transcript *transcript, const char *chr, const gff_t *gff, const int *constitutive) {
    regidx_t *exons = gff_get((gff_t *)gff, idx_exon);
    regitr_t *itr = regitr_init(exons);

    if (regidx_overlap(exons, chr, transcript->start, transcript->end, itr) == 0) {
        log_warn("No overlapping exons found for %s:%i-%i", chr, transcript->start, transcript->end);
        regitr_destroy(itr);
        return 0;
    }

    transcript->nexons = 0;
    transcript->mexons = transcript->mcds ? transcript->mcds : 4;
    transcript->exons = malloc(transcript->mexons * sizeof(Exon *));

    int ci = 0;
    while (regitr_overlap(itr)) {
        gf_exon_t *gf_e = regitr_payload(itr, gf_exon_t *);

        // Exons for other transcripts and genes can overlap, but we don't want them!
        if (gf_e->tr->id != transcript->id) continue;

        Range *cds = NULL;
        if (ci < transcript->ncds) cds = &transcript->cds[ci];

        Range *exon = malloc(sizeof(Range));
        *exon = (Range) {gf_e->beg, gf_e->end};

        // Resize as needed
        if (transcript->nexons >= transcript->mexons) {
            transcript->mexons *= 2;
            transcript->exons = realloc(transcript->exons, transcript->mexons * sizeof(Exon *));
        }


        transcript->exons[transcript->nexons] = malloc(sizeof(Exon));
        // Does this exon have a CDS?
        Range cds_intersection = {0};
        if (cds != NULL) intersect_ranges(*exon, *cds, &cds_intersection);
        if (range_len(cds_intersection) == 0) {
            *transcript->exons[transcript->nexons++] = (Exon) { NULL, exon, NULL, NULL, false, false};
            continue;
        }

        *transcript->exons[transcript->nexons++] = (Exon) { cds, exon, NULL, NULL, constitutive[ci], range_len(*cds) % 3 == 0};
        ci++;
    }


    regitr_destroy(itr);
    return 1;
}

static void shift_exons_to_haplotype_sequence(const Transcript *t, const Range gene, const Records *rs) {
    for (int i = 0; i < t->nexons; i++) {
        Exon *exon = t->exons[i];

        // Shift ranges to by start of gene to make them relative to gene start
        shift_range(exon->hap_bounds, -gene.start);
        exon->hap_bounds->start -= haplotype_offset(rs, (Range) { gene.start, exon->hap_bounds->start });
        exon->hap_bounds->end -= haplotype_offset(rs, (Range) { gene.start, exon->hap_bounds->end });

        if (exon->cds == NULL) continue; // Can't shift the CDS if it doesn't exist
        shift_range(exon->cds, -gene.start);
        exon->cds->start -= haplotype_offset(rs, (Range) { gene.start, exon->cds->start });
        exon->cds->end -= haplotype_offset(rs, (Range) { gene.start, exon->cds->end });
    }
}

int trigger_cryptic_splice(bcf1_t *rec, bcf_hdr_t *hdr) {
    bcf_info_t *info = bcf_get_info(hdr, rec, "SpliceAI");
    if (info == NULL) {
        log_info("SpliceAI tag not found in record");
        return 0;  // Cryptic splice trigger is never activated if SpliceAI scores don't exist in input VCF record
    }

    char *spliceai_string = (char*) info->vptr;
    int len = info->len;

    char *spliceai_copy = strndup(spliceai_string, len); // Make a copy of the SpliceAI string for safe tokenising
    char *token;
    int i = 0;
    double spliceai_scores[4] = {0.0};  // Store DS_AG, DS_AL, DS_DG, DS_DL

    token = strtok(spliceai_copy, "|");
    while (token != NULL) {
        if (i == SPLICEAI_DS_AG_INDEX) {
            spliceai_scores[0] = atof(token);
        } else if (i == SPLICEAI_DS_AL_INDEX) {
            spliceai_scores[1] = atof(token);
        } else if (i == SPLICEAI_DS_DG_INDEX) {
            spliceai_scores[2] = atof(token);
        } else if (i == SPLICEAI_DS_DL_INDEX) {
            spliceai_scores[3] = atof(token);
        }

        token = strtok(NULL, "|");
        i++;
    }

    free(spliceai_copy);

    for (int i = 0; i < 4; i++) {
        if (spliceai_scores[i] >= SPLICEAI_CRYPTIC_SPLICE_THRESHOLD) {
            // log_info("Cryptic Found!");
            // return CRYPTIC_SPLICE_VALUE;  // Set cryptic splice trigger to 1 if values are greater than CRYPTIC_SPLICE_THRESHOLD
        }
    }

    return 0;
}

float * predict_on_sequence(Model *models, const char *seq, const int width, uint32_t position, Range prediction_boundaries, int strand, Records *rs) {
    int64_t gene_start = prediction_boundaries.start;
    int64_t gene_end = prediction_boundaries.end;
    int64_t vposition = position;

    if (rs != NULL) {
        int start_pred_offset = haplotype_offset(rs, (Range) { vposition - CAI_DISTANCE, vposition });
        gene_start -= start_pred_offset;
        gene_end -= (haplotype_offset(rs, prediction_boundaries) - start_pred_offset);
        vposition -= haplotype_offset(rs, (Range) { prediction_boundaries.start, position });

    }
    int32_t start_offset = VARIANT_POS + (gene_start - vposition);
    int32_t end_offset = VARIANT_POS - (gene_end - vposition);

    const char *shifted_seq = seq - start_offset;
    if (start_offset < 0) start_offset = 0;
    if (end_offset < 0) end_offset = 0;

    Range padding_boundary = { start_offset, end_offset };
    char *padded_seq = pad_sequence(shifted_seq, padding_boundary, width);
    float *predictions; int num_predictions;
    predict_for_seq(models, padded_seq, width, strand, &predictions, &num_predictions);
    free(padded_seq);

    return predictions;
}

void set_record_flag(Record *r, int transcript_index, int phase_index, uint32_t flag_value) {
     r->flags[transcript_index * NUM_PHASES + phase_index] |= flag_value;
}

void get_prediction_pairs_for_site(int site, float *ref_predictions, size_t ref_len, float *hap_alt_predictions, size_t hap_alt_len, Records *rs, Range pred_range, size_t *map_size, PairPrediction **map) {
    float ref_site_predictions[ref_len];
    for (int i = 0; i < ref_len; i++) ref_site_predictions[i] = ref_predictions[NUM_SCORES * i + site];

    float ref_site_prediction_posteriors[ref_len];
    int n_ref_sites = transform_to_posterior_probabilities(ref_len, ref_site_predictions, ref_site_prediction_posteriors);

    float halt_site_predictions[hap_alt_len];
    for (int i = 0; i < hap_alt_len; i++) halt_site_predictions[i] = hap_alt_predictions[NUM_SCORES * i + site];

    float halt_site_prediction_posteriors[hap_alt_len];
    int n_alt_sites = transform_to_posterior_probabilities(hap_alt_len, halt_site_predictions, halt_site_prediction_posteriors);

    PairPrediction *ref_alt_map = NULL;
    size_t ref_alt_map_size = 0;
    parse_prediction_pairs(ref_site_predictions, ref_site_prediction_posteriors, halt_site_predictions, halt_site_prediction_posteriors, rs, pred_range, &ref_alt_map, &ref_alt_map_size);

    *map = ref_alt_map;
    *map_size = ref_alt_map_size;
}

void kv_push_new_probability(SpliceProbabilities *ps, SpliceProbability p) {
    for (int i = 0; i < ps->n; i++) if (p.position == kv_A(*ps, i).position) return;

    kv_push(SpliceProbability, *ps, p);
}

void get_splice_probabilities_from_map(Record *r, enum SpliceType type, int record_phase, int *gene_mask, Range pred_range, size_t map_size, PairPrediction *map, SpliceProbabilities *probabilities) {
    for (int t = 0; t < r->tscripts->n; t++) {
        if (gene_mask[t] == 0) continue;

        Gene *transcript_gene = kv_A(*r->transcript_genes, t);
        Transcript *transcript = kv_A(*r->tscripts, t);

        int canonical_sites[map_size];
        memset(canonical_sites, 0, sizeof(canonical_sites));

        for (int e = 0; e < transcript->nexons; e++) {
            Exon *exon = transcript->exons[e];

            int exon_boundary_position = type == DONOR ? exon->ref_bounds->end : exon->ref_bounds->start;

            if ((type == DONOR) & ((transcript->strand == STRAND_REV && e == 0) | (transcript->strand == STRAND_FWD && e >= transcript->nexons-1))) continue;
            if ((type == ACCEPTOR) & ((transcript->strand == STRAND_REV && e >= transcript->nexons-1) | (transcript->strand == STRAND_FWD && e == 0))) continue;
            if (!((exon_boundary_position >= pred_range.start) & (exon_boundary_position < pred_range.end))) continue;


            int index = (int32_t) exon_boundary_position - pred_range.start;

            for (int i = 0; i < map_size; i++) {
                PairPrediction pair = map[i];
                if (pair.ref.missing) continue;
                if (pair.ref.position != exon_boundary_position) continue;

                canonical_sites[i] = 1;

                if (pair.ref.prior < 0.01) {
                    type == DONOR ? set_record_flag(r, t, record_phase, WEAK_CANONICAL_DONOR) : set_record_flag(r, t, record_phase, WEAK_CANONICAL_ACCEPTOR);

                    if (pair.alt.missing | (pair.alt.prior < 0.01)) continue;

                    // type == DONOR ? set_record_flag(r, t, record_phase, HAPLOTYPE_SPECIFIC_DONOR) : set_record_flag(r, t, record_phase, HAPLOTYPE_SPECIFIC_ACCEPTOR);

                    SpliceProbability p = { r, type, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, exon_boundary_position };
                    kv_push_new_probability(probabilities, p);
                } else {
                    // if (!pair.alt.missing && pair.alt.prior >= 0.01) continue; // No change

                    // If splice site has weakened
                    if (pair.ref.prior - pair.alt.prior > 0.1) {
                        type == DONOR ? set_record_flag(r, t, record_phase, ABLATED_CANONICAL_DONOR) : set_record_flag(r, t, record_phase, ABLATED_CANONICAL_ACCEPTOR);
                    }

                    SpliceProbability p = { r, type, pair.ref.prior, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, exon_boundary_position};
                    kv_push_new_probability(probabilities, p);
                }
            }
        }

        for (int i = 0; i < map_size; i++) {
            if (canonical_sites[i]) continue; // Canonical sites have already been checked, only handle cryptics here

            PairPrediction pair = map[i];

            if (pair.ref.prior == 0) {
                if (pair.alt.prior == 0) continue;

                type == DONOR ? set_record_flag(r, t, record_phase, CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, CRYPTIC_ACCEPTOR);
                // type == DONOR ? set_record_flag(r, t, record_phase, HAPLOTYPE_SPECIFIC_DONOR) : set_record_flag(r, t, record_phase, HAPLOTYPE_SPECIFIC_ACCEPTOR);
            } else {
                if (pair.alt.prior > 0) {
                    type == DONOR ? set_record_flag(r, t, record_phase, CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, CRYPTIC_ACCEPTOR);
                } else {
                    type == DONOR ? set_record_flag(r, t, record_phase, ABLATED_CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, ABLATED_CRYPTIC_ACCEPTOR);
                }
            }

            SpliceProbability p = { r, type, pair.ref.prior, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, 0 };
            kv_push_new_probability(probabilities, p);
        }
    }
}

// int create_splice_probabilities(int vposition, Record *r, Transcripts *transcripts, int *transcript_mask, int num_hap_alt_predictions, float *hap_alt_predictions, SpliceProbability **splice_probabilities) {
//     float *posterior_probabilities = calloc(num_hap_alt_predictions, sizeof(float));
//     int num_acceptors = transform_to_posterior_probabilities(ACCEPTOR_POS, num_hap_alt_predictions, hap_alt_predictions, &posterior_probabilities);
//     int num_donors = transform_to_posterior_probabilities(DONOR_POS, num_hap_alt_predictions, hap_alt_predictions, &posterior_probabilities);
//
//     int num_splice_probabilities = num_acceptors + num_donors;
//
//     int pred_width = num_hap_alt_predictions / NUM_SCORES;
//     Range pred_range = { vposition - CAI_DISTANCE, vposition - CAI_DISTANCE + pred_width };
//     int canonicals[pred_width];
//     // log_info("Length: %i, # transcripts: %i", pred_width, transcripts->n);
//     for (int i = 0; i < pred_width; i++) canonicals[i] = 0;
//     for (int ti = 0; ti < transcripts->n; ti++) {
//         Transcript *t = kv_A(*transcripts, ti);
//
//         if (t == NULL || !transcript_mask[ti]) continue;
//         for (int ei = 0; ei < t->nexons; ei++) {
//             Exon *e = t->exons[ei];
//
//             if (t->strand == STRAND_UNK) {
//                 log_error("UNKNOWN STRAND!");
//                 continue;
//             }
//
//
//             if (((t->strand == STRAND_FWD && ei > 0) || (t->strand == STRAND_REV && ei < t->nexons-1)) && e->bounds->start >= pred_range.start && e->bounds->start < pred_range.end) {
//                 int index = (int32_t) e->bounds->start - vposition + CAI_DISTANCE;
//
//                 // log_info("Acceptor: index=%i, value=%f, ti=%i", index, posterior_probabilities[index * NUM_SCORES + ACCEPTOR], ti);
//
//                 // Only add for loss of canonical site
//                 if(canonicals[index] == 0 && posterior_probabilities[index * NUM_SCORES + ACCEPTOR] == 0) num_splice_probabilities++;
//                 canonicals[index] = ACCEPTOR;
//
//             }
//             if (((t->strand == STRAND_FWD && ei < t->nexons-1) || (t->strand == STRAND_REV && ei > 0)) && e->bounds->end >= pred_range.start && e->bounds->end < pred_range.end) {
//                 int index = (int32_t) e->bounds->end - vposition + CAI_DISTANCE;
//
//                 // log_info("Donor: index=%i, value=%f, ti=%i", index, posterior_probabilities[index * NUM_SCORES + DONOR], ti);
//
//                 // Only add for loss of canonical site
//                 if(canonicals[index] == 0 && posterior_probabilities[index * NUM_SCORES + DONOR] == 0) num_splice_probabilities++;
//                 canonicals[index] = DONOR;
//             }
//         }
//     }
//     // log_info("Num splice probabilities: %i", num_splice_probabilities);
//
//     SpliceProbability *probabilities = malloc(num_splice_probabilities * sizeof(SpliceProbability));
//     int prob_idx = 0;
//     for (int i = 0; i < num_hap_alt_predictions; i++) {
//         if (i % NUM_SCORES == 0) continue; // Skip predictions for neither
//
//         float prior_probability = hap_alt_predictions[i];
//         float posterior_probability = posterior_probabilities[i];
//         if (canonicals[i / NUM_SCORES] == DONOR && i % NUM_SCORES == DONOR) {
//             int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
//             int prediction_position = vposition + offset_from_variant;
//             probabilities[prob_idx++] = (SpliceProbability) { r, DONOR, prior_probability, posterior_probability, prediction_position, prediction_position };
//             continue;
//         } else if (canonicals[i / NUM_SCORES] == ACCEPTOR && i % NUM_SCORES == ACCEPTOR) {
//             int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
//             int prediction_position = vposition + offset_from_variant;
//             probabilities[prob_idx++] = (SpliceProbability) { r, ACCEPTOR, prior_probability, posterior_probability, prediction_position, prediction_position };
//             continue;
//         }
//
//         if (posterior_probability == 0) continue; // Skip null predictions
//
//         int type = i % 3;
//         int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
//         int prediction_position = vposition + offset_from_variant;
//
//         // log_info("Type: %i, offset: %i, prediction_position: %i, prior: %f, posterior: %f", type, offset_from_variant, prediction_position, prior_probability, posterior_probability);
//         probabilities[prob_idx++] = (SpliceProbability) { r, type, prior_probability, posterior_probability, prediction_position, 0 };
//     }
//
//     *splice_probabilities = probabilities;
//     return num_splice_probabilities;
// }

void resize_splice_sites(SpliceSite **sites_array, size_t *current_size) {
    *current_size *= 2;
    *sites_array = realloc(*sites_array, (*current_size) * sizeof(SpliceSite));
}

int check_proximity(uint32_t predicted_position, uint32_t canonical_site_position) {
    return abs((int32_t)predicted_position - (int32_t)canonical_site_position) <= PROXIMITY_THRESHOLD;
}

void kv_push_new_record(Records *rs, Record *r) {
    for (int i = 0; i < rs->n; i++) if (r == kv_A(*rs, i)) return;

    kv_push(Record *, *rs, r);
}

void store_splice_sites(SpliceProbability prediction, int *predict_for_pair, int *predict_for_competitor, SpliceSites *splice_sites) {
    if (prediction.type == DONOR && prediction.posterior_probability == 0 && prediction.position == splice_sites->canonical_donors[0].position) {
        kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
        splice_sites->canonical_donors[0].predicted = 1;
        splice_sites->canonical_donors[0].ref_prior_probability = prediction.ref_prior_probability;
        splice_sites->canonical_donors[0].prior_probability = prediction.prior_probability;
        splice_sites->canonical_donors[0].posterior_probability = prediction.posterior_probability;
        splice_sites->canonical_donors[0].lost = 1;
        *predict_for_pair = 1;
        *predict_for_competitor = 1;
    } else if (prediction.type == ACCEPTOR && prediction.posterior_probability == 0 && prediction.position == splice_sites->canonical_acceptors[0].position) {
        kv_push_new_record((Records *) splice_sites->canonical_acceptors[0].records_ptr, prediction.record);
        splice_sites->canonical_acceptors[0].predicted = 1;
        splice_sites->canonical_acceptors[0].ref_prior_probability = prediction.ref_prior_probability;
        splice_sites->canonical_acceptors[0].prior_probability = prediction.prior_probability;
        splice_sites->canonical_acceptors[0].posterior_probability = prediction.posterior_probability;
        splice_sites->canonical_acceptors[0].lost = 1;
        *predict_for_pair = 1;
        *predict_for_competitor = 1;
    } else if (prediction.type == DONOR & prediction.posterior_probability > 0) {
        if (splice_sites->canonical_donors[0].gff_canonical && prediction.posterior_probability == 1.0 && prediction.position == splice_sites->canonical_donors[0].position) { // If canonical exists and is the only splice site within the window
            kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
            splice_sites->canonical_donors[0].predicted = 1;
            splice_sites->canonical_donors[0].ref_prior_probability = prediction.ref_prior_probability;
            splice_sites->canonical_donors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_donors[0].posterior_probability = prediction.posterior_probability;
        } else if (splice_sites->canonical_donors[0].gff_canonical && check_proximity(prediction.position, splice_sites->canonical_donors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical donor site
            int canonical_donor_exists = 0;
            int canonical_donor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_donors; i++) {
                if (splice_sites->canonical_donors[i].position == prediction.position) {
                    canonical_donor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->canonical_donors[i].records_ptr, prediction.record);
                    if (splice_sites->canonical_donors[i].predicted == 1) {
                        canonical_donor_predicted = 1;
                    }
                }
            }

            if (canonical_donor_exists & canonical_donor_predicted) return;
            if (canonical_donor_exists & !canonical_donor_predicted) {
                *predict_for_pair = 1;
                for (int i = 0; i < splice_sites->n_canonical_donors; i++) {
                    if (splice_sites->canonical_donors[i].position == prediction.position) {
                        kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
                        splice_sites->canonical_donors[i].canonical_position = prediction.canonical_position;
                        splice_sites->canonical_donors[i].predicted = 1;
                        splice_sites->canonical_donors[i].ref_prior_probability = prediction.ref_prior_probability;
                        splice_sites->canonical_donors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_donors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_donors >= splice_sites->m_canonical_donors) {
                    resize_splice_sites(&splice_sites->canonical_donors, &splice_sites->m_canonical_donors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_donors[splice_sites->n_canonical_donors++] = site;
            }
        // Store cryptic donor site
        } else {
            // Check if prediction position already exists as a cryptic donor site
            int cryptic_donor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_donors; i++) {
                if (splice_sites->cryptic_donors[i].position == prediction.position) {
                    cryptic_donor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->cryptic_donors[i].records_ptr, prediction.record);
                }
            }

            if (!cryptic_donor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_donors >= splice_sites->m_cryptic_donors) {
                    resize_splice_sites(&splice_sites->cryptic_donors, &splice_sites->m_cryptic_donors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);

                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_donors[splice_sites->n_cryptic_donors++] = site;
            }
        }
    } else if (prediction.type == ACCEPTOR & prediction.posterior_probability > 0) {
        // Store canonical acceptor site
        if (splice_sites->canonical_acceptors[0].gff_canonical && prediction.posterior_probability == 1.0 && prediction.position == splice_sites->canonical_acceptors[0].position) { // If canonical exists and is the only splice site within the window
            kv_push_new_record((Records *) splice_sites->canonical_acceptors[0].records_ptr, prediction.record);
            splice_sites->canonical_acceptors[0].predicted = 1;
            splice_sites->canonical_acceptors[0].ref_prior_probability = prediction.ref_prior_probability;
            splice_sites->canonical_acceptors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_acceptors[0].posterior_probability = prediction.posterior_probability;
        } else if (splice_sites->canonical_acceptors[0].gff_canonical & check_proximity(prediction.position, splice_sites->canonical_acceptors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical acceptor site
            int canonical_acceptor_exists = 0;
            int canonical_acceptor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_acceptors; i++) {
                if (splice_sites->canonical_acceptors[i].position == prediction.position) {
                    canonical_acceptor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->canonical_acceptors[i].records_ptr, prediction.record);
                    if (splice_sites->canonical_acceptors[i].predicted == 1) {
                        canonical_acceptor_predicted = 1;
                    }
                }
            }

            if (canonical_acceptor_exists & canonical_acceptor_predicted) return;
            if (canonical_acceptor_exists & !canonical_acceptor_predicted) {
                *predict_for_pair = 1;
                for (int i = 0; i < splice_sites->n_canonical_acceptors; i++) {
                    if (splice_sites->canonical_acceptors[i].position == prediction.position) {
                        kv_push_new_record((Records *) splice_sites->canonical_acceptors[i].records_ptr, prediction.record);
                        // splice_sites->canonical_acceptors[i].canonical_position = prediction.canonical_position;
                        splice_sites->canonical_acceptors[i].predicted = 1;
                        splice_sites->canonical_acceptors[i].ref_prior_probability = prediction.ref_prior_probability;
                        splice_sites->canonical_acceptors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_acceptors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_acceptors >= splice_sites->m_canonical_acceptors) {
                    resize_splice_sites(&splice_sites->canonical_acceptors, &splice_sites->m_canonical_acceptors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_acceptors[splice_sites->n_canonical_acceptors++] = site;
            }
        // Store cryptic acceptor site
        } else {
            // Check if prediction position already exists as a cryptic acceptor site
            int cryptic_acceptor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_acceptors; i++) {
                if (splice_sites->cryptic_acceptors[i].position == prediction.position) {
                    cryptic_acceptor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->cryptic_acceptors[i].records_ptr, prediction.record);
                }
            }

            if (!cryptic_acceptor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_acceptors >= splice_sites->m_cryptic_acceptors) {
                    resize_splice_sites(&splice_sites->cryptic_acceptors, &splice_sites->m_cryptic_acceptors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_acceptors[splice_sites->n_cryptic_acceptors++] = site;
            }
        }
    }
}

int has_af_trigger(const bcf1_t *record, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold) {
    float *af_value = NULL;
    int ndst = 0;

    // Set based on INFO/AF
    if (bcf_get_info_float(hdr, (bcf1_t *) record, af_tag, &af_value, &ndst) <= 0) return 0;

    return *af_value <= af_threshold;
}

int has_mnv_trigger(bcf1_t *record, bcf_hdr_t *hdr) {
    int32_t *mnv = NULL;
    int ndst = 0;

    // Get based on FORMAT/MNV
    if (bcf_get_format_int32(hdr, record, "MNV", &mnv, &ndst) <= 0) return 0;

    return *mnv;
}

int is_purine(char base) {
    return base == 'A' || base == 'G';
}

int calculate_kozak_strength(char* sequence, int pos, int length) {
    int score = 0;

    if (pos >= 3 && is_purine(sequence[pos - 3])) {
        score += 1;  // Give a point for purine at -3
    }

    if (pos + 3 < length && sequence[pos + 3] == 'G') {
        score += 1;  // Give a point for G at +4
    }

    return score;
}

int32_t* find_possible_reinitiation_positions(char* reinitiation_sequence, int* count, int** kozak_strengths, int offset) {
    int length = strlen(reinitiation_sequence);
    char codon[4] = {0}; // Buffer to hold codon (3 bases + null terminator)
    int32_t* reinitiation_positions = NULL;
    *kozak_strengths = NULL;
    int allocated_size = 0;
    *count = 0;

    for (int i = 0; i < length; i += 3) {
        codon[0] = reinitiation_sequence[i];
        codon[1] = reinitiation_sequence[i + 1];
        codon[2] = reinitiation_sequence[i + 2];
        codon[3] = '\0';

        if (strcmp(codon, START_CODON) == 0) {
            int strength = calculate_kozak_strength(reinitiation_sequence, i, length);

            if (*count >= allocated_size) {
                allocated_size = (allocated_size == 0) ? 1 : allocated_size * 2;
                reinitiation_positions = realloc(reinitiation_positions, allocated_size * sizeof(int32_t));
                *kozak_strengths = realloc(*kozak_strengths, allocated_size * sizeof(int));
            }

            reinitiation_positions[*count] = i + 1 + offset; // Store 1-based position in the actual sequence (i + 1 + offset)
            (*kozak_strengths)[*count] = strength;
            (*count)++;
        }
    }
    return reinitiation_positions;
}

void get_prediction_flags(uint32_t flags, kstring_t *concatenated_flags) {
    const char *flag_names[] = {
        "unphased_variant_in_prediction",
        "skipped_variant_in_prediction",
        "variant_larger_than_splice_window",
        "rescued_by_splice_frameshift",
        "possible_reinitiation",
        "manufactured_splice_rescue",
        "haplotype_specific_donor",
        "haplotype_specific_acceptor",
        "cryptic_donor",
        "cryptic_acceptor",
        "weak_canonical_donor",
        "weak_canonical_acceptor",
        "ablated_canonical_donor",
        "ablated_canonical_acceptor",
        "ablated_cryptic_donor",
        "ablated_cryptic_acceptor",
        "upstream_ptv",
        "5utr_variant",
        "3utr_variant",
        "last_exon",
        "single_exon",
        "non_constitutive",
        "in_frame"
    };

    size_t flag_count = sizeof(flag_names) / sizeof(flag_names[0]);
    int is_first = 1;

    for (size_t i = 0; i < flag_count; i++) {
        if (flags & (1 << i)) {
            if (!is_first) {
                kputs(":", concatenated_flags); // append delimiter if not the first flag
            }
            kputs(flag_names[i], concatenated_flags); // append the flag name
            is_first = 0;
        }
    }
}

int buffer_record(htsFile *fp, bcf_hdr_t *hdr, RecordBuffer *buffer) {
    if (buffer->unused) return 0; // Record hasn't been processed yet, no need to load in a new one

    int next = bcf_read(fp, hdr, buffer->v);

    if (next < 0) {
        // No next line, but no errorcode means end of file has been reached.
        if (buffer->v->errcode == 0) return 1;

        log_error("Unhandled critical error encountered during parsing of line: %s", fp->line.s);
        return 1;
    }
    if (next > 0) {
        // TODO: These should be handled properly, but no time to implement this.
        log_warn("Unhandled non-critical error encountered during parsing of line: %s", fp->line.s);
        return 0;
    }

    buffer->unused = true;

    return 0;
}

int get_next_record_block(htsFile *fp, bcf_hdr_t *hdr, HdrFormat *vep_fmt, gff_t *gff, RecordBuffer *buffer, RecordBlock2 *block_out) {
    RecordBlock2 block = {0};
    kv_init(block.records);

    IRegion batch_region = {0};
    while(buffer_record(fp, hdr, buffer) == 0) {

        // Check if record is coding
        IRegion record_region = {0};
        int is_non_coding = find_gene_region_for_variant(buffer->v, hdr, gff, &record_region);
        if (is_non_coding) {
            if (block.records.n == 0) {
                kv_push(Record *, block.records, init_record(buffer->v, NULL));
                buffer->unused = false;
            }

            *block_out = block;
            return 0;
        }


        // Check if the coding region the record falls in overlaps the current batch
        block.is_coding = 1;
        int overlaps_batch = 1;
        if (block.records.n == 0) {
            batch_region = record_region; // First record in batch has nothing to overlap with
        } else {
            overlaps_batch = iregion_intersects(batch_region, record_region);
            if (!overlaps_batch) {
                *block_out = block;
                return 0;
            }

            // Extends the batch region to allow for overlapping genes
            if (overlaps_batch & (record_region.end > batch_region.end)) batch_region.end = record_region.end; 
        }

        VepCsqs *vep_csqs = parse_vep_transcript_consequences(buffer->v, vep_fmt, hdr);
        kv_push(Record *, block.records, init_record(buffer->v, vep_csqs));
        buffer->unused = false;
    }

    if (block.records.n > 0) {
        *block_out = block;
        return 0;
    }

    return 1;
}

void destroy_record_block(RecordBlock2 batch) {
    for (int i = 0; i < batch.records.n; i++) {
        destroy_record(kv_A(batch.records, i));
    }
    kv_destroy(batch.records);
}

void set_coding_transcripts_genes(Record *r, bcf_hdr_t *hdr, gff_t *gff, faidx_t *fai, kh_TranscriptHash_t *transcripts, kh_GeneHash_t *genes) {
    gf_tscripts *tscripts = NULL;
    find_gf_tscripts_for_record(r, hdr, gff, &tscripts);
    if (tscripts == NULL || tscripts->n == 0) return; // A variant in a coding gene can still have no overlapping coding transcripts

    r->tscripts = malloc(sizeof(Transcripts));
    kv_init(*r->tscripts);
    for (int j = 0; j < tscripts->n; j++) {
        gf_tscript_t *gf_t = kv_A(*tscripts, j);
        if (gf_t == NULL) { kv_push(Transcript *, *r->tscripts, NULL); continue; }

        int key = gf_t->id;
        // khint32_t hash_index = kh_get_TranscriptHash(transcripts, key);
        // if (hash_index != kh_end(transcripts)) {
        //     Transcript *t = kh_val(transcripts, hash_index);
        //     kv_push(Transcript *, *r->tscripts, t);
        //     continue;
        // }

        Transcript *t = malloc(sizeof(Transcript));
        t->id = key;
        t->start = gf_t->beg;
        t->end = gf_t->end;
        t->strand = gf_t->strand;

        t->mcds = gf_t->mcds;
        t->ncds = gf_t->ncds;
        t->cds = malloc(t->mcds * sizeof(Range));
        int constitutives[t->mcds];
        for (int k = 0; k < t->ncds; k++) {
            uint32_t start = gf_t->cds[k]->beg;
            uint32_t end = gf_t->cds[k]->beg + (gf_t->cds[k]->len - 1);
            t->cds[k] = (Range) { start, end };

            constitutives[k] = gf_t->cds[k]->constitutive;
        }
        t->triggers[0] = 0;
        t->triggers[1] = 0;
        t->triggers[2] = 0;

        const char *chr = bcf_hdr_id2name(hdr, r->v->rid);
        set_transcript_exons(t, chr, gff, (int *) constitutives);
        kv_push(Transcript *, *r->tscripts, t);

        int is_new;
        khint32_t hash_index = kh_put_TranscriptHash(transcripts, key, &is_new);
        if (!is_new) continue;
        kh_val(transcripts, hash_index) = t;
    }

    // Create Genes
    r->transcript_genes = malloc(sizeof(Genes));
    kv_init(*r->transcript_genes);

    r->genes = malloc(sizeof(Genes));
    kv_init(*r->genes);
    for (int j = 0; j < tscripts->n; j++) {
        gf_tscript_t *t = kv_A(*tscripts, j);
        if (t == NULL) { kv_push(Gene *, *r->transcript_genes, NULL); continue; }

        khint_t hash_index;
        hash_index = kh_get_GeneHash(genes, (khint32_t) t->gene->id);
        if (hash_index != kh_end(genes)) {
            Gene *gene = kh_val(genes, hash_index);

            int has_gene = 0;
            for (int g = 0; g < r->genes->n; g++) {
                if (gene->gff_id != kv_A(*r->genes, g)->gff_id) continue;

                has_gene = 1;
                break;
            }

            if (!has_gene) kv_push(Gene *, *r->genes, gene);

            kv_push(Gene *, *r->transcript_genes, gene);
            continue;
        }

        // No match found, create and store a new gene
        Range gene_range = {t->gene->beg, t->gene->end};

        int slen = 0;
        const char *chr = bcf_hdr_id2name(hdr, r->v->rid);
        char *seq = faidx_fetch_seq(fai, chr, gene_range.start, gene_range.end, &slen);
        seq[slen] = '\0';
        kstring_t s = { slen, slen+1, seq };

        Gene *gene = malloc(sizeof(Gene));
        gene->gff_id = t->gene->id;
        gene->name = t->gene->name;
        gene->coords = gene_range;
        gene->seq = s;
        gene->strand = t->gene->strand;

        // Add gene to hashmap
        khint32_t key = gene->gff_id;
        int is_new;
        hash_index = kh_put_GeneHash(genes, key, &is_new);
        if (!is_new) {
            log_error("Failed to add %s due to hash conflict.", gene->name);
        }
        kh_val(genes, hash_index) = gene;

        kv_push(Gene *, *r->genes, gene);
        kv_push(Gene *, *r->transcript_genes, gene);
    }

    free(tscripts);

    // Initialize scores, isoforms and flags
    // TODO: This should be moved, doesn't belong in this function. This should be done only once per record, whereas transcripts/exons carry phase specific data. It's too messy.
    if (r->haplotype_scores != NULL) return;

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

int has_matching_phase(const bcf1_t *v, const bcf_hdr_t *hdr, const int phase, const bool include_unphased) {
    int *gt = NULL; int ngts = 0;
    bcf_get_genotypes(hdr, (bcf1_t *) v, &gt, &ngts);

    // Homozygous
    if (bcf_gt_allele(gt[0]) == bcf_gt_allele(gt[1])) {
        free(gt);
        return 1;
    }

    // Unphased
    if (!bcf_gt_is_phased(gt[0]) & !bcf_gt_is_phased(gt[1])) {
        free(gt);
        return include_unphased;
    }

    int is_present = bcf_gt_allele(gt[phase]);
    free(gt);
    return is_present;
}

void set_vep_triggers(Records records, int transcript_id) {
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);

        for (int vi = 0; vi < r->vep->n; vi++) {
            TscriptCsqs *transcript_consequence = kv_A(*r->vep, vi);
            Transcript *transcript = kv_A(*r->tscripts, vi);

            if (transcript == NULL || transcript->id != transcript_id) continue;

            for (int ci = 0; ci < transcript_consequence->num_consequences; ci++) {
                const char *consequence = transcript_consequence->consequences[ci];
                if ((strcmp(consequence, "non_coding_transcript_variant") == 0) | (strcmp(consequence, "NMD_transcript_variant") == 0)) continue;
                if ((strcmp(consequence, SPLICE_DONOR_VARIANT) == 0) | (strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0)) transcript_consequence->triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
                if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) transcript_consequence->triggers[FRAMESHIFT_TRIGGER] = 1;
                if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) transcript_consequence->triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;
            }
        }
    }
}

void set_af_triggers(Records records, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold) {
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);
        int triggered = has_af_trigger(r->v, (bcf_hdr_t *) hdr, af_tag, af_threshold);

        if (!triggered) continue;

        for (int vi = 0; vi < r->vep->n; vi++) {
            TscriptCsqs *transcript_consequences = kv_A(*r->vep, vi);

            if (transcript_consequences->triggers[SPLICE_TRIGGER]) continue;

            transcript_consequences->triggers[SPLICE_TRIGGER] = AF_SPLICE_VALUE;
        }
    }
}

// TODO: Rename this to something other than MNV, as it really just checks for possible reading frame conflicts not checked by VEP
void set_mnv_triggers(Records records, Transcript *transcript) {
    size_t n_exons = transcript->nexons;

    int intron_sizes[n_exons-1];
    for (int ei = 1; ei < n_exons; ei++) {
        Range *prev_exon = transcript->exons[ei-1]->ref_bounds;
        Range *curr_exon = transcript->exons[ei]->ref_bounds;
        intron_sizes[ei-1] = curr_exon->start - prev_exon->end;
    }

    // Figure out the positions in the coding transcript
    int64_t transcript_start = transcript->start;
    int64_t positions[records.n];
    int position_record_map[records.n];
    size_t n_positions = 0;
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);
        int64_t record_pos = r->v->pos;

        bool is_coding = true;
        int skipped_intron_bases = 0;
        for (int ei = 0; ei < n_exons; ei++) {
            Range exon = *transcript->exons[ei]->ref_bounds;

            // Either intron or UTR, neither is coding
            if (record_pos < exon.start) {
                is_coding = false;
                break;
            }

            if (record_pos > exon.end) {
                // UTR if beyond last exon
                if (ei == n_exons-1) {
                    is_coding = false;
                    break;
                }

                skipped_intron_bases += intron_sizes[ei];
                continue;
            }

            break; // Falls within current exon, we're done!
        }

        if (!is_coding) continue;

        position_record_map[n_positions] = ri;
        positions[n_positions++] = record_pos - transcript_start - skipped_intron_bases;
    }

    if (!n_positions) return; // No coding positions, nothing to be done!

    for (int pi = 0; pi < n_positions-1; pi++) { // n_positions is size_t, so without the if size 0 -> return this can overflow
        int64_t curr = positions[pi];
        int64_t next = positions[pi+1];

        // Cannot be inframe!
        if (next - curr > 2) continue;

        int32_t transcript_id = transcript->id;

        // Find the record for the current position
        Record *curr_record = kv_A(records, position_record_map[pi]);
        for (int ti = 0; ti < curr_record->tscripts->n; ti++) {
            Transcript *record_transcript = kv_A(*curr_record->tscripts, ti);

            if (record_transcript == NULL) continue;
            if (transcript_id != record_transcript->id) continue;
            if (record_transcript->triggers[STOP_TRIGGER]) continue; // Already set!

            transcript->triggers[STOP_TRIGGER] = STOP_GAIN_VALUE;
        }

        Record *next_record = kv_A(records, position_record_map[pi+1]);
        for (int ti = 0; ti < next_record->tscripts->n; ti++) {
            Transcript *record_transcript = kv_A(*next_record->tscripts, ti);

            if (record_transcript == NULL) continue;
            if (transcript_id != record_transcript->id) continue;
            if (record_transcript->triggers[STOP_TRIGGER]) continue; // Already set!

            transcript->triggers[STOP_TRIGGER] = STOP_GAIN_VALUE;
        }
    }
}

void set_triggers(Record *r) {
    int n_tscripts = r->tscripts->n;
    for (int i = 0; i < n_tscripts; i++) {
        Gene *g = kv_A(*r->genes, i);
        Transcript *t = kv_A(*r->tscripts, i);

    }


    log_info("");
    return;
}

void write_predictions(FILE *scores_tsv, const char *chr, const int64_t pos, const char *ref, const char *alt, const char *gene_name, const char *transcript_name, const int phase, const char *type, const char *pred_tags[], const int n_predictions, const float *predictions[]) {
    for (int i = 0; i < 101; i++) {
        const int offset = i - 50;

        bool should_print = false;
        for (int j = 0; j < n_predictions; j++) {
            const float acceptor = predictions[j][i * 3 + ACCEPTOR];
            const float donor = predictions[j][i * 3 + DONOR];

            if ((acceptor > SPLICEAI_RESCUE_THRESHOLD) | (donor > SPLICEAI_RESCUE_THRESHOLD)) {
                should_print = true;
            }

        }
        for (int j = 0; j < n_predictions; j++) {
            if (!should_print) continue;

            const float neither = predictions[j][i * 3];
            const float acceptor = predictions[j][i * 3 + ACCEPTOR];
            const float donor = predictions[j][i * 3 + DONOR];

            fprintf(scores_tsv, "%s\t%li\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
                    chr, pos + 1, ref, alt, gene_name, transcript_name, phase+1, offset, pred_tags[j], type, donor, acceptor, neither);
        }
    }
}

void print_usage(const char *prog_name) {
    printf("Usage: %s -m <models> -g <gff> -r <reference> -i <in_vcf> -o <out_vcf> [options]\n", prog_name);
    printf("\nRequired arguments:\n");
    printf("  -m, --models        Path to SpliceAI models\n");
    printf("  -g, --gff           Path to GFF file\n");
    printf("  -r, --reference     Path to reference FASTA file\n");
    printf("  -i, --in            Path to input VCF file\n");
    printf("  -o, --out           Path to output VCF file\n");
    printf("  -s, --scores        Path to splice scores file\n");

    printf("\nOptional arguments:\n");
    printf("      --use-af        Use non-default field for allele frequency thresholding (default: INFO/AF)\n");
    printf("      --use-gt        Use non-default field for phased genotypes (default: FORMAT/GT)\n");
    printf("      --af-threshold  Allele frequency threshold for splice score annotation (default: 0.01)\n");
    printf("      --include-unphased  Include unphased variants (flag)\n");
    printf("      --max-unphased  Maximum number of unphased variants allowed (requires --include-unphased) (default: 2)\n");
    printf("  -h, --help          Show this help message\n");
}

int parse_args(int argc, char *argv[], Arguments *args) {

    // Set default values
    args->use_af = AF_FIELD;
    args->use_gt = "GT";
    args->af_threshold = AF_THRESHOLD;
    // args->max_unphased = 2;
    args->include_unphased = 0;

    struct option long_options[] = {
        {"models", required_argument, 0, 'm'},
        {"gff", required_argument, 0, 'g'},
        {"reference", required_argument, 0, 'r'},
        {"in", required_argument, 0, 'i'},
        {"out", required_argument, 0, 'o'},
        {"scores", required_argument, 0, 's'},
        {"use-af", required_argument, 0, 0},
        {"use-gt", required_argument, 0, 0},
        {"af-threshold", required_argument, 0, 0},
        {"include-unphased", no_argument, &args->include_unphased, 1},
        // {"max-unphased", required_argument, 0, 0},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    if (argc == 1) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int option_index = 0;
    int c;

    while ((c = getopt_long(argc, argv, "m:g:r:i:o:s:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'm': args->models = optarg; break;
            case 'g': args->gff = optarg; break;
            case 'r': args->reference = optarg; break;
            case 'i': args->in = optarg; break;
            case 'o': args->out = optarg; break;
            case 's': args->scores = optarg; break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            case 0:
                if (strcmp("use-af", long_options[option_index].name) == 0)
                    args->use_af = optarg;
                else if (strcmp("use-gt", long_options[option_index].name) == 0)
                    args->use_gt = optarg;
                else if (strcmp("af-threshold", long_options[option_index].name) == 0)
                    args->af_threshold = strtof(optarg, NULL);
                // else if (strcmp("max-unphased", long_options[option_index].name) == 0)
                //     args->max_unphased = atoi(optarg);
                break;
            default:
                print_usage(argv[0]);
                return EXIT_FAILURE;
        }
    }

    // Check required arguments
    if (!args->models || !args->gff || !args->reference || !args->in || !args->out || !args->scores) {
        fprintf(stderr, "Error: Missing required arguments\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    log_set_level(3);
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);

    // Parse arguments
    Arguments args = {0};
    int status = parse_args(argc, argv, &args);

    if (status != EXIT_SUCCESS) {
        return status;
    }

    // log_debug("Loading SpliceAI models...");
    Model *models = load_models(args.models);

    // Parse gff file
    const char *gff_path = args.gff;
    log_debug("Loading GFF from: %s", gff_path);
    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);

    // Open fasta file
    const char *fasta_path = args.reference;
    log_debug("Loading FASTA from: %s", fasta_path);
    faidx_t *fai = fai_load(fasta_path);

    // Open input VCF
    const char *input_vcf_path = args.in;
    log_debug("Loading VCF from: %s", input_vcf_path);
    htsFile *fp_in = hts_open(input_vcf_path, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp_in);
    HdrFormat *vep_fmt = parse_csq_hdr_fmt(hdr);

    // Open output VCF
    const char *output_vcf_path = args.out;
    log_debug("Saving output to: %s", output_vcf_path);
    htsFile *fp_out = hts_open(output_vcf_path, "w");
    add_rescueranger_tags(hdr);
    bcf_hdr_write(fp_out, hdr);

    // Open raw scores file
    const char *output_tsv_path = args.scores;
    log_debug("Saving raw scores to: %s", output_tsv_path);
    FILE *scores_tsv = fopen(output_tsv_path, "w");
    fprintf(scores_tsv, "CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tPHASE\tOFFSET\tSOURCE\tTYPE\tDONOR\tACCEPTOR\tNEITHER\n");

    RecordBuffer buffer = { bcf_init(), false };
    RecordBlock2 block = { 0 };

    while(get_next_record_block(fp_in, hdr, vep_fmt, gff, &buffer, &block) == 0) {
        if (!block.is_coding) {
            for (int i = 0; i < block.records.n; i++) {
                bcf_write(fp_out, hdr, kv_A(block.records, i)->v);
            }
            destroy_record_block(block);
            continue;
        }

        // Split on phase
        // This is done before setting transcripts and genes so that each has their own copy for isoforms and flags
        for (int phase = 0; phase < NUM_PHASES; phase++) {

            // Find matching gff transcripts/genes for records in batch
            khash_t(GeneHash) *genes = kh_init_GeneHash();
            khash_t(TranscriptHash) *transcripts = kh_init_TranscriptHash();
            for (int i = 0; i < block.records.n; i++) {
                Record *r = kv_A(block.records, i);
                set_coding_transcripts_genes(r, hdr, gff, fai, transcripts, genes);
            }


            // Triggers need to checked per transcript
            for (khiter_t k = kh_begin(transcripts); k != kh_end(transcripts); ++k) {
                if (!kh_exist(transcripts, k)) continue;

                Transcript *transcript = kh_val(transcripts, k);

                // Find all records for current transcript
                Records records;
                kv_init(records);
                for (int ri = 0; ri < block.records.n; ri++) {
                    Record *r = kv_A(block.records, ri);

                    if (r->tscripts == NULL) continue;
                    if (!has_matching_phase(r->v, hdr, phase, args.include_unphased)) continue;

                    for (int ti = 0; ti < r->tscripts->n; ti++) {
                        Transcript *record_transcript = kv_A(*r->tscripts, ti);
                        if (record_transcript == NULL || record_transcript->id != transcript->id) continue;

                        kv_push(Record *, records, r);
                    }
                }

                set_vep_triggers(records, transcript->id);
                set_af_triggers(records, hdr, args.use_af, args.af_threshold);
                set_mnv_triggers(records, transcript);

                kv_destroy(records);
            }

            kh_destroy_GeneHash(genes);
            kh_destroy_TranscriptHash(transcripts);

            // Get records for phase
            Records records;
            kv_init(records);

            for (int ri = 0; ri < block.records.n; ri++) {
                Record *r = kv_A(block.records, ri);

                if (r->tscripts == NULL) continue;
                if (!has_matching_phase(r->v, hdr, phase, args.include_unphased)) continue;

                kv_push(Record *, records, r);
            }

            filter_overlapping_variants(&records, hdr);

            // TODO: Refactor everything below this point into something with actual structure
            for (int i = 0; i < records.n; i++) {
                Record *r = kv_A(records, i);
                // log_info("%d:%d Initialize splice sites", r->v->rid, r->v->pos);

                for (int ti = 0; ti < r->tscripts->n; ti++) {
                    Transcript *t = kv_A(*r->tscripts, ti);
                    Gene *gene = kv_A(*r->transcript_genes, ti);
                    if (t == NULL || gene == NULL) continue;


                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;

                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    for (int j = 0; j < records.n; j++) {
                        Record *hap_r = kv_A(records, j);
                        if (hap_r->genes == NULL) continue;
                        Range record_range = variant_to_range(hap_r->v, 0);

                        for (int gi2 = 0; gi2 < hap_r->genes->n; gi2++) {
                            Gene *gene2 = kv_A(*hap_r->genes, gi2);
                            if (gene2 == NULL) continue;

                            if (gene->gff_id == gene2->gff_id) {
                                Range intersection;
                                intersect_ranges(gene->coords, record_range, &intersection);
                                if (range_len(intersection) == 0) continue;

                                kv_push(Record *, hap_in, hap_r);

                            }
                        }
                    }

                    for (int ei = 0; ei < t->nexons; ei++) {
                        Exon *e = t->exons[ei];
                        if (e->splice_sites != NULL) continue;
                        int kept_start = 1, kept_end = 1;

                        // Checks if donor/acceptor sites were lost due to a deletion  TODO: Check effect of +1 here in gdb
                        Range canonical_start = (Range) { e->ref_bounds->start, e->ref_bounds->start }; // +1?
                        Range canonical_end = (Range) { e->ref_bounds->end, e->ref_bounds->end }; // +1?
                        for (int ri = 0; ri < hap_in.n; ri++) {
                            Record *hap_r = kv_A(hap_in, ri);
                            int is_deletion = allele_offset(hap_r->v, 0, 1, gene->coords) > 0;
                            if (!is_deletion) continue;

                            Range ref = variant_to_range(r->v, 0);

                            Range intersection;
                            if (intersect_ranges(ref, canonical_start, &intersection) && range_len(intersection) > 0) kept_start = 0;
                            if (intersect_ranges(ref, canonical_end, &intersection) && range_len(intersection) > 0) kept_end = 0;
                        }

                        uint32_t exon_hap_start = e->ref_bounds->start - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->ref_bounds->start - 1});
                        uint32_t exon_hap_end = e->ref_bounds->end - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->ref_bounds->end - 2});

                        SpliceSites *sites = malloc(sizeof(SpliceSites));
                        init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
                        init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
                        init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
                        init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

                        e->hap_bounds = malloc(sizeof(Range));

                        if (t->strand == STRAND_FWD) {
                            if (ei > 0) {
                                Records *splice_records = malloc(sizeof(Records));
                                kv_init(*splice_records);
                                sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { splice_records, e->ref_bounds->start, e->ref_bounds->start, 0, kept_start, 0, 1.0, 1.0, 1.0 };
                            }
                            if (ei < t->nexons-1) {
                                Records *splice_records = malloc(sizeof(Records));
                                kv_init(*splice_records);
                                sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { splice_records, e->ref_bounds->end, e->ref_bounds->end, 0, kept_end, 0, 1.0, 1.0, 1.0 };
                            }

                            *e->hap_bounds = (Range) {exon_hap_start, exon_hap_end};
                        } else if (t->strand == STRAND_REV) {
                            if (ei < t->nexons-1) {
                                Records *splice_records = malloc(sizeof(Records));
                                kv_init(*splice_records);
                                sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { splice_records, e->ref_bounds->end, e->ref_bounds->end, 0, kept_end, 0, 1.0, 1.0, 1.0 };
                            }
                            if (ei > 0) {
                                Records *splice_records = malloc(sizeof(Records));
                                kv_init(*splice_records);
                                sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { splice_records, e->ref_bounds->start, e->ref_bounds->start, 0, kept_start, 0, 1.0, 1.0, 1.0 };
                            }

                            *e->hap_bounds = (Range) {exon_hap_end, exon_hap_start};
                            *e->ref_bounds = (Range) {e->ref_bounds->end, e->ref_bounds->start};
                        }

                        e->splice_sites = sites;

                        // Convert CDS boundaries to haplotype positions
                        if (e->cds == NULL) continue;
                        // uint32_t cds_hap_start = e->cds->start - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->cds->start });
                        // uint32_t cds_hap_end = e->cds->end - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->cds->end });
                        if (t->strand == STRAND_FWD) *e->cds = (Range) {e->cds->start, e->cds->end};
                        if (t->strand == STRAND_REV) *e->cds = (Range) {e->cds->end, e->cds->start};
                    }

                    kv_destroy(hap_in);
                }
            }

            for (int i = 0; i < records.n; i++) {
                Record *r = kv_A(records, i);

                // log_info("%d:%d", r->v->rid, r->v->pos);
                int splice_triggered = 0;
                for (int j = 0; j < r->vep->n; j++) {
                    if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) splice_triggered = 1;
                }
                if (!splice_triggered) continue;

                log_info("%s:%li (splice)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos + 1); 

                if (strlen(r->v->d.allele[1]) > CAI_DISTANCE) {
                    r->skipped = 1;
                    log_warn("Skipping on %s:%li, insertion is too long", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos + 1);
                    continue;
                }
                if (r->v->rlen > CAI_DISTANCE) {
                    r->skipped = 1;
                    log_warn("Skipping prediction on %s:%li, deletion is too long", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos + 1);
                    continue;
                }

                for (int g = 0; g < r->genes->n; g++) {
                    Gene *gene = kv_A(*r->genes, g);
                    if (gene == NULL) continue;

                    int n_unphased = 0;
                    for (int j = 0; j < records.n; j++) {
                        Record *r = kv_A(records, j);
                        if (has_matching_phase(r->v, hdr, phase, false)) {
                            continue;
                        }

                        r->unphased = 1;
                        n_unphased++;
                    }

                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;

                    // Find records for the current allele
                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    Records hap_ex; kv_init(hap_ex); // Variant for current allele excluding current variant

                    for (int j = 0; j < records.n; j++) {
                        Record *hap_r = kv_A(records, j);
                        if (hap_r->genes == NULL) continue;
                        if (hap_r->unphased & (hap_r != r)) continue;

                        Range record_range = variant_to_range(hap_r->v, 0);

                        for (int g2 = 0; g2 < hap_r->genes->n; g2++) {
                            Gene *gene2 = kv_A(*hap_r->genes, g2);

                            if (gene2 == NULL) continue;

                            if (gene->gff_id == gene2->gff_id) {
                                Range intersection;
                                intersect_ranges(gene->coords, record_range, &intersection);
                                if (range_len(intersection) == 0) continue;

                                kv_push(Record *, hap_in, hap_r);

                                if (hap_r == r) continue;
                                kv_push(Record *, hap_ex, hap_r);
                            }
                        }
                    }

                    Records just_alt;
                    kv_init(just_alt);
                    kv_push(Record *, just_alt, r);

                    kstring_t alt;
                    get_haplotype_sequence(gene->seq, &just_alt, gene->coords, &alt);

                    // Haplotype ref/alt sequences
                    kstring_t hap_ref;
                    get_haplotype_sequence(gene->seq, &hap_ex, gene->coords, &hap_ref);

                    kstring_t hap_alt;
                    get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_alt);


                    // REF
                    int32_t start = VARIANT_POS - (r->v->pos - gene->coords.start);
                    int32_t end = VARIANT_POS - (gene->coords.end - (r->v->pos+1)); // Open ended, so + 1

                    char *ref_seq = gene->seq.s - start; // Pointer to existing memory, no free necessary
                    if (start < 0) start = 0;
                    if (end < 0) end = 0;

                    add_phasing_for_prediction_records(&hap_in, (Range) {start, end});

                    // log_info("%s:%li (ref)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos); 
                    Range ref_boundary = { start, end };
                    char *ref_padded = pad_sequence(ref_seq, ref_boundary, WIDTH);
                    float *ref_predictions = NULL; int num_ref_predictions = 0;
                    int ref_len = strlen(ref_padded);
                    // log_info("ref_len: %i", ref_len);
                    predict_for_seq(models, ref_padded, ref_len, gene->strand, &ref_predictions, &num_ref_predictions);
                    // log_info("num_ref_predictions: %i", num_ref_predictions);
                    free(ref_padded);

                    // ALT
                    int alt_pred_offset = allele_offset(r->v, 0, 1, (Range) { r->v->pos - 50, r->v->pos + 51});

                    int32_t alt_start = VARIANT_POS - (r->v->pos - gene->coords.start);
                    int32_t alt_end = VARIANT_POS - ((gene->coords.end - alt_pred_offset) - r->v->pos);

                    char *alt_seq = alt.s - alt_start;
                    if (start < 0) start = 0;
                    if (end < 0) end = 0;

                    Range alt_boundary = { start, end };
                    char *alt_padded = pad_sequence(alt_seq, alt_boundary, WIDTH - alt_pred_offset);
                    float *alt_predictions = NULL; int num_alt_predictions = 0;
                    int alt_len = strlen(alt_padded);

                    predict_for_seq(models, alt_padded, alt_len, gene->strand, &alt_predictions, &num_alt_predictions);

                    float *corrected_alt_predictions = match_predictions_to_reference(alt_predictions, (Range) { r->v->pos - 50, r->v->pos + 51 }, &just_alt, num_ref_predictions);
                    free(alt_padded);

                    // HAP REF
                    // int pred_offset = prediction_offset(r->v, &hap_ex);
                    //
                    // int64_t gene_start = gene->coords.start;
                    // int64_t gene_end = gene->coords.end - haplotype_offset(&hap_ex, gene->coords);
                    // int64_t vposition = r->v->pos - haplotype_offset(&hap_ex, (Range) {gene->coords.start, r->v->pos});
                    //
                    // int start_offset = haplotype_offset(&hap_ex, (Range) { r->v->pos - CAI_DISTANCE, r->v->pos});
                    // int32_t hap_ref_start = VARIANT_POS + (gene_start - vposition) - start_offset;
                    // int32_t hap_ref_end = VARIANT_POS - (gene_end - vposition) - (pred_offset - start_offset);
                    //
                    // char *hap_ref_seq = hap_ref.s - hap_ref_start;
                    // if (hap_ref_start < 0) hap_ref_start = 0;
                    // if (hap_ref_end < 0) hap_ref_end = 0;
                    //
                    // // log_info("%s:%li (hap_ref)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos); 
                    // Range hap_ref_boundary = {hap_ref_start, hap_ref_end};
                    // char *hap_ref_padded = pad_sequence(hap_ref_seq, hap_ref_boundary, WIDTH - pred_offset);
                    // float *hap_ref_predictions = NULL; int num_hap_ref_predictions = 0;
                    // int hap_ref_len = strlen(hap_ref_padded);
                    // predict_for_seq(models, hap_ref_padded, hap_ref_len, gene->strand, &hap_ref_predictions, &num_hap_ref_predictions);
                    // float *corrected_hap_ref_predictions = match_predictions_to_reference(hap_ref_predictions, r->v, &hap_ex, num_ref_predictions);
                    //
                    // free(hap_ref_padded);

                    // HAP ALT
                    int pred_offset = prediction_offset(r->v, &hap_in);

                    int64_t gene_start = gene->coords.start;
                    int64_t gene_end = gene->coords.end - haplotype_offset(&hap_in, gene->coords);
                    int64_t vposition = r->v->pos - haplotype_offset(&hap_ex, (Range) {gene->coords.start, r->v->pos});

                    int64_t start_offset = haplotype_offset(&hap_in, (Range) { r->v->pos - CAI_DISTANCE, r->v->pos});
                    int32_t hap_alt_start = VARIANT_POS + (gene_start - vposition) - start_offset;
                    int32_t hap_alt_end = VARIANT_POS - (gene_end - vposition) - (pred_offset - start_offset);

                    char *hap_alt_seq = hap_alt.s - hap_alt_start;
                    if (hap_alt_start < 0) hap_alt_start = 0;
                    if (hap_alt_end < 0) hap_alt_end = 0;

                    // log_info("%s:%li (hap_alt)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos); 
                    Range hap_alt_boundary = {hap_alt_start, hap_alt_end};
                    char *hap_alt_padded = pad_sequence(hap_alt_seq, hap_alt_boundary, WIDTH - pred_offset);
                    float *hap_alt_predictions = NULL; int num_hap_alt_predictions = 0;
                    int hap_alt_len = strlen(hap_alt_padded);
                    predict_for_seq(models, hap_alt_padded, hap_alt_len, gene->strand, &hap_alt_predictions, &num_hap_alt_predictions);
                    float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, (Range) { r->v->pos - 50, r->v->pos + 51 }, &hap_in, num_ref_predictions);
                    free(hap_alt_padded);


                    // Calculate scores between ref and alt
                    Score score = calculate_delta_scores(r->v->d.allele[1], gene->name, ref_predictions, corrected_alt_predictions, num_ref_predictions, CAI_DISTANCE);
                    r->native_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
                    *r->native_scores[g * NUM_PHASES + phase] = score;


                    // Calculate scores between ref and hap_alt
                    Score hap_score = calculate_delta_scores(r->v->d.allele[1], gene->name, ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);
                    r->haplotype_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
                    *r->haplotype_scores[g * NUM_PHASES + phase] = hap_score;

                    // Score diff;
                    r->diff_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
                    Score diff = {0};
                    diff.alt = r->v->d.allele[1];
                    diff.gene = gene->name;
                    diff.ag_idx = hap_score.ag_idx;
                    diff.ag = corrected_hap_alt_predictions[(diff.ag_idx + 50) * NUM_SCORES + ACCEPTOR] - corrected_alt_predictions[(diff.ag_idx + 50) * NUM_SCORES + ACCEPTOR];

                    diff.al_idx = hap_score.al_idx;
                    diff.al = corrected_hap_alt_predictions[(diff.al_idx + 50) * NUM_SCORES + ACCEPTOR] - corrected_alt_predictions[(diff.al_idx + 50) * NUM_SCORES + ACCEPTOR];

                    diff.dg_idx = hap_score.dg_idx;
                    diff.dg = corrected_hap_alt_predictions[(diff.dg_idx + 50) * NUM_SCORES + DONOR] - corrected_alt_predictions[(diff.dg_idx + 50) * NUM_SCORES + DONOR];

                    diff.dl_idx = hap_score.dl_idx;
                    diff.dl = corrected_hap_alt_predictions[(diff.dl_idx + 50) * NUM_SCORES + DONOR] - corrected_alt_predictions[(diff.dl_idx + 50) * NUM_SCORES + DONOR];


                    *r->diff_scores[g * NUM_PHASES + phase] = diff;

                    if (
                        diff.ag > SPLICEAI_RESCUE_THRESHOLD |
                        diff.ag < -SPLICEAI_RESCUE_THRESHOLD |
                        diff.al > SPLICEAI_RESCUE_THRESHOLD |
                        diff.al < -SPLICEAI_RESCUE_THRESHOLD
                    ) {
                        for (int ti = 0; ti < r->transcript_genes->n; ti++) {
                            Gene *transcript_gene = kv_A(*r->transcript_genes, ti);

                            if (transcript_gene == NULL) continue;
                            if (transcript_gene->gff_id != gene->gff_id) continue;

                            set_record_flag(r, ti, phase, HAPLOTYPE_SPECIFIC_ACCEPTOR);
                        }
                    }

                    if (
                        diff.dg > SPLICEAI_RESCUE_THRESHOLD |
                        diff.dg < -SPLICEAI_RESCUE_THRESHOLD |
                        diff.dl > SPLICEAI_RESCUE_THRESHOLD |
                        diff.dl < -SPLICEAI_RESCUE_THRESHOLD
                    ) {
                        for (int ti = 0; ti < r->transcript_genes->n; ti++) {
                            Gene *transcript_gene = kv_A(*r->transcript_genes, ti);

                            if (transcript_gene == NULL) continue;
                            if (transcript_gene->gff_id != gene->gff_id) continue;

                            set_record_flag(r, ti, phase, HAPLOTYPE_SPECIFIC_DONOR);
                        }
                    }

                    const char *chr = bcf_seqname(hdr, r->v);

                    const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
                    const float *predictions[3] = { ref_predictions, corrected_alt_predictions, corrected_hap_alt_predictions };
                    write_predictions(scores_tsv, chr, r->v->pos, r->v->d.allele[0], r->v->d.allele[1], gene->name, ".", phase+1, "canonical", prediction_names, 3, predictions);


                // void write_predictions(const char *chr, const int64_t pos, const char *ref, const char *alt, const char *gene_name, const char *transcript_name, const int phase) {
                //     const char *chr = bcf_hdr_id2name(hdr, r->v->rid);
                //     for (int i = 0; i < 101; i++) {
                //         const int pos = i - 50;
                //
                //         const float ref_neither = ref_predictions[i * 3];
                //         const float ref_acceptor = ref_predictions[i * 3 + ACCEPTOR];
                //         const float ref_donor = ref_predictions[i * 3 + DONOR];
                //
                //         const float alt_neither = corrected_alt_predictions[i * 3];
                //         const float alt_acceptor = corrected_alt_predictions[i * 3 + ACCEPTOR];
                //         const float alt_donor = corrected_alt_predictions[i * 3 + DONOR];
                //
                //         const float hap_alt_neither = corrected_hap_alt_predictions[i * 3];
                //         const float hap_alt_acceptor = corrected_hap_alt_predictions[i * 3 + ACCEPTOR];
                //         const float hap_alt_donor = corrected_hap_alt_predictions[i * 3 + DONOR];
                //
                //         if (
                //             (ref_acceptor > SPLICEAI_RESCUE_THRESHOLD) | (ref_donor > SPLICEAI_RESCUE_THRESHOLD) |
                //             (alt_acceptor > SPLICEAI_RESCUE_THRESHOLD) | (alt_donor > SPLICEAI_RESCUE_THRESHOLD) | 
                //             (hap_alt_acceptor > SPLICEAI_RESCUE_THRESHOLD) | (hap_alt_donor > SPLICEAI_RESCUE_THRESHOLD)
                //         ) {
                //             fprintf(scores_tsv, "%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%.6f\t%.6f\t%.6f\n",
                //                     chr, r->v->pos + 1, r->v->d.allele[0], r->v->d.allele[1], gene->name, ".", phase+1, pos, "ref", ref_donor, ref_acceptor, ref_neither);
                //             fprintf(scores_tsv, "%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%.6f\t%.6f\t%.6f\n",
                //                     chr, r->v->pos + 1, r->v->d.allele[0], r->v->d.allele[1], gene->name, ".", phase+1, pos, "alt", alt_donor, alt_acceptor, alt_neither);
                //             fprintf(scores_tsv, "%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%.6f\t%.6f\t%.6f\n",
                //                     chr, r->v->pos + 1, r->v->d.allele[0], r->v->d.allele[1], gene->name, ".", phase+1, pos, "hap_alt", hap_alt_donor, hap_alt_acceptor, hap_alt_neither);
                //         }
                //
                //     }
                // }

                    // r->cs[g * NUM_PHASES + phase] = malloc(sizeof(SpliceChange));
                    // int ag_position = (score.ag_idx + CAI_DISTANCE) * NUM_SCORES + ACCEPTOR_POS;
                    // int al_position = (score.al_idx + CAI_DISTANCE) * NUM_SCORES + ACCEPTOR_POS;
                    // int dg_position = (score.dg_idx + CAI_DISTANCE) * NUM_SCORES + DONOR_POS;
                    // int dl_position = (score.dl_idx + CAI_DISTANCE) * NUM_SCORES + DONOR_POS;
                    // float cs_ag = corrected_hap_alt_predictions[ag_position] - corrected_hap_ref_predictions[ag_position];
                    // float cs_al = corrected_hap_alt_predictions[al_position] - corrected_hap_ref_predictions[al_position];
                    // float cs_dg = corrected_hap_alt_predictions[dg_position] - corrected_hap_ref_predictions[dg_position];
                    // float cs_dl = corrected_hap_alt_predictions[dl_position] - corrected_hap_ref_predictions[dl_position];
                    // *r->cs[g * NUM_PHASES + phase] = (SpliceChange) {cs_ag, cs_al, cs_dg, cs_dl, fabs(cs_ag) + fabs(cs_al) + fabs(cs_dg) + fabs(cs_dl) };

                    // If there's less than 0.1 change between ref and hap alt, don't predict for competitors and pairs
                    if (hap_score.ag < SPLICEAI_RESCUE_THRESHOLD && hap_score.al < SPLICEAI_RESCUE_THRESHOLD && hap_score.dg < SPLICEAI_RESCUE_THRESHOLD && hap_score.dl < SPLICEAI_RESCUE_THRESHOLD) continue;

                    // Remove predictions that are most likely to be neutral
                    filter_to_most_likely_prediction(num_ref_predictions, ref_predictions);
                    filter_to_most_likely_prediction(num_hap_alt_predictions, hap_alt_predictions);


                    int current_gene[r->tscripts->n];
                    for (int i = 0; i < r->tscripts->n; i++) {
                        Transcript *t = kv_A(*r->tscripts, i);
                        Gene *tgene = kv_A(*r->transcript_genes, i);

                        current_gene[i] = (t != NULL && tgene->gff_id == gene->gff_id);
                    }

                    Range pred_range = { r->v->pos - 50, r->v->pos + 51 };
                    SpliceProbabilities splice_probabilities;
                    kv_init(splice_probabilities);

                    // Get acceptor pairs
                    PairPrediction *ref_alt_acceptor_map = NULL;
                    size_t ref_alt_acceptor_map_size = 0;
                    get_prediction_pairs_for_site(ACCEPTOR_POS, ref_predictions, 101, hap_alt_predictions, 101-pred_offset, &hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
                    get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, &splice_probabilities);


                    // Get donor pairs
                    PairPrediction *ref_alt_donor_map = NULL;
                    size_t ref_alt_donor_map_size = 0;
                    get_prediction_pairs_for_site(DONOR_POS, ref_predictions, 101, hap_alt_predictions, 101-pred_offset, &hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
                    get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, &splice_probabilities);

                    for (int i = 0; i < r->tscripts->n; i++) {
                        Transcript *t = kv_A(*r->tscripts, i);
                        Gene *tgene = kv_A(*r->transcript_genes, i);

                        if (t == NULL || tgene->gff_id != gene->gff_id) continue;

                        int current_gene[r->tscripts->n];
                        for (int i = 0; i < r->tscripts->n; i++) current_gene[i] = 0;
                        current_gene[i] = 1;

                        for (int si = 0; si < splice_probabilities.n; si++) {
                            SpliceProbability prediction = splice_probabilities.a[si];
                            Range splice_position = (Range) { prediction.position, prediction.position+1 };
                            for (int ei = 0; ei < t->nexons; ei++) {
                                Exon *e = t->exons[ei];

                                int direction = 0;
                                if (t->strand == STRAND_FWD) direction = 1;
                                if (t->strand == STRAND_REV) direction = -1;
                                if (t->strand == STRAND_UNK) continue;

                                Range positions;
                                if (prediction.type == DONOR) {
                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                    positions = (Range) { e->hap_bounds->start + 1, t->exons[ei + direction]->hap_bounds->start };
                                } else { // ACCEPTOR
                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                    positions = (Range) { t->exons[ei - direction]->hap_bounds->end, e->hap_bounds->end - 1 };
                                }

                                Range intersection;
                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                if (prediction.type == DONOR) {
                                    prediction.canonical_position = e->splice_sites->canonical_donors[0].position;
                                } else { // ACCEPTOR
                                    prediction.canonical_position = e->splice_sites->canonical_acceptors[0].position;
                                }

                                int predict_for_pair = 0;
                                int predict_for_competitor = 0;
                                store_splice_sites(prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);

                                if (predict_for_pair) {
                                    if (prediction.type == DONOR) {
                                        int pi = ei + direction;
                                        if (pi < 0 || pi >= t->nexons) continue;
                                        SpliceSite pair = t->exons[pi]->splice_sites->canonical_acceptors[0];

                                        kv_push_new_record((Records *) pair.records_ptr, prediction.record);

                                        if (pair.predicted) continue;

                                        Range pred_range = { pair.position - 50, pair.position + 51 };
                                        float *pair_predictions_ref = predict_on_sequence(models, gene->seq.s, WIDTH, pair.position, tgene->coords, tgene->strand, NULL);

                                        int pred_offset = haplotype_offset(&just_alt, pred_range);
                                        float *pair_predictions_alt = predict_on_sequence(models, alt.s, WIDTH-pred_offset, pair.position, tgene->coords, tgene->strand, &just_alt);
                                        float *corrected_pair_predictions_alt = match_predictions_to_reference(pair_predictions_alt, pred_range, &just_alt, 303);

                                        int hap_pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *pair_predictions_hap_alt = predict_on_sequence(models, hap_alt.s, WIDTH-hap_pred_offset, pair.position, tgene->coords, tgene->strand, &hap_in);
                                        float *corrected_pair_predictions_hap_alt = match_predictions_to_reference(pair_predictions_hap_alt, pred_range, &hap_in, 303);

                                        const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
                                        const float *predictions[3] = { ref_predictions, corrected_pair_predictions_alt, corrected_pair_predictions_hap_alt};
                                        write_predictions(scores_tsv, chr, pair.position, r->v->d.allele[0], r->v->d.allele[1], gene->name, gff_id2string(gff, 0, t->id), phase+1, "pair", prediction_names, 3, predictions);


                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, pair_predictions_ref);
                                        filter_to_most_likely_prediction((PRED_WIDTH-hap_pred_offset) * NUM_SCORES, pair_predictions_hap_alt);


                                        SpliceProbabilities pair_probabilities;
                                        kv_init(pair_probabilities);

                                        // Get acceptor pairs
                                        PairPrediction *ref_alt_acceptor_map = NULL;
                                        size_t ref_alt_acceptor_map_size = 0;
                                        get_prediction_pairs_for_site(ACCEPTOR_POS, pair_predictions_ref, 101, pair_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
                                        get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, &pair_probabilities);


                                        // Get donor pairs
                                        PairPrediction *ref_alt_donor_map = NULL;
                                        size_t ref_alt_donor_map_size = 0;
                                        get_prediction_pairs_for_site(DONOR_POS, pair_predictions_ref, 101, pair_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
                                        get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, &pair_probabilities);


                                        for (int pi = 0; pi < pair_probabilities.n; pi++) {
                                            SpliceProbability pair_prediction = pair_probabilities.a[pi];
                                            Range pair_splice_position = (Range) { pair_prediction.position, pair_prediction.position+1 };

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (pair_prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->hap_bounds->start + 1, t->exons[ei + direction]->hap_bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->hap_bounds->end, e->hap_bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(pair_splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                if (pair_prediction.type == DONOR) {
                                                    pair_prediction.canonical_position = e->splice_sites->canonical_donors[0].position;
                                                } else { // ACCEPTOR
                                                    pair_prediction.canonical_position = e->splice_sites->canonical_acceptors[0].position;
                                                }

                                                store_splice_sites(pair_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    } else { // ACCEPTOR
                                        int pi = ei - direction;
                                        if (pi < 0 || pi >= t->nexons) continue;
                                        SpliceSite pair = t->exons[pi]->splice_sites->canonical_donors[0];

                                        kv_push_new_record((Records *) pair.records_ptr, prediction.record);

                                        if (pair.predicted) continue;

                                        Range pred_range = { pair.position - 50, pair.position + 51 };
                                        float *pair_predictions_ref = predict_on_sequence(models, gene->seq.s, WIDTH, pair.position, tgene->coords, tgene->strand, NULL);
                                        // float *pair_predictions_alt = predict_on_sequence(models, hap_alt.s, WIDTH - pred_offset, pair.position, tgene->coords, tgene->strand, &hap_in);

                                        int pred_offset = haplotype_offset(&just_alt, pred_range);
                                        float *pair_predictions_alt = predict_on_sequence(models, alt.s, WIDTH-pred_offset, pair.position, tgene->coords, tgene->strand, &just_alt);
                                        float *corrected_pair_predictions_alt = match_predictions_to_reference(pair_predictions_alt, pred_range, &just_alt, 303);

                                        int hap_pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *pair_predictions_hap_alt = predict_on_sequence(models, hap_alt.s, WIDTH-hap_pred_offset, pair.position, tgene->coords, tgene->strand, &hap_in);
                                        float *corrected_pair_predictions_hap_alt = match_predictions_to_reference(pair_predictions_hap_alt, pred_range, &hap_in, 303);

                                        const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
                                        const float *predictions[3] = { ref_predictions, corrected_pair_predictions_alt, corrected_pair_predictions_hap_alt};
                                        write_predictions(scores_tsv, chr, pair.position, r->v->d.allele[0], r->v->d.allele[1], gene->name, gff_id2string(gff, 0, t->id), phase+1, "pair", prediction_names, 3, predictions);

                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, pair_predictions_ref);
                                        filter_to_most_likely_prediction((PRED_WIDTH - hap_pred_offset) * NUM_SCORES, pair_predictions_hap_alt);


                                        SpliceProbabilities pair_probabilities;
                                        kv_init(pair_probabilities);

                                        // Get acceptor pairs
                                        PairPrediction *ref_alt_acceptor_map = NULL;
                                        size_t ref_alt_acceptor_map_size = 0;
                                        get_prediction_pairs_for_site(ACCEPTOR_POS, pair_predictions_ref, 101, pair_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
                                        get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, &pair_probabilities);


                                        // Get donor pairs
                                        PairPrediction *ref_alt_donor_map = NULL;
                                        size_t ref_alt_donor_map_size = 0;
                                        get_prediction_pairs_for_site(DONOR_POS, pair_predictions_ref, 101, pair_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
                                        get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, &pair_probabilities);

                                        for (int pi = 0; pi < pair_probabilities.n; pi++) {
                                            SpliceProbability pair_prediction = pair_probabilities.a[pi];
                                            Range pair_splice_position = (Range) { pair_prediction.position, pair_prediction.position+1 };

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (pair_prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->hap_bounds->start + 1, t->exons[ei + direction]->hap_bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->hap_bounds->end, e->hap_bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(pair_splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                if (pair_prediction.type == DONOR) {
                                                    pair_prediction.canonical_position = e->splice_sites->canonical_donors[0].position;
                                                } else { // ACCEPTOR
                                                    pair_prediction.canonical_position = e->splice_sites->canonical_acceptors[0].position;
                                                }

                                                store_splice_sites(pair_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    }
                                }

                                if (predict_for_competitor) {
                                    if (prediction.type == DONOR) {
                                        int ci = ei - direction;
                                        if (ci < 0 || ci >= t->nexons) continue;
                                        SpliceSite competitor = t->exons[ci]->splice_sites->canonical_donors[0];

                                        kv_push_new_record((Records *) competitor.records_ptr, prediction.record);

                                        if (competitor.predicted) continue;

                                        Range pred_range = { competitor.position - 50, competitor.position + 51 };
                                        // int pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *competitor_predictions_ref = predict_on_sequence(models, gene->seq.s, WIDTH, competitor.position, tgene->coords, tgene->strand, NULL);
                                        // float *competitor_predictions_alt = predict_on_sequence(models, hap_alt.s, WIDTH - pred_offset, competitor.position, tgene->coords, tgene->strand, &hap_in);

                                        int pred_offset = haplotype_offset(&just_alt, pred_range);
                                        float *competitor_predictions_alt = predict_on_sequence(models, alt.s, WIDTH-pred_offset, competitor.position, tgene->coords, tgene->strand, &just_alt);
                                        float *corrected_competitor_predictions_alt = match_predictions_to_reference(competitor_predictions_alt, pred_range, &just_alt, 303);

                                        int hap_pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *competitor_predictions_hap_alt = predict_on_sequence(models, hap_alt.s, WIDTH-hap_pred_offset, competitor.position, tgene->coords, tgene->strand, &hap_in);
                                        float *corrected_competitor_predictions_hap_alt = match_predictions_to_reference(competitor_predictions_hap_alt, pred_range, &hap_in, 303);

                                        const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
                                        const float *predictions[3] = { ref_predictions, corrected_competitor_predictions_alt, corrected_competitor_predictions_hap_alt};
                                        write_predictions(scores_tsv, chr, competitor.position, r->v->d.allele[0], r->v->d.allele[1], gene->name, gff_id2string(gff, 0, t->id), phase+1, "competitor", prediction_names, 3, predictions);

                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, competitor_predictions_ref);
                                        filter_to_most_likely_prediction((PRED_WIDTH - hap_pred_offset) * NUM_SCORES, competitor_predictions_hap_alt);


                                        SpliceProbabilities competitor_probabilities;
                                        kv_init(competitor_probabilities);

                                        // Get acceptor pairs
                                        PairPrediction *ref_alt_acceptor_map = NULL;
                                        size_t ref_alt_acceptor_map_size = 0;
                                        get_prediction_pairs_for_site(ACCEPTOR_POS, competitor_predictions_ref, 101, competitor_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
                                        get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, &competitor_probabilities);


                                        // Get donor pairs
                                        PairPrediction *ref_alt_donor_map = NULL;
                                        size_t ref_alt_donor_map_size = 0;
                                        get_prediction_pairs_for_site(DONOR_POS, competitor_predictions_ref, 101, competitor_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
                                        get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, &competitor_probabilities);

                                        for (int csi = 0; csi < competitor_probabilities.n; csi++) {
                                            SpliceProbability competitor_prediction = competitor_probabilities.a[csi];
                                            Range competitor_splice_position = (Range) { competitor_prediction.position, competitor_prediction.position+1 };

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (competitor_prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->hap_bounds->start + 1, t->exons[ei + direction]->hap_bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->hap_bounds->end, e->hap_bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(competitor_splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                if (competitor_prediction.type == DONOR) {
                                                    competitor_prediction.canonical_position = e->splice_sites->canonical_donors[0].position;
                                                } else { // ACCEPTOR
                                                    competitor_prediction.canonical_position = e->splice_sites->canonical_acceptors[0].position;
                                                }

                                                store_splice_sites(competitor_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    } else { // ACCEPTOR
                                        int ci = ei + direction;
                                        if (ci < 0 || ci >= t->nexons) continue;
                                        SpliceSite competitor = t->exons[ci]->splice_sites->canonical_acceptors[0];

                                        kv_push_new_record((Records *) competitor.records_ptr, prediction.record);

                                        if (competitor.predicted) continue;

                                        Range pred_range = { competitor.position - 50, competitor.position + 51 };
                                        // int pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *competitor_predictions_ref = predict_on_sequence(models, gene->seq.s, WIDTH, competitor.position, tgene->coords, tgene->strand, NULL);
                                        // float *competitor_predictions_alt = predict_on_sequence(models, hap_alt.s, WIDTH - pred_offset, competitor.position, tgene->coords, tgene->strand, &hap_in);
                                        int pred_offset = haplotype_offset(&just_alt, pred_range);
                                        float *competitor_predictions_alt = predict_on_sequence(models, alt.s, WIDTH-pred_offset, competitor.position, tgene->coords, tgene->strand, &just_alt);
                                        float *corrected_competitor_predictions_alt = match_predictions_to_reference(competitor_predictions_alt, pred_range, &just_alt, 303);

                                        int hap_pred_offset = haplotype_offset(&hap_in, pred_range);
                                        float *competitor_predictions_hap_alt = predict_on_sequence(models, hap_alt.s, WIDTH-hap_pred_offset, competitor.position, tgene->coords, tgene->strand, &hap_in);
                                        float *corrected_competitor_predictions_hap_alt = match_predictions_to_reference(competitor_predictions_hap_alt, pred_range, &hap_in, 303);

                                        const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
                                        const float *predictions[3] = { ref_predictions, corrected_competitor_predictions_alt, corrected_competitor_predictions_hap_alt};
                                        write_predictions(scores_tsv, chr, competitor.position, r->v->d.allele[0], r->v->d.allele[1], gene->name, gff_id2string(gff, 0, t->id), phase+1, "competitor", prediction_names, 3, predictions);

                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, competitor_predictions_ref);
                                        filter_to_most_likely_prediction((PRED_WIDTH - hap_pred_offset) * NUM_SCORES, competitor_predictions_hap_alt);


                                        SpliceProbabilities competitor_probabilities;
                                        kv_init(competitor_probabilities);

                                        // Get acceptor pairs
                                        PairPrediction *ref_alt_acceptor_map = NULL;
                                        size_t ref_alt_acceptor_map_size = 0;
                                        get_prediction_pairs_for_site(ACCEPTOR_POS, competitor_predictions_ref, 101, competitor_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
                                        get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, &competitor_probabilities);


                                        // Get donor pairs
                                        PairPrediction *ref_alt_donor_map = NULL;
                                        size_t ref_alt_donor_map_size = 0;
                                        get_prediction_pairs_for_site(DONOR_POS, competitor_predictions_ref, 101, competitor_predictions_hap_alt, 101-hap_pred_offset, &hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
                                        get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, &competitor_probabilities);

                                        for (int ci = 0; ci < competitor_probabilities.n; ci++) {
                                            SpliceProbability competitor_prediction = competitor_probabilities.a[ci];
                                            Range competitor_splice_position = (Range) { competitor_prediction.position, competitor_prediction.position+1 };

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (competitor_prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->hap_bounds->start + 1, t->exons[ei + direction]->hap_bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->hap_bounds->end, e->hap_bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(competitor_splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                if (competitor_prediction.type == DONOR) {
                                                    competitor_prediction.canonical_position = e->splice_sites->canonical_donors[0].position;
                                                } else { // ACCEPTOR
                                                    competitor_prediction.canonical_position = e->splice_sites->canonical_acceptors[0].position;
                                                }

                                                store_splice_sites(competitor_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    free(ref_predictions);
                    free(alt_predictions);
                    free(corrected_alt_predictions);
                    // free(hap_ref_predictions);
                    // free(corrected_hap_ref_predictions);
                    free(hap_alt_predictions);
                    free(corrected_hap_alt_predictions);

                    free(hap_ref.s);
                    free(hap_alt.s);

                    kv_destroy(hap_in);
                    kv_destroy(hap_ex);
                }
            }

            /**
            *** ISOFORM FORMATION
            **/

            for (int i = 0; i < records.n; i++) {
                Record *r = kv_A(records, i);

                if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;


                for (int ti = 0; ti < r->tscripts->n; ti++) {
                    Transcript *t = kv_A(*r->tscripts, ti);
                    Gene *gene = kv_A(*r->transcript_genes, ti);
                    TscriptCsqs *csq = kv_A(*r->vep, ti);

                    if (t == NULL || gene == NULL) continue;
                    if (!csq->triggers[SPLICE_TRIGGER] && !csq->triggers[FRAMESHIFT_TRIGGER] && !csq->triggers[STOP_TRIGGER]) continue;

                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;


                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    for (int j = 0; j < records.n; j++) {
                        Record *hap_r = kv_A(records, j);
                        if (hap_r->genes == NULL) continue;
                        Range record_range = variant_to_range(hap_r->v, 0);

                        for (int gi2 = 0; gi2 < hap_r->genes->n; gi2++) {
                            Gene *gene2 = kv_A(*hap_r->genes, gi2);
                            if (gene2 == NULL) continue;

                            if (gene->gff_id == gene2->gff_id) {
                                Range intersection;
                                intersect_ranges(gene->coords, record_range, &intersection);
                                if (range_len(intersection) == 0) continue;

                                kv_push(Record *, hap_in, hap_r);

                            }
                        }
                    }

                    // Initialize 2d array for storing noncoding regions
                    int n_nc_regions = 0, m_nc_regions = t->nexons+10; // TODO: rework this, +10 is messy, and probably not needed in the vast majority of cases
                    NonCodingRegions *nc_regions = malloc(m_nc_regions * sizeof(NonCodingRegions));
                    for (int i2 = 0; i2 < m_nc_regions; i2++) {
                        kv_init(nc_regions[i2]);
                    }

                    int is_coding = 0;
                    int first_cds, last_cds;
                    for (int ei = 0; ei < t->nexons; ei++) {
                        if (t->exons[ei]->cds == NULL) continue;

                        is_coding = 1;
                        first_cds = ei;
                        break;
                    }
                    for (int ei = t->nexons-1; ei >= 0; ei--) {
                        if (t->exons[ei]->cds == NULL) continue;

                        last_cds = ei;
                        break;
                    }

                    int gene_end = gene->coords.end;// - haplotype_offset(&hap_in, gene->coords);

                    if (is_coding && t->strand == STRAND_FWD) {
                        for (int ei = first_cds; ei <= last_cds; ei++) {
                            Exon *e = t->exons[ei];

                            if (ei == first_cds) { // 5' UTR
                                NonCodingRegion ncr = {
                                    1,
                                    { gene->coords.start, e->cds->start }, 
                                    { gene->coords.start, e->cds->start },
                                    NULL,
                                    NULL,
                                    UTR5,
                                    0
                                };
                                kv_init(ncr.records);
                                kv_push(Record *, ncr.records, r);
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                            }

                            if (ei == last_cds) { // 3' UTR
                                NonCodingRegion ncr = {
                                    1,
                                    { e->cds->end, gene_end },
                                    { e->cds->end, gene_end },
                                    NULL,
                                    NULL,
                                    UTR3,
                                    0
                                };
                                kv_init(ncr.records);
                                kv_push(Record *, ncr.records, r);
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                break; // End of transcript
                            }

                            if (ei+1 > last_cds) continue;

                            SpliceSite *canonical_donor = t->exons[ei]->splice_sites->canonical_donors;
                            uint32_t lost_donor = canonical_donor->predicted && canonical_donor->posterior_probability == 0;
                            uint32_t ablated_donor = canonical_donor->predicted && (canonical_donor->ref_prior_probability - canonical_donor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;

                            SpliceSite *canonical_acceptor = t->exons[ei+1]->splice_sites->canonical_acceptors;
                            uint32_t lost_acceptor = canonical_acceptor->predicted && canonical_acceptor->posterior_probability == 0;
                            uint32_t ablated_acceptor = canonical_acceptor->predicted && (canonical_acceptor->ref_prior_probability - canonical_acceptor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;

                            if (lost_donor && lost_acceptor) {
                                NonCodingRegion ncr = {
                                    0,
                                    { canonical_donor->position, canonical_donor->position },
                                    { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                    canonical_donor,
                                    canonical_acceptor,
                                    0,
                                    0
                                }; // No splicing happens for this donor
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);

                                n_nc_regions++;
                                continue;
                            }

                            kvec_t(SpliceSite) donors; kv_init(donors);
                            SpliceSites target_donor_sites = *t->exons[ei]->splice_sites;
                            for (int i = 0; i < target_donor_sites.n_canonical_donors; i++) {
                                SpliceSite *donor = target_donor_sites.canonical_donors + i;
                                if (donor->posterior_probability == 0) continue;
                                kv_push(SpliceSite, donors, *donor);
                            }
                            for (int i = 0; i < target_donor_sites.n_cryptic_donors; i++) {
                                SpliceSite *donor = target_donor_sites.cryptic_donors + i;
                                if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;
                                kv_push(SpliceSite, donors, *donor);
                            }
                            if (ablated_donor && ei-1 >= first_cds) {
                                SpliceSites alternative_donor_sites = *t->exons[ei-1]->splice_sites;
                                for (int i = 0; i < alternative_donor_sites.n_canonical_donors; i++) {
                                    SpliceSite *donor = alternative_donor_sites.canonical_donors + i;
                                    if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;
                                    kv_push(SpliceSite, donors, *donor);
                                }
                                for (int i = 0; i < alternative_donor_sites.n_cryptic_donors; i++) {
                                    SpliceSite *donor = alternative_donor_sites.cryptic_donors + i;
                                    if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;
                                    kv_push(SpliceSite, donors, *donor);
                                }
                            }

                            kvec_t(SpliceSite) acceptors; kv_init(acceptors);
                            SpliceSites target_acceptor_sites = *t->exons[ei+1]->splice_sites;
                            for (int i = 0; i < target_acceptor_sites.n_canonical_acceptors; i++) {
                                SpliceSite *acceptor = target_acceptor_sites.canonical_acceptors + i;
                                if (acceptor->posterior_probability == 0) continue;

                                kv_push(SpliceSite, acceptors, *acceptor);
                            }
                            for (int i = 0; i < target_acceptor_sites.n_cryptic_acceptors; i++) {
                                SpliceSite *acceptor = target_acceptor_sites.cryptic_acceptors + i;
                                if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                kv_push(SpliceSite, acceptors, *acceptor);
                            }
                            if (ablated_acceptor && ei+2 < last_cds) {
                                SpliceSites alternative_acceptor_sites = *t->exons[ei+2]->splice_sites;
                                for (int i = 0; i < alternative_acceptor_sites.n_canonical_acceptors; i++) {
                                    SpliceSite *acceptor = alternative_acceptor_sites.canonical_acceptors + i;
                                    if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                    kv_push(SpliceSite, acceptors, *acceptor);
                                }
                                for (int i = 0; i < alternative_acceptor_sites.n_cryptic_acceptors; i++) {
                                    SpliceSite *acceptor = alternative_acceptor_sites.cryptic_acceptors + i;
                                    if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                    kv_push(SpliceSite, acceptors, *acceptor);
                                }
                            }

                            // If no pairs have been found, add as empty
                            if (donors.n == 0 || acceptors.n == 0) {
                                for (int i = 0; i < donors.n; i++) {
                                    SpliceSite *donor = donors.a + i;
                                    NonCodingRegion ncr = {
                                        0,
                                        { donor->position, donor->position },
                                        { donor->canonical_position, canonical_acceptor->canonical_position },
                                        donor,
                                        canonical_acceptor,
                                        0,
                                        0
                                    }; // No splicing happens for this donor
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }

                                for (int i = 0; i < acceptors.n; i++) {
                                    SpliceSite *acceptor = acceptors.a + i;
                                    NonCodingRegion ncr = {
                                        0,
                                        { canonical_donor->position, canonical_donor->position },
                                        { canonical_donor->canonical_position, acceptor->canonical_position },
                                        canonical_donor,
                                        acceptor,
                                        0,
                                        0
                                    }; // No splicing happens for this donor
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }

                                n_nc_regions++;
                                continue;
                            }

                            for (int i = 0; i < donors.n; i++) {
                                SpliceSite *donor = donors.a + i;
                                for (int j = 0; j < acceptors.n; j++) {
                                    SpliceSite *acceptor = acceptors.a + j;
                                    NonCodingRegion ncr = {
                                        1,
                                        { donor->position, acceptor->position },
                                        { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                        donor,
                                        acceptor,
                                        0,
                                        canonical_donor->canonical_position != donor->canonical_position || canonical_acceptor->canonical_position != acceptor->canonical_position
                                    };
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            n_nc_regions++;
                            continue;
                        }
                    } else if (is_coding && t->strand == STRAND_REV) {
                        for (int ei = last_cds; ei >= first_cds; ei--) {
                            Exon *e = t->exons[ei];

                            // Add UTRs
                            if (ei == last_cds) { // 5' UTR
                                NonCodingRegion ncr = {
                                    1,
                                    { gene_end, e->cds->start },
                                    { gene_end, e->cds->start },
                                    NULL,
                                    NULL,
                                    UTR5,
                                    0
                                };
                                kv_init(ncr.records);
                                kv_push(Record *, ncr.records, r);
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                            }
                            if (ei == first_cds) { // 3' UTR
                                NonCodingRegion ncr = {
                                    1,
                                    { e->cds->end,  gene->coords.start },
                                    { e->cds->end, gene->coords.start },
                                    NULL,
                                    NULL,
                                    UTR3,
                                    0
                                };
                                kv_init(ncr.records);
                                kv_push(Record *, ncr.records, r);
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                break; // End of transcript
                            }

                            if (ei-1 < first_cds) continue;

                            SpliceSite *canonical_donor = t->exons[ei]->splice_sites->canonical_donors;
                            uint32_t lost_donor = canonical_donor->predicted && canonical_donor->posterior_probability == 0;
                            uint32_t ablated_donor = canonical_donor->predicted && (canonical_donor->ref_prior_probability - canonical_donor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;

                            SpliceSite *canonical_acceptor = t->exons[ei-1]->splice_sites->canonical_acceptors;
                            uint32_t lost_acceptor = canonical_acceptor->predicted && canonical_acceptor->posterior_probability == 0;
                            uint32_t ablated_acceptor = canonical_acceptor->predicted && (canonical_acceptor->ref_prior_probability - canonical_acceptor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;



                            // If there are no valid splicing pairs, in a +- 15 bp window around the canonicals
                            if (lost_donor && lost_acceptor) {
                                NonCodingRegion ncr = {
                                    0,
                                    { canonical_donor->position, canonical_donor->position },
                                    { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                    canonical_donor,
                                    canonical_acceptor,
                                    0,
                                    0
                                }; // No splicing happens for this donor
                                kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);

                                n_nc_regions++;
                                continue;
                            }

                            kvec_t(SpliceSite) donors; kv_init(donors);
                            SpliceSites target_donor_sites = *t->exons[ei]->splice_sites;
                            for (int i = 0; i < target_donor_sites.n_canonical_donors; i++) {
                                SpliceSite *donor = target_donor_sites.canonical_donors + i;
                                if (donor->posterior_probability == 0) continue;
                                kv_push(SpliceSite, donors, *donor);
                            }
                            for (int i = 0; i < target_donor_sites.n_cryptic_donors; i++) {
                                SpliceSite *donor = target_donor_sites.cryptic_donors + i;
                                if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;
                                kv_push(SpliceSite, donors, *donor);
                            }
                            if (ablated_donor && ei+1 < last_cds-1) {
                                SpliceSites alternative_donor_sites = *t->exons[ei+1]->splice_sites;
                                for (int i = 0; i < alternative_donor_sites.n_canonical_donors; i++) {
                                    SpliceSite *donor = alternative_donor_sites.canonical_donors + i;
                                    if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;

                                    kv_push(SpliceSite, donors, *donor);
                                }
                                for (int i = 0; i < alternative_donor_sites.n_cryptic_donors; i++) {
                                    SpliceSite *donor = alternative_donor_sites.cryptic_donors + i;
                                    if (donor->posterior_probability == 0 || (labs((int64_t) donor->position - (int64_t) canonical_donor->position)) > 500) continue;
                                    kv_push(SpliceSite, donors, *donor);
                                }
                            }

                            kvec_t(SpliceSite) acceptors; kv_init(acceptors);
                            SpliceSites target_acceptor_sites = *t->exons[ei-1]->splice_sites;
                            for (int i = 0; i < target_acceptor_sites.n_canonical_acceptors; i++) {
                                SpliceSite *acceptor = target_acceptor_sites.canonical_acceptors + i;
                                if (acceptor->posterior_probability == 0) continue;
                                kv_push(SpliceSite, acceptors, *acceptor);
                            }
                            for (int i = 0; i < target_acceptor_sites.n_cryptic_acceptors; i++) {
                                SpliceSite *acceptor = target_acceptor_sites.cryptic_acceptors + i;
                                if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                kv_push(SpliceSite, acceptors, *acceptor);
                            }
                            if (ablated_acceptor && ei-2 > first_cds) {
                                SpliceSites alternative_acceptor_sites = *t->exons[ei-2]->splice_sites;
                                for (int i = 0; i < alternative_acceptor_sites.n_canonical_acceptors; i++) {
                                    SpliceSite *acceptor = alternative_acceptor_sites.canonical_acceptors + i;
                                    if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                    kv_push(SpliceSite, acceptors, *acceptor);
                                }
                                for (int i = 0; i < alternative_acceptor_sites.n_cryptic_acceptors; i++) {
                                    SpliceSite *acceptor = alternative_acceptor_sites.cryptic_acceptors + i;
                                    if (acceptor->posterior_probability == 0 || (labs((int64_t) acceptor->position - (int64_t) canonical_acceptor->position)) > 500) continue;
                                    kv_push(SpliceSite, acceptors, *acceptor);
                                }
                            }

                            // If no pairs have been found, add as empty
                            if (donors.n == 0 || acceptors.n == 0) {
                                for (int i = 0; i < donors.n; i++) {
                                    SpliceSite *donor = donors.a + i;
                                    NonCodingRegion ncr = {
                                        0,
                                        { donor->position, donor->position },
                                        { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                        donor,
                                        canonical_acceptor,
                                        0,
                                        0
                                    }; // No splicing happens for this donor
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }

                                for (int i = 0; i < acceptors.n; i++) {
                                    SpliceSite *acceptor = acceptors.a + i;
                                    NonCodingRegion ncr = {
                                        0,
                                        { canonical_donor->position, canonical_donor->position },
                                        { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                        canonical_donor,
                                        acceptor,
                                        0,
                                        0
                                    }; // No splicing happens for this donor
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }

                                n_nc_regions++;
                                continue;
                            }

                            for (int i = 0; i < donors.n; i++) {
                                SpliceSite *donor = donors.a + i;
                                for (int j = 0; j < acceptors.n; j++) {
                                    SpliceSite *acceptor = acceptors.a + j;
                                    NonCodingRegion ncr = {
                                        1,
                                        { donor->position, acceptor->position },
                                        { canonical_donor->canonical_position, canonical_acceptor->canonical_position },
                                        donor,
                                        acceptor,
                                        0,
                                        canonical_donor->canonical_position != donor->canonical_position || canonical_acceptor->canonical_position != acceptor->canonical_position
                                    };
                                    kv_push(NonCodingRegion, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            n_nc_regions++;
                            continue;
                        }
                    }

                    if (!is_coding) continue; // Transcript is not protein-coding after all!

                    int nisoforms = 0, misoforms = 1;
                    for (int i = 0; i < n_nc_regions; i++) {
                        NonCodingRegions regions = nc_regions[i];
                        if (regions.n == 0) continue;

                        misoforms *= regions.n;
                    }

                    NonCodingRegions **intron_sets = malloc(misoforms * sizeof(NonCodingRegions *));
                    for (int i = 0; i < misoforms; i++) {
                        intron_sets[i] = malloc(sizeof(NonCodingRegions));
                        kv_init(*intron_sets[i]);
                    }

                    for (int i = 0; i < n_nc_regions; i++) {
                        NonCodingRegions ranges = nc_regions[i];
                        if (ranges.n == 0) continue;

                        if (nisoforms == 0) { // First isoform
                            for (int j = 0; j < ranges.n; j++) {
                                NonCodingRegions *first = intron_sets[j];
                                NonCodingRegion to_push = kv_A(ranges, j);
                                kv_push(NonCodingRegion, *first, to_push);
                                nisoforms++;
                            }
                        } else { // Extend/duplicate isoforms already created
                            // Copy over isoforms to add new endings
                            int prior_n_isoforms = nisoforms;
                            for (int j = 1; j < ranges.n; j++) {
                                for (int k = 0; k < prior_n_isoforms; k++) {
                                    NonCodingRegions *original = intron_sets[k];
                                    NonCodingRegions *new = intron_sets[j * prior_n_isoforms + k];
                                    kv_copy(NonCodingRegion, *new, *original);
                                    nisoforms++;
                                }
                            }

                            for (int j = 0; j < ranges.n; j++) {
                                for (int k = 0; k < prior_n_isoforms; k++) {
                                    NonCodingRegions *target = intron_sets[j * prior_n_isoforms + k];
                                    NonCodingRegion to_push = kv_A(ranges, j);
                                    kv_push(NonCodingRegion, *target, to_push);
                                }
                            }
                        }
                    }

                    kstring_t hap_seq;
                    get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_seq);
                    log_info("%s:%li (num_isoforms=%i)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos, nisoforms); 

                    Isoforms isoforms = {
                        1,
                        0,
                        nisoforms,
                        0,
                        1,
                        0,
                        0.0,
                        0,
                        0.0,
                        NO_PTV,
                        malloc(nisoforms * sizeof(Isoform *)),
                        NULL,
                        NULL
                    };

                    // TODO: On rework the order of this should be changed, ALL splice effects should be resolved before haplotype effects get applied.
                    // Right now the back and forth on the coordinates makes a mess of things.
                    for (int i = 0; i < nisoforms; i++) {
                        NonCodingRegions *ncrs = intron_sets[i];
                        kstring_t spliced_seq = {0};
                        int seq_idx = 0;
                        log_debug("Processing isoform: %i", i);
                        log_debug("Gene length is: %i", hap_seq.l);

                        NonCodingRegions hap_ncrs;
                        kv_init(hap_ncrs);

                        if (t->strand == STRAND_FWD) {
                            for (int i = 0; i < ncrs->n; i++) {
                                NonCodingRegion ncr = kv_A(*ncrs, i);

                                ncr.coords.start -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.coords.start + 1 });
                                ncr.coords.end -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.coords.end });

                                ncr.canonical_coords.start -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.canonical_coords.start + 1});
                                ncr.canonical_coords.end -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.canonical_coords.end });

                                kv_push(NonCodingRegion, hap_ncrs, ncr);
                            }
                        } else {
                            for (int i = 0; i < ncrs->n; i++) {
                                NonCodingRegion ncr = kv_A(*ncrs, i);

                                ncr.coords.start -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.coords.start - 1 });
                                ncr.coords.end -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.coords.end });

                                ncr.canonical_coords.start -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.canonical_coords.start - 1});
                                ncr.canonical_coords.end -= haplotype_offset(&hap_in, (Range) { gene->coords.start, ncr.canonical_coords.end });

                                kv_push(NonCodingRegion, hap_ncrs, ncr);
                            }
                        }

                        Range *isoform_bounds = NULL;

                        int gene_start = gene->coords.start;
                        if (t->strand == STRAND_FWD) {
                            int length = 0;
                            seq_idx = 1;
                            for (int j = 0; j < hap_ncrs.n; j++) {
                                NonCodingRegion ncr = kv_A(hap_ncrs, j);
                                int intron_start = ncr.coords.start - gene_start + 1;
                                int intron_end = ncr.coords.end - gene_start;

                                if (intron_start < seq_idx) {
                                    continue;
                                }

                                int bases_shifted = range_len(ncr.canonical_coords) - range_len(ncr.coords);
                                length += intron_start - seq_idx;

                                if (length > 0 && isoform_bounds == NULL) {
                                    isoform_bounds = malloc(sizeof(Range));
                                    isoform_bounds->start = seq_idx;
                                }
                                if (length > 0) isoform_bounds->end = intron_start;

                                kputsn(hap_seq.s + seq_idx, intron_start - seq_idx, &spliced_seq);
                                seq_idx = intron_end;
                            }

                            if (isoform_bounds == NULL) {
                                isoform_bounds = malloc(sizeof(Range));
                                isoform_bounds->start = 0;
                            }
                            if (hap_seq.l != seq_idx) isoform_bounds->end = hap_seq.l;
                        } else if (t->strand == STRAND_REV) {
                            int length = 0;
                            for (int j = hap_ncrs.n-1; j >= 0; j--) {
                                NonCodingRegion ncr = kv_A(hap_ncrs, j);
                                int intron_start = ncr.coords.end - gene_start;
                                int intron_end = ncr.coords.start - gene_start;

                                if (intron_start < seq_idx) {
                                    continue;
                                }

                                int bases_shifted = range_len(ncr.canonical_coords) - range_len(ncr.coords);
                                length += intron_start - seq_idx;

                                if (length > 0 && isoform_bounds == NULL) {
                                    isoform_bounds = malloc(sizeof(Range));
                                    isoform_bounds->end = seq_idx;
                                }

                                if (length > 0) { 
                                    isoform_bounds->start = intron_start;
                                    kputsn(hap_seq.s + seq_idx + 1, intron_start - seq_idx, &spliced_seq);
                                }
                                seq_idx = intron_end - 1;
                            }

                            if (isoform_bounds == NULL) {
                                isoform_bounds = malloc(sizeof(Range));
                                isoform_bounds->end = 0;
                            }
                            if (hap_seq.l != seq_idx) isoform_bounds->start = hap_seq.l;
                        }

                        // Convert strand to complementary, if necessary (keeps direction the same)
                        if (t->strand == STRAND_REV) {
                            for (int si = 0; si < spliced_seq.l; si++) {
                                switch (spliced_seq.s[si]) {
                                    case 'A':
                                        spliced_seq.s[si] = 'T';
                                        break;
                                    case 'C':
                                        spliced_seq.s[si] = 'G';
                                        break;
                                    case 'G':
                                        spliced_seq.s[si] = 'C';
                                        break;
                                    case 'T':
                                        spliced_seq.s[si] = 'A';
                                        break;
                                }
                            }
                        }


                        Range shifted_isoform = { isoform_bounds->start + gene->coords.start, isoform_bounds->end + gene->coords.end };
                        if (t->strand == STRAND_REV) {
                            int64_t tmp = shifted_isoform.start;
                            shifted_isoform.start = shifted_isoform.end;
                            shifted_isoform.end = tmp;
                        }

                        if (t->strand == STRAND_FWD) {
                            if (hap_ncrs.a[0].is_utr) {
                                shifted_isoform.start = hap_ncrs.a[0].coords.end + 1;
                            }
                            if (hap_ncrs.a[hap_ncrs.n-1].is_utr) {
                                shifted_isoform.end = hap_ncrs.a[hap_ncrs.n-1].coords.start - 1;
                            }
                        } else if (t->strand == STRAND_REV) {
                            if (hap_ncrs.a[0].is_utr) {
                                shifted_isoform.end = hap_ncrs.a[0].coords.end + 1;
                            }
                            if (hap_ncrs.a[hap_ncrs.n-1].is_utr) {
                                shifted_isoform.start = hap_ncrs.a[hap_ncrs.n-1].coords.start - 1;
                            }
                        }

                        // Get isoform position for each record
                        kvec_t(IsoformVariant) variant_shifts; kv_init(variant_shifts);
                        for (int i = 0; i < hap_in.n; i++) {
                            Record *hap_r = kv_A(hap_in, i);
                            int64_t position = hap_r->v->pos;

                            int64_t isoform_position;
                            int allele_lengths[2];
                            const int coding = !get_range_in_isoform(hap_r->v, ncrs, &hap_ncrs, &hap_in, t->strand, gene->coords, &isoform_position, allele_lengths);
                            // isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, variant_to_range(hap_r->v, 0).start });


                            const int offset = coding ? (allele_lengths[0] - allele_lengths[1]) : 0; // Ref_len - alt_len
                            if (t->strand == STRAND_REV && offset != 0) { // Get position of the other end of the inserted sequence in the case of a reverse strand transcript
                                isoform_position += offset < 0 ? abs(offset) + 1 : 1;
                            }

                            IsoformVariant isoshift = {coding, 1, INDEL, isoform_position, offset, offset, hap_r};
                            kv_push(IsoformVariant, variant_shifts, isoshift);
                        }

                        // Get isoform position for splice site changes
                        kvec_t(IsoformVariant) splice_shifts; kv_init(splice_shifts);
                        int64_t **splice_positions = malloc(2 * ncrs->n * sizeof(int64_t *));
                        int64_t splice_offsets[2*ncrs->n];
                        int num_splice_sites = 0;
                        int strand_correction = t->strand == STRAND_FWD ? 1 : -1;
                        for (int i = 0; i < ncrs->n; i++) {
                            NonCodingRegion ncr = kv_A(*ncrs, i);
                            int64_t isoform_position;

                            for (int j = 0; j < ncr.records.n; j++) {
                                Record *r = kv_A(ncr.records, j);
                                enum UTRType utr_type = is_utr_variant(r, ncrs, &hap_in, t->strand, gene->coords);

                                if (utr_type == UTR5) set_record_flag(r, ti, phase, UTR5_VARIANT);
                                if (utr_type == UTR3) set_record_flag(r, ti, phase, UTR3_VARIANT);
                            }

                            if (ncr.donor != NULL && ncr.acceptor != NULL && ncr.donor->posterior_probability == 0 && ncr.acceptor->posterior_probability == 0) {
                                {
                                    int64_t position = ncr.canonical_coords.start;
                                    int offset = -(range_len(ncr.canonical_coords) - 1); // INTRON_RETENTION

                                    int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                    isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                    Records *records = (Records *) ncr.donor->records_ptr;
                                    for (int j = 0; j < records->n; j++) {
                                        Record *rr = kv_A(*records, j);

                                        IsoformVariant variant = {coding, 0, SPLICE, isoform_position, offset, offset, rr };
                                        kv_push(IsoformVariant, splice_shifts, variant);
                                        num_splice_sites++;
                                    }
                                }

                                {
                                    int64_t position = ncr.canonical_coords.end;
                                    int offset = -(range_len(ncr.canonical_coords) - 1); // INTRON_RETENTION

                                    int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                    isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                    Records *records = (Records *) ncr.acceptor->records_ptr;
                                    for (int j = 0; j < records->n; j++) {
                                        Record *rr = kv_A(*records, j);

                                        IsoformVariant variant = {coding, 0, SPLICE, isoform_position, offset, offset, rr };
                                        kv_push(IsoformVariant, splice_shifts, variant);
                                        num_splice_sites++;
                                    }
                                }
                                continue;
                            } else if (ncr.donor != NULL && ncr.donor->posterior_probability == 0) {
                                int64_t position = ncr.canonical_coords.start;
                                int offset = -(range_len(ncr.canonical_coords) - 1); // INTRON_RETENTION

                                int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                Records *records = (Records *) ncr.donor->records_ptr;
                                for (int j = 0; j < records->n; j++) {
                                    Record *rr = kv_A(*records, j);

                                    IsoformVariant variant = {coding, 0, SPLICE, isoform_position, offset, offset, rr };
                                    kv_push(IsoformVariant, splice_shifts, variant);
                                    num_splice_sites++;
                                }
                                continue;
                            } else if (ncr.acceptor != NULL && ncr.acceptor->posterior_probability == 0) {
                                int64_t position = ncr.canonical_coords.end;
                                int offset = -(range_len(ncr.canonical_coords) - 1); // INTRON_RETENTION

                                int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                Records *records = (Records *) ncr.acceptor->records_ptr;
                                for (int j = 0; j < records->n; j++) {
                                    Record *rr = kv_A(*records, j);

                                    IsoformVariant variant = {coding, 0, SPLICE, isoform_position, offset, offset, rr };
                                    kv_push(IsoformVariant, splice_shifts, variant);
                                    num_splice_sites++;
                                }
                                continue;
                            }

                            int donor_modified_bases = 0;
                            if (ncr.donor) {
                                donor_modified_bases = strand_correction * (ncr.donor->position - ncr.donor->canonical_position);
                            }
                            int donor_shift = strand_correction * (ncr.coords.start - ncr.canonical_coords.start);
                            if (ncr.donor != NULL) {
                                int64_t position = ncr.coords.start;
                                int64_t isoform_position;

                                int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                Records *records = (Records *) ncr.donor->records_ptr;
                                for (int j = 0; j < records->n; j++) {
                                    Record *rr = kv_A(*records, j);

                                    IsoformVariant variant = {coding, !ncr.uses_competitor, SPLICE, isoform_position, donor_shift, donor_modified_bases, rr };
                                    kv_push(IsoformVariant, splice_shifts, variant);
                                    num_splice_sites++;
                                }
                            }


                            int acceptor_modified_bases = 0;
                            if (ncr.acceptor) {
                                acceptor_modified_bases = strand_correction * (ncr.acceptor->position - ncr.acceptor->canonical_position);
                            }
                            int acceptor_shift = strand_correction * (ncr.coords.end - ncr.canonical_coords.end);
                            if (ncr.acceptor != NULL) {
                                int64_t position = ncr.coords.end;
                                int64_t isoform_position;

                                int coding = !get_position_in_isoform(position, ncrs, &hap_in, t->strand, gene->coords, &isoform_position);
                                isoform_position -= haplotype_offset(&hap_in, (Range) { gene->coords.start, position - 1 });

                                Records *records = (Records *) ncr.acceptor->records_ptr;
                                for (int j = 0; j < records->n; j++) {
                                    Record *rr = kv_A(*records, j);

                                    IsoformVariant variant = {coding, !ncr.uses_competitor, SPLICE, isoform_position, acceptor_shift, acceptor_modified_bases, rr };
                                    kv_push(IsoformVariant, splice_shifts, variant);
                                    num_splice_sites++;
                                }
                            }
                        }


                        kvec_t(IsoformVariant) shifts;
                        kv_init(shifts);

                        for (int i = 0; i < variant_shifts.n; i++) kv_push(IsoformVariant, shifts, kv_A(variant_shifts, i));
                        for (int i = 0; i < splice_shifts.n; i++) kv_push(IsoformVariant, shifts, kv_A(splice_shifts, i));
                        if (shifts.n > 0 && t->strand == STRAND_FWD) ks_introsort_frameshift_sort(shifts.n, shifts.a);
                        if (shifts.n > 0 && t->strand == STRAND_REV) ks_introsort_frameshift_reverse_sort(shifts.n, shifts.a);

                        // for (int i = 0; i < shifts.n; i++) {
                        //     IsoformVariant iso_shift = kv_A(shifts, i);
                        //     if (iso_shift.r != r) continue;
                        //
                        //     for (int j = 0; j < i; j++) {
                        //         IsoformVariant upstream_shift = kv_A(shifts, j);
                        //         if (upstream_shift.r == iso_shift.r) continue;
                        //
                        //         for (int ti2 = 0; ti2 < upstream_shift.r->tscripts->n; ti2++) {
                        //             Transcript *tscript = kv_A(*upstream_shift.r->tscripts, ti2);
                        //             if (tscript == NULL || tscript->id != t->id) continue;
                        //
                        //             TscriptCsqs *csq = kv_A(*upstream_shift.r->vep, ti2);
                        //             if (!upstream_shift.r->isoforms->valid | 
                        //                 (!(
                        //                  (upstream_shift.r->isoforms->rescue == 0) & 
                        //                  ((csq->triggers[SPLICE_TRIGGER] == CANONICAL_SPLICE_VALUE) |
                        //                  (csq->triggers[FRAMESHIFT_TRIGGER]) |
                        //                  (csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE))
                        //                  ) & 
                        //                 !upstream_shift.r->isoforms->ptv)) continue;
                        //
                        //             set_record_flag(r, ti, phase, UPSTREAM_PTV);
                        //         }
                        //     }
                        // }

                        // Filter to shifts that overlap the isoform
                        size_t j = 0; // Index for the next valid position
                        for (size_t i = 0; i < kv_size(shifts); i++) {
                            IsoformVariant shift = kv_A(shifts, i);
                            if (shift.r == NULL) continue;

                            int vposition = shift.r->v->pos - haplotype_offset(&hap_in, (Range) { gene->coords.start, shift.r->v->pos });
                            if ((shift.shift_type == SPLICE && shift.shift != 0) || (shift.shift_type == INDEL && shift.coding && vposition >= shifted_isoform.start && vposition < shifted_isoform.end)) {
                                kv_A(shifts, j++) = kv_A(shifts, i);
                            }
                        }
                        shifts.n = j; // Adjust the size of the vector


                        // Reverse direction if necessary, 3' - 5' --> 5' - 3'
                        log_info("%s:%li (isoform_length=%i)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos, spliced_seq.l);
                        if (t->strand == STRAND_REV) {
                            for (int first = 0, last = spliced_seq.l-1; first < last; first++, last--) {
                                const char tmp = spliced_seq.s[first];
                                spliced_seq.s[first] = spliced_seq.s[last];
                                spliced_seq.s[last] = tmp;
                            }
                        }

                        // Forward
                        int shift_index = 0;

                        const char stop_codons[3][3] = STOP_CODONS;

                        int affected_bases = 0;
                        int rescued = 0;
                        int stop_checked = 0, frameshift_checked = 0;
                        int stop_rescued = 0;
                        enum RescueType rescue_type = NO_RESCUE;
                        enum PTVType ptv_type = NO_PTV;
                        int touched = 0;
                        for (int ri = 0; ri < shifts.n; ri++) {
                            IsoformVariant iso_shift = kv_A(shifts, ri);

                            if (iso_shift.r != r) continue;
                            if (stop_checked) continue;

                            // Check if constutive and in frame
                            for (int ei = 0; ei < t->nexons; ei++) {
                                Exon exon = *t->exons[ei];
                                if (exon.cds == NULL) continue;
                                if ((iso_shift.r->v->pos < exon.cds->start) | (iso_shift.r->v->pos > exon.cds->end)) continue;

                                if (!exon.constitutive) set_record_flag(r, ti, phase, NON_CONSTITUTIVE);
                                if (exon.in_frame) set_record_flag(r, ti, phase, INFRAME);
                                break;
                            }

                            touched = 1;
                            if (csq->triggers[STOP_TRIGGER]) {
                                int codon_start;
                                if (t->strand == STRAND_FWD) {
                                    codon_start = (int) ((iso_shift.seq_pos - shifted_isoform.start + 1) / 3) * 3;
                                } else {
                                    codon_start = (int) ((shifted_isoform.end - iso_shift.seq_pos - 1) / 3) * 3;
                                }

                                const char *codon = spliced_seq.s + codon_start;

                                Range variant = variant_to_range(iso_shift.r->v, 1);

                                // Check if the current codon codes for a stop
                                int is_stop = 0;
                                for (int i = 0; i < range_len(variant); i+=3) {
                                    const char *variant_codon = codon + i;
                                    for (int j = 0; j < 3; j++) {
                                        const char *stop = stop_codons[j];
                                        if (!(variant_codon[0] == stop[0] && variant_codon[1] == stop[1] && variant_codon[2] == stop[2])) continue;

                                        is_stop = 1;
                                        break;
                                    }

                                    if (is_stop) break;
                                }

                                if (!is_stop && csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE) {
                                    // log_error("%d:%d t=%i if=%i STOP RESCUE", r->v->rid, r->v->pos, ti, i);
                                    if (rescue_type == NO_RESCUE) {
                                        rescue_type = COMPLETE_RESCUE;
                                        rescued = 1;
                                        stop_rescued = 1;
                                    }
                                } else if (csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE) {
                                    stop_checked = 1;
                                }

                                if (is_stop && (csq->triggers[STOP_TRIGGER] == STOP_GAIN_VALUE)) {
                                    ptv_type = STOP_GAINED;
                                }
                            }

                            int should_check_induced = 0;
                            Range possible_induced_sequence = {0};

                            if (csq->triggers[SPLICE_TRIGGER] && !stop_checked && iso_shift.shift_type == SPLICE) {
                                if (iso_shift.shift == 0) {
                                    rescue_type = COMPLETE_RESCUE;
                                    rescued = 1;
                                    continue;
                                }

                                if (iso_shift.partial) {
                                    rescue_type = iso_shift.shift < 0 ? PARTIAL_INTRON_RETENTION : PARTIAL_EXON_SKIPPING;
                                } else {
                                    rescue_type = iso_shift.shift < 0 ? INTRON_RETENTION : EXON_SKIPPING;
                                }
                                affected_bases = abs(iso_shift.shift);

                                if (iso_shift.shift % 3 == 0 && abs(iso_shift.shift) < MAX_SPLICE_RADIUS) {
                                    rescued = 1;
                                    int splice_insertion_length = iso_shift.modified_bases < 0 ? abs(iso_shift.modified_bases): 0;
                                    should_check_induced = 1;

                                    if (t->strand == STRAND_FWD) {
                                        possible_induced_sequence.start = (int) ((iso_shift.seq_pos - shifted_isoform.start + 1) / 3) * 3;
                                        possible_induced_sequence.end = (iso_shift.seq_pos - shifted_isoform.start + 1) + splice_insertion_length + 1;
                                    } else {
                                        possible_induced_sequence.start = (int) ((shifted_isoform.end - iso_shift.seq_pos - 1) / 3) * 3;
                                        possible_induced_sequence.end = (shifted_isoform.end - iso_shift.seq_pos - 1) + splice_insertion_length + 1;
                                    }
                                } else if (csq->triggers[SPLICE_TRIGGER] && iso_shift.shift % 3 != 0) {
                                    rescued = 0;
                                    ptv_type = FRAMESHIFT;
                                } else if (csq->triggers[SPLICE_TRIGGER] && abs(iso_shift.shift) > MAX_SPLICE_RADIUS) {
                                    ptv_type = LONG_INFRAME_INDEL;

                                    int splice_insertion_length = iso_shift.modified_bases < 0 ? abs(iso_shift.modified_bases) : 0;
                                    should_check_induced = 1;

                                    if (t->strand == STRAND_FWD) {
                                        possible_induced_sequence.start = (int) ((iso_shift.seq_pos - shifted_isoform.start + 1) / 3) * 3;
                                        possible_induced_sequence.end = (iso_shift.seq_pos - shifted_isoform.start + 1) + splice_insertion_length + 1;
                                    } else {
                                        possible_induced_sequence.start = (int) ((shifted_isoform.end - iso_shift.seq_pos - 1) / 3) * 3;
                                        possible_induced_sequence.end = (shifted_isoform.end - iso_shift.seq_pos - 1) + splice_insertion_length + 1;
                                    }
                                }
                            }

                            if ((csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[SPLICE_TRIGGER]) && (ptv_type == FRAMESHIFT || ptv_type == NO_RESCUE || rescued == 1) && !stop_checked) {
                                if (iso_shift.shift == 0) continue;
                                frameshift_checked = 1;


                                if (iso_shift.shift % 3 == 0 && abs(iso_shift.shift) < MAX_FRAMESHIFT_RADIUS) {
                                    rescued = 1;
                                    if (iso_shift.partial) {
                                        rescue_type = iso_shift.shift < 0 ? PARTIAL_INTRON_RETENTION : PARTIAL_EXON_SKIPPING;
                                    } else {
                                        rescue_type = iso_shift.shift < 0 ? INTRON_RETENTION : EXON_SKIPPING;
                                    }
                                    // rescue_type = FRAME_RESTORING_INDEL;
                                    affected_bases = abs(iso_shift.shift);
                                } else {
                                    int32_t summed_shift = 0;
                                    int32_t shift_start = 0;
                                    int32_t starts_with_forward_indel = 0;
                                    int included_current_shift = 0;
                                    int induced_frameshift_rescue = 0;
                                    // int inserted_bases = 0;
                                    // int deleted_bases = 0;
                                    for (int i = 0; i < shifts.n; i++) {
                                        IsoformVariant shift = kv_A(shifts, i);
                                        if (shift.shift == 0) continue; // Rescued splice, do nothing
                                        if (abs(shift.shift) > MAX_FRAMESHIFT_RADIUS) continue;

                                        // If outside of rescue radius
                                        if (abs(((int32_t) shift_start) - ((int32_t) shift.seq_pos)) > MAX_FRAMESHIFT_RADIUS) {
                                            summed_shift = shift.shift;
                                            shift_start = shift.seq_pos;
                                            included_current_shift = 0 | (shift.r == r);
                                            // inserted_bases += shift.shift < 0 ? abs(shift.shift) : 0;
                                            // deleted_bases += shift.shift > 0 ? abs(shift.shift) : 0;
                                            induced_frameshift_rescue = shift.shift_type == SPLICE;
                                            starts_with_forward_indel = shift.shift_type == INDEL;
                                            continue;
                                        }

                                        // Keep track of whether the set of frameshifts includes the target 
                                        if (!included_current_shift && shift.r == r) included_current_shift = 1;

                                        summed_shift += shift.shift;
                                        // inserted_bases += shift.shift < 0 ? abs(shift.shift) : 0;
                                        // deleted_bases += shift.shift > 0 ? abs(shift.shift) : 0;

                                        induced_frameshift_rescue = shift.shift_type == SPLICE && shift.shift % 3 != 0;
                                        if (summed_shift % 3 != 0) continue;

                                        if (included_current_shift && shift.seq_pos + MAX_FRAMESHIFT_RADIUS > iso_shift.seq_pos && shift.seq_pos - MAX_FRAMESHIFT_RADIUS < iso_shift.seq_pos) {
                                            rescued = 1;
                                            rescue_type = (rescue_type == NO_RESCUE) ? FRAME_RESTORING_INDEL : rescue_type; // Only set to rescue type
                                            int32_t shift_end = shift.seq_pos;
                                            affected_bases = abs(shift_end - shift_start) + (shift.modified_bases < 0 ? abs(shift.modified_bases) : 0);

                                            if (induced_frameshift_rescue) r->flags[ti * NUM_PHASES + phase] |= RESCUED_BY_SPLICE_FRAMESHIFT;

                                            if (ptv_type == FRAMESHIFT) ptv_type = NO_PTV;

                                            should_check_induced = 1;
                                            if (t->strand == STRAND_FWD) {
                                                possible_induced_sequence.start = (int) ((shift_start - shifted_isoform.start + 1) / 3) * 3;
                                                possible_induced_sequence.end = (shift_start - shifted_isoform.start + 1) + affected_bases;
                                            } else {
                                                possible_induced_sequence.start = (int) ((shifted_isoform.end - shift_start - 1) / 3) * 3;
                                                possible_induced_sequence.end = (shifted_isoform.end - shift_start - 1) + affected_bases;
                                            }

                                        }

                                        summed_shift = 0;
                                        shift_start = 0;
                                        included_current_shift = 0;
                                        starts_with_forward_indel = 0;
                                        // inserted_bases = 0;
                                        // deleted_bases = 0;
                                    }
                                }
                            }

                            if (!should_check_induced) continue;

                            const char *codon = spliced_seq.s + possible_induced_sequence.start;

                            // Check if the current codon codes for a stop
                            int is_stop = 0;
                            for (int i = 0; i < range_len(possible_induced_sequence); i+=3) {
                                const char *variant_codon = codon + i;
                                for (int j = 0; j < 3; j++) {
                                    const char *stop = stop_codons[j];
                                    if (!(variant_codon[0] == stop[0] && variant_codon[1] == stop[1] && variant_codon[2] == stop[2])) continue;

                                    is_stop = 1;
                                    break;
                                }
                            }

                            if (is_stop) {
                                // log_error("%d:%d t=%i if=%i STOP GAIN", r->v->rid, r->v->pos, ti, i);
                                ptv_type = STOP_GAINED;
                                continue;
                            }

                        }

                        if (!touched && !stop_checked && !frameshift_checked && rescue_type == NO_RESCUE) {
                            rescued = 1;
                            rescue_type = COMPLETE_RESCUE;
                        }

                        for (int i = 0; i < shifts.n; i++) {
                            IsoformVariant iso_shift = kv_A(shifts, i);
                            if (iso_shift.r != r) continue;

                            if (
                                (csq->triggers[SPLICE_TRIGGER] == AF_SPLICE_VALUE || csq->triggers[STOP_TRIGGER] == STOP_GAIN_VALUE) &&
                                (csq->triggers[SPLICE_TRIGGER] != CANONICAL_SPLICE_VALUE && !csq->triggers[FRAMESHIFT_TRIGGER] && csq->triggers[STOP_TRIGGER] != STOP_RESCUE_VALUE) &&
                                (ptv_type == NO_PTV || ptv_type == LONG_INFRAME_INDEL)
                            ) continue;
                            if (rescued && !(ptv_type == STOP_GAINED || ptv_type == FRAMESHIFT)) continue;

                            for (int j = i; j < shifts.n; j++) {
                                IsoformVariant upstream_shift = kv_A(shifts, j);
                                if (upstream_shift.r == iso_shift.r) continue;

                                for (int ti2 = 0; ti2 < upstream_shift.r->tscripts->n; ti2++) {
                                    Transcript *tscript = kv_A(*upstream_shift.r->tscripts, ti2);
                                    if (tscript == NULL || tscript->id != t->id) continue;

                                    TscriptCsqs *csq = kv_A(*upstream_shift.r->vep, ti2);
                                    if (!(upstream_shift.r->isoforms->valid && (csq->triggers[SPLICE_TRIGGER] || csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[STOP_TRIGGER]))) continue;

                                    set_record_flag(upstream_shift.r, ti2, phase, UPSTREAM_PTV);
                                }
                            }
                        }

                        float isoform_prior = 1.0;
                        for (int i = 0; i < hap_ncrs.n; i++) {
                            NonCodingRegion ncr = kv_A(hap_ncrs, i);
                            if (ncr.coords.start == ncr.coords.end) continue; // Skips should not be included

                            if (ncr.donor != NULL) {
                                float prob = ncr.donor->prior_probability;
                                if (ncr.donor->posterior_probability == 0) prob = 1 - prob;
                                isoform_prior *= prob;
                            }
                            if (ncr.acceptor != NULL) {
                                float prob = ncr.acceptor->prior_probability;
                                if (ncr.acceptor->posterior_probability == 0) prob = 1 - prob;
                                isoform_prior *= prob;
                            }
                        }

                        Isoform *isoform = malloc(sizeof(Isoform));
                        isoform->prior_probability = isoform_prior;
                        isoform->posterior_probability = 0.0;
                        isoform->rescue = rescued;
                        isoform->rescue_type = rescue_type;
                        isoform->affected_bases = affected_bases;

                        isoform->ptv = (ptv_type == STOP_GAINED || ptv_type == FRAMESHIFT || ptv_type == LONG_INFRAME_INDEL);
                        isoform->ptv_type = ptv_type;

                        isoforms.isoforms[isoforms.n_isoforms++] = isoform;
                    }

                    // Calculate isoforms posterior
                    if (isoforms.n_isoforms != isoforms.m_isoforms) {
                        log_error("Less isoforms found than expected. This is a bug!");
                        exit(1);
                    }

                    float probability_sum = 0.0;
                    for (int i = 0; i < isoforms.n_isoforms; i++) probability_sum += isoforms.isoforms[i]->prior_probability;
                    for (int i = 0; i < isoforms.n_isoforms; i++) isoforms.isoforms[i]->posterior_probability = isoforms.isoforms[i]->prior_probability / probability_sum;

                    calculate_rescue_probabilities(&isoforms);
                    calculate_ptv_probabilities(&isoforms);

                    r->isoforms[ti * NUM_PHASES + phase] = isoforms;

                    kv_destroy(hap_in);
                    free(hap_seq.s);
                }

            }



            // Clean up records
            kv_destroy(records);
        }


        /**
            * FORMATTING AND ADDING TAGS
        **/

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);

            int splice_triggered = 0;
            for (int j = 0; j < r->vep->n; j++) {
                if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) splice_triggered = 1;
            }

            if (!splice_triggered) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            kstring_t splice_score = {0};

            int has_scores = 0;
            for (int g = 0; g < r->genes->n; g++) {
                if (g > 0) {
                    kputc(',', &splice_score);
                }

                for (int p = 0; p < NUM_PHASES; p++) {

                    if (p > 0) {
                        kputc('&', &splice_score);
                    }

                    Score *score = r->native_scores[g * NUM_PHASES + p];
                    kputs(r->v->d.allele[1], &splice_score);
                    if (score == NULL) {
                        kputc('|', &splice_score); kputs(kv_A(*r->genes, g)->name, &splice_score);
                        kputs("|.|.|.|.|.|.|.|.", &splice_score);
                        continue;
                    }
                    has_scores = 1;

                    kputc('|', &splice_score); kputs(score->gene, &splice_score);

                    // Format and append DS scores
                    char ds_scores[20];
                    snprintf(ds_scores, 20, "%.2f|%.2f|%.2f|%.2f", score->ag, score->al, score->dg, score->dl);
                    kputc('|', &splice_score); kputsn(ds_scores, 19, &splice_score);

                    kputc('|', &splice_score); kputw(score->ag_idx, &splice_score);
                    kputc('|', &splice_score); kputw(score->al_idx, &splice_score);
                    kputc('|', &splice_score); kputw(score->dg_idx, &splice_score);
                    kputc('|', &splice_score); kputw(score->dl_idx, &splice_score);

                }
            }

            if (has_scores) log_info(splice_score.s);

            bcf_update_info_string(hdr, r->v, SPLICEAI_TAG, splice_score.s);
        }

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);

            int splice_triggered = 0;
            for (int j = 0; j < r->vep->n; j++) {
                if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) splice_triggered = 1;
            }

            if (!splice_triggered) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            kstring_t diff_score = {0};

            int has_scores = 0;
            for (int g = 0; g < r->genes->n; g++) {
                if (g > 0) {
                    kputc(',', &diff_score);
                }

                for (int p = 0; p < NUM_PHASES; p++) {

                    if (p > 0) {
                        kputc('&', &diff_score);
                    }

                    Score *score = r->diff_scores[g * NUM_PHASES + p];
                    kputs(r->v->d.allele[1], &diff_score);
                    if (score == NULL) {
                        kputc('|', &diff_score); 
                        kputs(kv_A(*r->genes, g)->name, &diff_score);
                        kputs("|.|.|.|.|.|.|.|.", &diff_score);
                        continue;
                    }
                    has_scores = 1;

                    kputc('|', &diff_score); kputs(score->gene, &diff_score);

                    // Format and append DS scores
                    char ds_scores[24];
                    snprintf(ds_scores, 23, "%.2f|%.2f|%.2f|%.2f", score->ag, score->al, score->dg, score->dl);

                    kputc('|', &diff_score); kputs(ds_scores, &diff_score);

                    kputc('|', &diff_score); kputw(score->ag_idx, &diff_score);
                    kputc('|', &diff_score); kputw(score->al_idx, &diff_score);
                    kputc('|', &diff_score); kputw(score->dg_idx, &diff_score);
                    kputc('|', &diff_score); kputw(score->dl_idx, &diff_score);

                }
            }

            if (has_scores) log_info(diff_score.s);

            bcf_update_info_string(hdr, r->v, CS_TAG, diff_score.s);
        }

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);

            int splice_triggered = 0;
            for (int j = 0; j < r->vep->n; j++) {
                if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) splice_triggered = 1;
            }

            if (!splice_triggered) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            kstring_t haplotype_splice_scores = {0};

            int has_scores = 0;
            for (int g = 0; g < r->genes->n; g++) {
                if (g > 0) {
                    kputc(',', &haplotype_splice_scores);
                }

                for (int p = 0; p < NUM_PHASES; p++) {

                    if (p > 0) {
                        kputc('&', &haplotype_splice_scores);
                    }

                    Score *score = r->haplotype_scores[g * NUM_PHASES + p];
                    kputs(r->v->d.allele[1], &haplotype_splice_scores);
                    if (score == NULL) {
                        kputc('|', &haplotype_splice_scores); kputs(kv_A(*r->genes, g)->name, &haplotype_splice_scores);
                        kputs("|.|.|.|.|.|.|.|.", &haplotype_splice_scores);
                        continue;
                    }
                    has_scores = 1;

                    kputc('|', &haplotype_splice_scores); kputs(score->gene, &haplotype_splice_scores);

                    // Format and append DS scores
                    char ds_scores[20];
                    snprintf(ds_scores, 20, "%.2f|%.2f|%.2f|%.2f", score->ag, score->al, score->dg, score->dl);
                    kputc('|', &haplotype_splice_scores); kputsn(ds_scores, 19, &haplotype_splice_scores);

                    kputc('|', &haplotype_splice_scores); kputw(score->ag_idx, &haplotype_splice_scores);
                    kputc('|', &haplotype_splice_scores); kputw(score->al_idx, &haplotype_splice_scores);
                    kputc('|', &haplotype_splice_scores); kputw(score->dg_idx, &haplotype_splice_scores);
                    kputc('|', &haplotype_splice_scores); kputw(score->dl_idx, &haplotype_splice_scores);

                }
            }

            if (has_scores) log_info(haplotype_splice_scores.s);

            bcf_update_info_string(hdr, r->v, SPLICEAI_HAP_TAG, haplotype_splice_scores.s);

            free(haplotype_splice_scores.s);
        }

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);

            int triggered = 0;
            int has_ptv_effect = 0;
            for (int j = 0; j < r->vep->n; j++) {
                TscriptCsqs *csq = kv_A(*r->vep, j);
                if (csq->triggers[SPLICE_TRIGGER] == CANONICAL_SPLICE_VALUE || csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE) triggered = 1;
                if (r->isoforms->valid && r->isoforms->ptv_probability > 0) has_ptv_effect = 1;
            }

            if (!triggered && !has_ptv_effect) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            kstring_t rescue_tag = {0};
            kstring_t rescue_prob_tag = {0};
            kstring_t rescue_type_tag = {0};

            for (int ti = 0; ti < r->transcript_genes->n; ti++) {
                if (ti > 0) {
                    kputc(',', &rescue_tag);
                    kputc(',', &rescue_prob_tag);
                    kputc(',', &rescue_type_tag);
                }

                if (kv_A(*r->transcript_genes, ti) == NULL) {
                    kputc('.', &rescue_tag);
                    kputc('.', &rescue_prob_tag);
                    kputc('.', &rescue_type_tag);
                    continue;
                }

                TscriptCsqs *csq = kv_A(*r->vep, ti);

                // No Rescue triggers at all
                if ((csq->triggers[SPLICE_TRIGGER] != CANONICAL_SPLICE_VALUE) & (!csq->triggers[FRAMESHIFT_TRIGGER]) & (csq->triggers[STOP_TRIGGER] != STOP_RESCUE_VALUE)) {
                    kputc('.', &rescue_tag);
                    kputc('.', &rescue_prob_tag);
                    kputc('.', &rescue_type_tag);
                    continue;
                }

                for (int p = 0; p < NUM_PHASES; p++) {
                    if (p > 0) {
                        kputc('&', &rescue_tag);
                        kputc('&', &rescue_prob_tag);
                        kputc('&', &rescue_type_tag);
                    }

                    Isoforms isoforms = r->isoforms[ti * NUM_PHASES + p];
                    if (!isoforms.valid | (
                        (((csq->triggers[STOP_TRIGGER] == STOP_GAIN_VALUE) | (csq->triggers[SPLICE_TRIGGER] == AF_SPLICE_VALUE)) & (isoforms.rescue == COMPLETE_RESCUE)) &
                        !((csq->triggers[SPLICE_TRIGGER] == CANONICAL_SPLICE_VALUE) | (csq->triggers[FRAMESHIFT_TRIGGER]) | (csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE)))
                    ) {
                        kputc('.', &rescue_tag);
                        kputc('.', &rescue_prob_tag);
                        kputs(".|.|.|.", &rescue_type_tag);
                        continue;
                    }

                    kputw(isoforms.rescue, &rescue_tag);

                    char prob_s[7];
                    sprintf(prob_s, "%.4f", isoforms.rescue_probability);
                    kputsn(prob_s, 6, &rescue_prob_tag);

                    kstring_t transcript_phase_tag = generate_rescue_type_tag(&isoforms);
                    kputsn(transcript_phase_tag.s, transcript_phase_tag.l, &rescue_type_tag);
                }
            }

            bcf_update_info_string(hdr, r->v, RESCUE_TAG, rescue_tag.s);
            bcf_update_info_string(hdr, r->v, RESCUE_PROB_TAG, rescue_prob_tag.s);
            bcf_update_info_string(hdr, r->v, RESCUE_TYPE_TAG, rescue_type_tag.s);
        }

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);

            int triggered = 0;
            for (int j = 0; j < r->vep->n; j++) {
                TscriptCsqs *csq = kv_A(*r->vep, j);
                if (csq->triggers[SPLICE_TRIGGER] || csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[STOP_TRIGGER]) triggered = 1;
            }

            if (!triggered) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            int has_ptv = 0;
            kstring_t ptv_tag = {0};
            kstring_t ptv_prob_tag = {0};
            kstring_t ptv_type_tag = {0};

            for (int ti = 0; ti < r->transcript_genes->n; ti++) {
                if (ti > 0) {
                    kputc(',', &ptv_tag);
                    kputc(',', &ptv_prob_tag);
                    kputc(',', &ptv_type_tag);
                }


                for (int p = 0; p < NUM_PHASES; p++) {
                    if (p > 0) {
                        kputc('&', &ptv_tag);
                        kputc('&', &ptv_prob_tag);
                        kputc('&', &ptv_type_tag);
                    }

                    Isoforms isoforms = r->isoforms[ti * NUM_PHASES + p];

                    if (!isoforms.valid) {
                        kputc('.', &ptv_tag);
                        kputc('.', &ptv_prob_tag);
                        kputc('.', &ptv_type_tag);

                        continue;
                    }

                    kputw(isoforms.ptv, &ptv_tag);

                    char prob_s[7];
                    sprintf(prob_s, "%.4f", isoforms.ptv_probability);
                    kputsn(prob_s, 6, &ptv_prob_tag);

                    switch(isoforms.ptv_type) {
                        case STOP_GAINED:
                            kputs("stop_gained", &ptv_type_tag);
                            has_ptv = 1;
                            break;
                        case FRAMESHIFT:
                            kputs("frameshift", &ptv_type_tag);
                            has_ptv = 1;
                            break;
                        case LONG_INFRAME_INDEL:
                            kputs("long_inframe_indel", &ptv_type_tag);
                            has_ptv = 1;
                            break;
                        case NO_PTV:
                            kputs("no_ptv", &ptv_type_tag);
                            break;
                    }
                }

            }

            // if (!has_ptv) continue;

            bcf_update_info_string(hdr, r->v, PTV_TAG, ptv_tag.s);
            bcf_update_info_string(hdr, r->v, PTV_PROB_TAG, ptv_prob_tag.s);
            bcf_update_info_string(hdr, r->v, PTV_TYPE_TAG, ptv_type_tag.s);

        }

        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);
            if (!block.is_coding || r->transcript_genes == NULL) continue;

            int triggered = 0;
            for (int j = 0; j < r->vep->n; j++) {
                TscriptCsqs *csq = kv_A(*r->vep, j);
                if (csq->triggers[SPLICE_TRIGGER] || csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[STOP_TRIGGER]) triggered = 1;
            }

            if (!triggered) continue;
            if (strlen(r->v->d.allele[1]) > CAI_DISTANCE || r->v->rlen > CAI_DISTANCE) continue;

            kstring_t tags = {0};

            for (int t = 0; t < r->transcript_genes->n; t++) {
                if (t > 0) { kputc(',', &tags); }

                for (int p = 0; p < NUM_PHASES; p++) {
                    if (p > 0) { kputc('&', &tags); }

                    uint32_t flags = r->flags[t * NUM_PHASES + p];
                    if (!flags) { 
                        kputc('.', &tags); 
                        continue;
                    }

                    kstring_t s = {0};
                    get_prediction_flags(flags, &s);
                    kputs(s.s, &tags);
                }
            }

            bcf_update_info_string(hdr, r->v, PREDICTION_FLAGS_TAG, tags.s);
        }

        // Write out records
        for (int i = 0; i < block.records.n; i++) {
            Record *r = kv_A(block.records, i);
            bcf_write(fp_out, hdr, r->v);
        }

        destroy_record_block(block);
    }

    // Clean up remaining structures
    destroy_hdr_fmt(vep_fmt);
    destroy_models(models);
    bcf_destroy(buffer.v);

    // Clean up remaining file handles
    bcf_close(fp_in);
    bcf_close(fp_out);
    bcf_hdr_destroy(hdr);
    fai_destroy(fai);
    gff_destroy(gff);

    // Free raw scores file
    fclose(scores_tsv);
}

