#include <cpliceai/range.h>
#include <cpliceai/utils.h>
#include <htslib/kstring.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <klib/kvec.h>
#include <klib/kstring.h>

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
#define CRYPTIC_SPLICE_VALUE 2

#define PROXIMITY_THRESHOLD 15

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

#define RESCUE_TYPE_TAG "RESCUE_TYPE"
#define RESCUE_TYPE_DESC "##INFO=<ID=RESCUE_TYPE,Number=.,Type=String,Description=\"Rescue type\">"

#define PTV_TAG "PTV"
#define PTV_DESC "##INFO=<ID=PTV,Number=.,Type=String,Description=\"Newly predicted PTV consequence\">"

#define PTV_TYPE_TAG "PTV_TYPE"
#define PTV_TYPE_DESC "##INFO=<ID=PTV_TYPE,Number=.,Type=String,Description=\"Type of PTV\">"

#define PREDICTION_FLAGS_TAG "PREDICTION_FLAGS"
#define PREDICTION_FLAGS_DESC "##INFO=<ID=PREDICTION_FLAGS,Number=.,Type=String,Description=\"Additional flags for the context in which the prediction happened.\">"

#define INFO_TAGS { SPLICEAI_HAP_TAG, CS_TAG, RESCUE_TAG, RESCUE_TYPE_TAG, PTV_TAG, PTV_TYPE_TAG, PREDICTION_FLAGS_TAG }

enum SpliceType {
    ACCEPTOR = 1,
    DONOR = 2
};

typedef struct SpliceProbability {
    enum SpliceType type;
    float prior_probability;
    float posterior_probability;
    uint32_t position;
} SpliceProbability;

typedef struct RecordBlock {
    int is_coding;
    Records *records;
} RecordBlock;

typedef kvec_t(Records *) PhasedRecords;


Range variant_to_range(const bcf1_t *v, const int allele) {
    return (Range) { v->pos, v->pos + strlen(v->d.allele[allele]) };
}

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

            if (strcmp(tid, kv_A(*r->vep, i).tid) != 0) continue;

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

    if (predict(models, enc_len, 1, &ohe, num_predictions, predictions)) {
        log_info("Predict failed.");
        return;
    }
    free(ohe);

    if (strand == STRAND_REV) reverse_prediction(*predictions, *num_predictions, 3);
}

void get_best_insertion_score(const float *insert_start, const int length, float *best) {
    for (int i = 0; i < length; i+=3) {
        for (int j = 0; j < 3; j++) {
            if (insert_start[i + j] > best[j]) best[j] = insert_start[i+j];
        }
    }
}

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

void remove_unworkable_variants(Records *vs, const bcf_hdr_t *hdr) {
    Range last_variant = {0};
    IRegion last_deletion = {-1, -1, -1};
    int *should_delete = calloc(vs->n, sizeof(int));
    for (int j = 0; j < vs->n; j++) {
        bcf1_t *v = kv_A(*vs, j)->v;

        Range bounds = variant_to_range(v, 0);
        Range intersection;
        int overlaps = intersect_ranges(bounds, last_variant, &intersection);

        if (overlaps && range_len(intersection) > 0) {
            log_warn("Skipping on %s:%li, overlaps previous variant", bcf_hdr_id2name(hdr, v->rid), v->pos);
            should_delete[j] = 1;
        }

        // if (strlen(v->d.allele[1]) > CAI_DISTANCE) {
        //     log_warn("Skipping on %s:%li, insertion is too long", bcf_hdr_id2name(hdr, v->rid), v->pos);
        //     should_delete[j] = 1;
        // }
        // if (v->rlen > CAI_DISTANCE) {
        //     log_warn("Skipping on %s:%li, deletion is too long", bcf_hdr_id2name(hdr, v->rid), v->pos);
        //     should_delete[j] = 1;
        // }

        if (range_len(intersection) > 0) {
            last_variant = (Range) { last_variant.start, last_variant.end > bounds.end ? last_variant.end : bounds.end };
        } else {
            last_variant = bounds;
        }
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

int prediction_offset(bcf1_t *v, Records *rs) {
    Range pred_range = {v->pos - (CAI_DISTANCE), v->pos + CAI_DISTANCE+1};

    int offset = 0;
    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        Range ref;
        if (intersect_ranges((Range) {r->v->pos, r->v->pos + r->v->rlen}, pred_range, &ref) == 0) continue;
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt;
        int alt_len = 0;
        if (intersect_ranges((Range) {r->v->pos, r->v->pos + strlen(r->v->d.allele[1])}, pred_range, &alt) > 0) {
            alt_len = range_len(alt);
        }
        offset += ref_len - alt_len;
    }

    return offset;
}

int find_gene_region_for_variant(bcf1_t *v, bcf_hdr_t *hdr, gff_t *gff, IRegion *region) {
    regidx_t *genes = gff_get(gff, idx_gene);
    regitr_t *itr = regitr_init(genes);
    const char *chr = bcf_hdr_id2name(hdr, v->rid);

    if (!regidx_overlap(genes, chr, v->pos, v->pos + v->rlen, itr)) { regitr_destroy(itr); return -1; };

    IRegion r = {0};
    r.rid = v->rid;
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
    if (r.start == 0 && r.end == 0) return -1; // No protein coding gene was found

    *region = r;

    return 0;
}

int allele_offset(const bcf1_t *v, const int ref_index, const int alt_index, const Range bounds) {
    Range ref = variant_to_range(v, ref_index);
    Range ref_range;
    if (intersect_ranges(ref, bounds, &ref_range) == 0) {
        return 0;
    }

    Range alt = variant_to_range(v, alt_index);
    Range alt_range;
    int alt_len = 0;
    if (intersect_ranges(alt, bounds, &alt_range) != 0) {
        alt_len = range_len(alt_range);
    }

    return range_len(ref_range) - alt_len;
}

int haplotype_offset(const Records *rs, const Range bounds) {
    int offset = 0;

    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);
        offset += allele_offset(r->v, 0, 1, bounds);
    }

    return offset;
}

float *match_predictions_to_reference(const float *predictions, const bcf1_t *v, const Records *rs, const size_t num_reference_predictions) {
    float *corrected_predictions = malloc(num_reference_predictions * sizeof(float));
    Range pred_range = {v->pos - (CAI_DISTANCE), v->pos + CAI_DISTANCE+1}; // +1 for open ended coordinates

    int ref_idx = 0, hap_idx = 0;
    for (int l = 0; l < rs->n; l++) {
        Record *r = kv_A(*rs, l);

        Range ref;
        if (intersect_ranges(variant_to_range(r->v, 0), pred_range, &ref) == 0) continue;
        shift_range(&ref, -pred_range.start);
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt;
        int alt_len = 0;
        if (intersect_ranges(variant_to_range(r->v, 1), pred_range, &alt) > 0) {
            shift_range(&alt, -pred_range.start);
            alt_len = range_len(alt);
        }

        Range ref_segment = scale_range(ref, NUM_SCORES);

        // Copy segments between variants/boundaries
        int segment_len = (ref_segment.start - ref_idx);
        memcpy(corrected_predictions + ref_idx, predictions + hap_idx, segment_len * sizeof(float));

        int offset = ref_len - alt_len; // Naive offset in end position vs ref and alt
        if (offset == 0) {
            memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, range_len(ref_segment) * sizeof(float));
        } else if (offset < 0) {
            float best_scores[NUM_SCORES] = {0};
            get_best_insertion_score(predictions + hap_idx + segment_len, alt_len, best_scores);
            memcpy(corrected_predictions + ref_segment.start, best_scores, NUM_SCORES * sizeof(float));
         } else { // offset > 0
            memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, NUM_SCORES * sizeof(float));
            memset(corrected_predictions + (ref_segment.start + NUM_SCORES), 0, offset * NUM_SCORES * sizeof(float));
        }

        hap_idx += segment_len + alt_len*NUM_SCORES;
        ref_idx = ref_segment.end;
    }

    memcpy(corrected_predictions + ref_idx, predictions + hap_idx, (num_reference_predictions - ref_idx) * sizeof(float));

    return corrected_predictions;
}

static void get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap) {
    const int sequence_offset = haplotype_offset(rs, gene);

    // Prepare kstring
    kstring_t haplotype_sequence = {0};
    haplotype_sequence.l = reference_sequence.l - sequence_offset;
    haplotype_sequence.m = haplotype_sequence.l + 1; 
    haplotype_sequence.s = malloc(haplotype_sequence.m);

    int ref_idx = 0;
    int hap_idx = 0;

    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        // Get reference position, intersected to gene boundaries to manage start overlap
        Range ref = variant_to_range(r->v, 0);
        Range intersected_ref = {0}; 
        if (intersect_ranges(ref, gene, &intersected_ref) == 0) continue;
        shift_range(&ref, -gene.start);
        shift_range(&intersected_ref, -gene.start);
        int ref_segment = intersected_ref.start - ref_idx;

        // Get alt position, intersected to gene boundaries
        Range alt = variant_to_range(r->v, 1);
        Range intersected_alt = {0}; 
        int alt_segment = 0, skipped_bases = 0;
        if (intersect_ranges(alt, gene, &intersected_alt)) {
            shift_range(&alt, -gene.start);
            shift_range(&intersected_alt, -gene.start);

            alt_segment = range_len(intersected_alt);
            skipped_bases = intersected_alt.start - alt.start; // Only ever non-zero in the case of a deletion overlapping start
        }

        // Copy over reference sequence without any variation
        strncpy(haplotype_sequence.s + hap_idx, reference_sequence.s + ref_idx, ref_segment);

        ref_idx += ref_segment;
        hap_idx += ref_segment;

        // Copy over variant
        strncpy(haplotype_sequence.s + hap_idx, r->v->d.allele[1] + skipped_bases, alt_segment);

        ref_idx += range_len(intersected_ref);
        hap_idx += alt_segment;
    }

    // Copy over remaining sequence from reference and complete sequence string
    strncpy(haplotype_sequence.s + hap_idx, reference_sequence.s + ref_idx, reference_sequence.l - ref_idx);
    haplotype_sequence.s[haplotype_sequence.l] = '\0';
 
    *out_hap = haplotype_sequence;
}

// int64_t get_position_in_spliced_haplotype(bcf1_t *v, Variants *vs, Transcript *t) {
//     int positions_skipped = 0;
//     int previous_exon_end = t->definition->beg;
//
//     for (int i = 0; i < t->definition->ncds; i++) {
//         int64_t start = t->definition->cds[i]->beg;
//         if (start > v->pos) break;
//         int64_t end = start + t->definition->cds[i]->len;
//
//         int64_t start_in_hap = get_position_in_haplotype(start, vs, t->definition->beg, t->definition->end);
//         int64_t end_in_hap = get_position_in_haplotype(end, vs, t->definition->beg, t->definition->end);
//
//         positions_skipped += start_in_hap - previous_exon_end; // Size of intron between current and previous exon
//
//         if (iregion_intersects((IRegion) {-1, v->pos, v->pos + strlen(v->d.allele[1])}, (IRegion) {-1, start, end})) {
//             return get_position_in_haplotype(v->pos, vs, t->definition->beg, t->definition->end) - t->definition->beg - positions_skipped;
//         }
//
//         previous_exon_end = end_in_hap;
//     }
//
//     return -1;
// }

// void process_phased_variants(Records *vs, HdrFormat *vep_fmt, gff_t *gff, faidx_t *fai, bcf_hdr_t *hdr, Model *models) {
//     kvec_t(int64_t) frameshifts_rescued; kv_init(frameshifts_rescued);
//
//
//     for (int j = 0; j < vs->n; j++) {
//         Record *r = kv_A(*vs, j);
//         bcf1_t *vv = kv_A(*vs, j)->v;
//         VepTranscripts *transcript_csqs = get_vep_transcript_consequences(vv, vep_fmt, hdr);
//
//         // Filter transcripts for overlap with current variant
//         Region variant_region = iregion_to_region(variant_to_iregion(vv), hdr);
//         Transcripts *variant_transcripts; find_region_transcripts(&variant_region, hdr, gff, fai, &variant_transcripts);
//
//         Transcript *v_ts[transcript_csqs->n];
//         Variants *v_vs[transcript_csqs->n];
//         Variants *v_vs2[transcript_csqs->n];
//         for (int k = 0; k < transcript_csqs->n; k++) {
//             v_ts[k] = NULL;
//             v_vs[k] = NULL;
//             v_vs2[k] = NULL;
//         }
//
//         for (int k = 0; k < variant_transcripts->n; k++) {
//             Transcript t = kv_A(*variant_transcripts, k);
//             const char *tid = gff_id2string(gff, 0, t.definition->id);
//
//             for (int l = 0; l < transcript_csqs->n; l++) {
//                 if (strcmp(tid, kv_A(*transcript_csqs, l).tr_id) == 0) {
//                     v_ts[l] = &variant_transcripts->a[k];
//
//                     // Find region for matched transcript
//                     IRegion tr_region = transcript_to_iregion(&t);
//
//                     // Create subset of variants overlapping the current transcript, excluding the current variant
//                     Variants *vs_j = malloc(sizeof(Variants)); kv_init(*vs_j);
//                     Variants *vs_j2 = malloc(sizeof(Variants)); kv_init(*vs_j2);
//                     for (int m = 0; m < vs->n; m++) {
//
//                         bcf1_t *t_v = kv_A(*vs, m)->v;
//                         if (iregion_intersects(tr_region, variant_to_iregion(t_v))) {
//                             kv_push(bcf1_t *, *vs_j2, t_v);
//                             if (m == j) continue;
//                             kv_push(bcf1_t *, *vs_j, t_v);
//                         }
//                     }
//
//                     /* if (vs_j->n > 0)*/ v_vs[l] = vs_j;
//                     /*if (vs_j->n > 0)*/ v_vs2[l] = vs_j2;
//                     break;
//                 }
//             }
//         }
//
//         // ConsequencePredictions *csq_preds; init_consequence_predictions(csq_preds, transcript_csqs->n);
//         for (int k = 0; k < transcript_csqs->n; k++) {
//             if (v_ts[k] == NULL) continue; // There is no coding transcript overlapping the current variant;
//             if (v_vs[k] == NULL) continue; // There is no difference between the haplotype and reference 
//
//             Transcript *t = v_ts[k];
//             Variants *vv_vs = v_vs[k];
//             Variants *vv_vs2 = v_vs2[k];
//
//             // Calculate offsets between ref, hap_ref and hap_alt
//             int hap_ref_offset = 0;
//             for (int l = 0; l < vv_vs->n; l++) hap_ref_offset += get_allele_offset_variant(kv_A(*vv_vs, l), t->definition->beg, t->definition->end);
//
//             kstring_t hap_ref; get_haplotype_sequence(t->seq, *vv_vs, hap_ref_offset, t->definition->beg, &hap_ref);
//
//             int alt_allele_offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//             int hap_alt_offset = hap_ref_offset + alt_allele_offset;
//             kstring_t hap_alt; get_haplotype_sequence(t->seq, *vv_vs2, hap_alt_offset, t->definition->beg, &hap_alt);
//
//
//             // Check vep consequences for triggers for rechecking consequences in haplotype context
//             int trigger = 0;
//             for (int l = 0; l < transcript_csqs->a[k].num_consequences; l++) {
//                 const char *consequence = transcript_csqs->a[k].consequences[l];
//                 if ((strcmp(consequence, SPLICE_DONOR_VARIANT) == 0) ||
//                     (strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0)) {
//                     trigger = SPLICE_TRIGGER;
//                     // log_warn("Splice trigger %s:%li = %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, consequence);
//                 } else if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) {
//                     trigger = FRAMESHIFT_TRIGGER;
//                     // log_warn("Frameshift trigger %s:%li = %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, consequence);
//                 } else if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) {
//                     trigger = STOP_TRIGGER;
//                 }
//             }
//
//             int is_rescued = 0;
//             if (trigger == SPLICE_TRIGGER) {
//
//                 // Run reference
//                 int32_t start = VARIANT_POS + (t->definition->beg - vv->pos);
//                 int32_t end = VARIANT_POS - (t->definition->end - (vv->pos+1)); // Open ended, so + 1
//                 char *ref_seq = t->seq.s - start;
//
//                 if (start < 0) start = 0;
//                 if (end < 0) end = 0;
//
//                 Range ref_boundary = { start, end };
//                 char *ref = pad_sequence(ref_seq, ref_boundary, WIDTH);
//                 float *ref_predictions; int num_ref_predictions;
//                 int ref_len = strlen(ref);
//                 predict_for_seq(models, ref, ref_len, t->definition->strand, &ref_predictions, &num_ref_predictions);
//                 free(ref);
//
//                 // Run haplotype reference
//                 int pred_offset = prediction_offset(vv, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t transcript_start = get_position_in_haplotype(t->definition->beg, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t transcript_end = get_position_in_haplotype(t->definition->end, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t vposition = get_position_in_haplotype(vv->pos, vv_vs, t->definition->beg, t->definition->end);
//                 int32_t hap_ref_start = VARIANT_POS + (transcript_start - vposition);
//                 int32_t hap_ref_end = VARIANT_POS - (transcript_end - vposition) + pred_offset;
//                 char *hap_ref_seq = hap_ref.s - hap_ref_start;
//                 int hap_ref_seq_len = strlen(hap_ref_seq);
//                 if (hap_ref_start < 0) hap_ref_start = 0;
//                 if (hap_ref_end < 0) hap_ref_end = 0;
//                 Range hap_ref_boundary = {hap_ref_start, hap_ref_end};
//                 char *hap_ref_padded = pad_sequence(hap_ref_seq, hap_ref_boundary, WIDTH + pred_offset);
//                 float *hap_ref_predictions; int num_hap_ref_predictions;
//                 int hap_ref_len = strlen(hap_ref_padded);
//                 predict_for_seq(models, hap_ref_padded, hap_ref_len, t->definition->strand, &hap_ref_predictions, &num_hap_ref_predictions);
//                 float *corrected_hap_ref_predictions = match_predictions_to_reference(hap_ref_predictions, vv, vv_vs, t, num_ref_predictions);
//                 Score score = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, ref_predictions, corrected_hap_ref_predictions, num_ref_predictions, CAI_DISTANCE);
//
//
//                 // // Run haplotype alt
//                 transcript_end += get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                 int hap_alt_pred_offset = prediction_offset(vv, vv_vs2, t->definition->beg, t->definition->end);
//                 int32_t hap_alt_start = VARIANT_POS + (transcript_start - vposition);
//                 int32_t hap_alt_end = VARIANT_POS - (transcript_end - vposition) + hap_alt_pred_offset;
//                 char *hap_alt_seq = hap_alt.s - hap_alt_start;
//                 int hap_alt_seq_len = strlen(hap_alt_seq);
//                 if (hap_alt_start < 0) hap_alt_start = 0;
//                 if (hap_alt_end < 0) hap_alt_end = 0;
//                 Range hap_alt_boundary = {hap_alt_start, hap_alt_end};
//                 char *hap_alt_padded = pad_sequence(hap_alt_seq, hap_alt_boundary, WIDTH + hap_alt_pred_offset); //+ get_allele_offset(vv, t->definition->beg, t->definition->end));
//                 float *hap_alt_predictions; int num_hap_alt_predictions;
//                 int hap_alt_len = strlen(hap_alt_padded);
//                 predict_for_seq(models, hap_alt_padded, hap_alt_len, t->definition->strand, &hap_alt_predictions, &num_hap_alt_predictions);
//                 float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, vv, vv_vs2, t, num_ref_predictions);
//                 Score score_alt = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, corrected_hap_ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);
//
//                 // Calculate scores between ref and hap_alt
//                 Score score_hap = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);
//
//
//                 if (r->c->hap_score[k]->l > 0) {
//                     kputc('&', r->c->hap_score[k]);
//                 }
//
//                 if (score_hap.dl > SPLICEAI_THRESHOLD && score_hap.dg > SPLICEAI_THRESHOLD && (score_hap.dg_idx - score_hap.dl_idx) % 3 == 0) is_rescued = 1;
//
//                 if (score_hap.al > SPLICEAI_THRESHOLD && score_hap.ag > SPLICEAI_THRESHOLD && (score_hap.ag_idx - score_hap.al_idx) % 3 == 0) is_rescued = 1;
//
//                 if (score_hap.al < SPLICEAI_RESCUE_THRESHOLD && score_hap.dl < SPLICEAI_RESCUE_THRESHOLD) is_rescued = 1;
//
//                 if (is_rescued) {
//                     log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//                 }
//
//                 ksprintf(r->c->hap_score[k], "%s|%s|%.2f|%.2f|%.2f|%.2f|%d|%d|%d|%d", score_hap.alt, score_hap.gene, score_hap.ag, score_hap.al, score_hap.dg, score_hap.dl, score_hap.ag_idx, score_hap.al_idx, score_hap.dg_idx, score_hap.dl_idx);
//                 log_info(r->c->hap_score[k]->s);
//
//             free(ref_predictions);
//             free(hap_ref_predictions);
//             free(corrected_hap_ref_predictions);
//             free(hap_alt_predictions);
//             free(corrected_hap_alt_predictions);
//
//             }
//
//             if (trigger) {
//
//                 kstring_t spliced_transcript = {0};
//                 // Splice transcript
//                 int length = 0;
//                 int hap_length = 0;
//                 for (int l = 0; l < t->definition->ncds; l++) {
//                     int64_t start = t->definition->cds[l]->beg;
//                     int64_t end = start + t->definition->cds[l]->len;
//                     int64_t start_in_hap = get_position_in_haplotype(start, vv_vs2, t->definition->beg, t->definition->end);
//                     int64_t end_in_hap = get_position_in_haplotype(end, vv_vs2, t->definition->beg, t->definition->end);
//
//                     length += (end - start);
//                     hap_length += (end_in_hap - start_in_hap);
//                     kputsn(hap_alt.s + (start_in_hap - t->definition->beg), end_in_hap - start_in_hap, &spliced_transcript);
//                 }
//
//                 int v_pos = get_position_in_spliced_haplotype(vv, vv_vs2, t);
//                 if (trigger == FRAMESHIFT_TRIGGER) {
//                     int offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                     if (v_pos == -1) continue;
//
//                     kvec_t(int64_t) frameshifts_used; kv_init(frameshifts_used);
//                     kv_push(int64_t, frameshifts_used, vv->pos);
//
//                     if (t->definition->strand == STRAND_FWD) {
//                         for (int l = 0; l < vv_vs->n; l++) {
//                             bcf1_t *vvv = kv_A(*vv_vs, l);
//
//                             int already_rescued = 0;
//                             for (int m = 0; m < frameshifts_rescued.n; m++) {
//                                 if (kv_A(frameshifts_rescued, m) == vvv->pos) {
//                                     already_rescued = 1;
//                                     break;
//                                 }
//                             }
//                             if (already_rescued) continue;
//
//                             int v2_pos = get_position_in_spliced_haplotype(vvv, vv_vs2, t);
//                             int variant_offset = get_allele_offset_variant(vvv, t->definition->beg, t->definition->end);
//                             if (variant_offset != 0 && v2_pos != -1 && v2_pos > v_pos && v2_pos <= v_pos + MAX_FRAMESHIFT_RADIUS) {
//                                 offset += variant_offset;
//                                 kv_push(int64_t, frameshifts_used, vvv->pos);
//                             }
//                             if (v2_pos > v_pos + MAX_FRAMESHIFT_RADIUS) break;
//                         }
//                     }
//                     else if (t->definition->strand == STRAND_REV) {
//                         for (int l = vv_vs->n-1; l >= 0; l--) {
//                             bcf1_t *vvv = kv_A(*vv_vs, l);
//
//                             int already_rescued = 0;
//                             for (int m = 0; m < frameshifts_rescued.n; m++) {
//                                 if (kv_A(frameshifts_rescued, m) == vvv->pos) {
//                                     already_rescued = 1;
//                                     break;
//                                 }
//                             }
//                             if (already_rescued) continue;
//
//                             int v2_pos = get_position_in_spliced_haplotype(vvv, vv_vs2, t);
//                             int variant_offset = get_allele_offset_variant(vvv, t->definition->beg, t->definition->end);
//                             if (variant_offset != 0 && v2_pos != -1 && v2_pos < v_pos && v2_pos >= v_pos - MAX_FRAMESHIFT_RADIUS) {
//                                 offset += variant_offset;
//                                 kv_push(int64_t, frameshifts_used, vvv->pos);
//                             }
//                             if (v2_pos < v_pos - MAX_FRAMESHIFT_RADIUS) break;
//                         }
//                     }
//
//                     if (offset % 3 == 0) {
//                         log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//
//                         is_rescued = 1;
//
//                         for (int l = 0; l < frameshifts_used.n; l++) {
//                             kv_push(int64_t, frameshifts_rescued, kv_A(frameshifts_used, l));
//                         }
//                     }
//
//                     kv_destroy(frameshifts_used);
//                 }
//
//                 if (trigger == STOP_TRIGGER) {
//                     int offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                     int start = v_pos / 3;
//                     int end = offset > 0 ? (v_pos + offset) / 3 + 1 : v_pos / 3 + 1;
//                     int stop_hit = 0;
//                     if (t->definition->strand == STRAND_FWD) {
//                         v_pos = get_position_in_spliced_haplotype(vv, vv_vs2, t);
//                         const char *stop_codons[3] = STOP_CODONS;
//                         for (int l = start; l < end; l++) {
//                             for (int m = 0; m < 3; m++) {
//                                 const char *stop = stop_codons[m];
//                                 char *frame = spliced_transcript.s + l * 3;
//                                 if (stop[0] == frame[0] && stop[1] == frame[1] && stop[2] == frame[2]) {
//                                     stop_hit = 1;
//                                     break;
//                                 }
//                             }
//                             if (stop_hit) break;
//                         }
//                     } else {
//                         stop_hit = 1;
//                     }
//                     // else if (t->definition->strand == STRAND_REV) {
//                     //     const char *stop_codons[3] = COMPLEMENTARY_STOP_CODONS;
//                     //     for (int l = end-1; l >= start; l--) {
//                     //         for (int m = 0; m < 3; m++) {
//                     //             const char *stop = stop_codons[m];
//                     //             char *frame = spliced_transcript.s + l * 3;
//                     //             if (stop[0] == frame[2] && stop[1] == frame[1] && stop[2] == frame[0]) {
//                     //                 stop_hit = 1;
//                     //                 break;
//                     //             }
//                     //         }
//                     //
//                     //         if (stop_hit) break;
//                     //     }
//                     // }
//
//                     if (!stop_hit) {
//                         log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//                         is_rescued = 1;
//                     }
//                 }
//
//
//                 free(spliced_transcript.s);
//             }
//
//             if (r->c->is_rescue[k]->l > 0) {
//                 kputc('&', r->c->is_rescue[k]);
//             }
//             kputw(is_rescued, r->c->is_rescue[k]);
//
//             if (r->c->rescue_type[k]->l > 0) {
//                 kputc('&', r->c->rescue_type[k]);
//             }
//             if (is_rescued) {
//                 if (trigger == SPLICE_TRIGGER) {
//                     kputs("rescued_splice", r->c->rescue_type[k]);
//                 } else if (trigger == FRAMESHIFT_TRIGGER) {
//                     kputs("rescued_frameshift", r->c->rescue_type[k]);
//                 } else {
//                     kputs("rescued_stop", r->c->rescue_type[k]);
//                 }
//             } else {
//                 kputc('.', r->c->rescue_type[k]);
//             }
//
//             if (trigger != SPLICE_TRIGGER) {
//                 if (r->c->hap_score[k]->l > 0) kputc('&', r->c->hap_score[k]);
//                 kputc('.', r->c->hap_score[k]);
//             }
//
//             // Free per transcript predicted for
//             free(hap_ref.s);
//             free(hap_alt.s);
//             kv_destroy(*vv_vs);
//             kv_destroy(*vv_vs2);
//             free(vv_vs);
//             free(vv_vs2);
//         }
//
//         // Free per variant predicted for
//         for (int i = 0; i < variant_transcripts->n; i++) { free(kv_A(*variant_transcripts, i).seq.s); }
//         kv_destroy(*variant_transcripts); free(variant_transcripts);
//
//         for (int i = 0; i < transcript_csqs->n; i++) { 
//             VepTranscript tcs = kv_A(*transcript_csqs, i);
//             free(tcs.tr_id);
//             for (int j = 0; j < transcript_csqs->a[i].num_consequences; j++) {
//                 free(transcript_csqs->a[i].consequences[j]);
//             }
//             free(tcs.consequences);
//         }
//         kv_destroy(*transcript_csqs);
//         free(transcript_csqs);
//
//     }
//
//     kv_destroy(frameshifts_rescued);
// }

static int add_rescueranger_tags(bcf_hdr_t *hdr) {
    const char *tags[7] = INFO_TAGS;
    const char *descriptions[7] = { SPLICEAI_HAP_DESC, CS_DESC, RESCUE_DESC, RESCUE_TYPE_DESC, PTV_DESC, PTV_TYPE_DESC, PREDICTION_FLAGS_DESC };

    for (int i = 0; i < 7; i++) {
        if (bcf_hdr_append(hdr, descriptions[i]) != 0) {
            log_error("Failed to append description for tag %s to vcf header.", tags[i]);
            return -1;
        }
    }

    return 0;
}

static void filter_to_most_likely_prediction(size_t num_probabilities, float **probabilities) {
    for (int i = 0; i < num_probabilities; i += NUM_SCORES) {
        float max = 0;
        int max_idx = 0;

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if ((*probabilities)[pred_idx] <= max) continue;

            max = (*probabilities)[pred_idx];
            max_idx = pred_idx;
        }

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if (pred_idx == max_idx) continue;

            (*probabilities)[pred_idx] = 0;
        }
    }
}

static int transform_to_posterior_probabilities(const int target, const size_t num_probabilities, const float *probabilities, float **posterior_probabilities) {
    int i = 0;
    float sum = 0;
    int num_values = 0;

    for (i = target; i < num_probabilities; i+=NUM_SCORES) {
        if (probabilities[i] == 0) continue;

        num_values++;
        sum += probabilities[i];
    }

    if (num_values == 0) return 0;

    for (i = target; i < num_probabilities; i+=NUM_SCORES) {
        if (probabilities[i] == 0) continue;

        (*posterior_probabilities)[i] = probabilities[i] / sum;
    }

    return num_values;
}

static int set_transcript_exons(Transcript *transcript, const char *chr, gff_t *gff) {
    regidx_t *exons = gff_get(gff, idx_exon);
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
        if (ci < transcript->ncds) cds = &transcript->cds[ci++];

        Range *exon = malloc(sizeof(Range));
        *exon = (Range) {gf_e->beg, gf_e->end};

        // Resize as needed
        if (transcript->nexons >= transcript->mexons) {
            transcript->mexons *= 2;
            transcript->exons = realloc(transcript->exons, transcript->mexons * sizeof(Exon *));
        }


        // Does this exon have a CDS?
        Range cds_intersection = {0};
        if (cds != NULL) intersect_ranges(*exon, *cds, &cds_intersection);

        transcript->exons[transcript->nexons] = malloc(sizeof(Exon));
        *transcript->exons[transcript->nexons] = (Exon) { cds, exon, NULL };
        transcript->nexons++;
    }


    regitr_destroy(itr);
    return 1;
}

static void shift_exons_to_haplotype_sequence(const Transcript *t, const Range gene, const Records *rs) {
    for (int i = 0; i < t->nexons; i++) {
        Exon *exon = t->exons[i];

        // Shift ranges to by start of gene to make them relative to gene start
        shift_range(exon->bounds, -gene.start);
        exon->bounds->start -= haplotype_offset(rs, (Range) { gene.start, exon->bounds->start });
        exon->bounds->end -= haplotype_offset(rs, (Range) { gene.start, exon->bounds->end });

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
            log_info("Cryptic Found!");
            return CRYPTIC_SPLICE_VALUE;  // Set cryptic splice trigger to 1 if values are greater than CRYPTIC_SPLICE_THRESHOLD
        }
    }

    return 0;
}

float * predict_on_sequence(Model *models, char *seq, uint32_t position, Range prediction_boundaries, int strand, Records rs) {
    int64_t gene_start = prediction_boundaries.start;
    int64_t gene_end = prediction_boundaries.end - haplotype_offset(&rs, prediction_boundaries);
    int64_t vposition = position - haplotype_offset(&rs, (Range) { prediction_boundaries.start, position });

    int32_t start_offset = VARIANT_POS + (gene_start - vposition);
    int32_t end_offset = VARIANT_POS - (gene_end - vposition);

    char *shifted_seq = seq - start_offset;
    if (start_offset < 0) start_offset = 0;
    if (end_offset < 0) end_offset = 0;

    Range padding_boundary = { start_offset, end_offset };
    char *padded_seq = pad_sequence(shifted_seq, padding_boundary, WIDTH);
    float *predictions; int num_predictions;
    predict_for_seq(models, padded_seq, WIDTH, strand, &predictions, &num_predictions);
    free(padded_seq);

    return predictions;
}

int create_splice_probabilities(int vposition, Transcripts *transcripts, int num_hap_alt_predictions, float *hap_alt_predictions, SpliceProbability **splice_probabilities) {
    float *posterior_probabilities = calloc(num_hap_alt_predictions, sizeof(float));
    int num_acceptors = transform_to_posterior_probabilities(ACCEPTOR_POS, num_hap_alt_predictions, hap_alt_predictions, &posterior_probabilities);
    int num_donors = transform_to_posterior_probabilities(DONOR_POS, num_hap_alt_predictions, hap_alt_predictions, &posterior_probabilities);

    int num_splice_probabilities = num_acceptors + num_donors;

    int canonicals[PRED_WIDTH] = {0};
    for (int ti = 0; ti < transcripts->n; ti++) {
        Transcript *t = kv_A(*transcripts, ti);

        if (t == NULL) continue;
        for (int ei = 0; ei < t->nexons; ei++) {
            Exon *e = t->exons[ei];

            if (t->strand == STRAND_UNK) continue;

            if (abs(vposition - (int32_t) e->bounds->start) <= CAI_DISTANCE) {
                int index = abs(vposition - (int32_t) e->bounds->start) + CAI_DISTANCE;

                if(canonicals[index] == 0) num_splice_probabilities++;
                canonicals[index] = ACCEPTOR;
            } else if (abs(vposition - (int32_t) e->bounds->end) < CAI_DISTANCE) {
                int index = abs(vposition - (int32_t) e->bounds->end) + CAI_DISTANCE;

                if(canonicals[index] == 0) num_splice_probabilities++;
                canonicals[index] = DONOR;
            }
        }
    }

    SpliceProbability *probabilities = malloc(num_splice_probabilities * sizeof(SpliceProbability));
    int prob_idx = 0;
    for (int i = 0; i < num_hap_alt_predictions; i++) {
        if (i % NUM_SCORES == 0) continue; // Skip predictions for neither

        float prior_probability = hap_alt_predictions[i];
        float posterior_probability = posterior_probabilities[i];
        if (canonicals[i / NUM_SCORES] == DONOR && i % NUM_SCORES == DONOR) {
            int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
            int prediction_position = vposition + offset_from_variant;
            probabilities[prob_idx++] = (SpliceProbability) { DONOR, prior_probability, posterior_probability, prediction_position };
            continue;
        } else if (canonicals[i / NUM_SCORES] == ACCEPTOR && i % NUM_SCORES == ACCEPTOR) {
            int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
            int prediction_position = vposition + offset_from_variant;
            probabilities[prob_idx++] = (SpliceProbability) { ACCEPTOR, prior_probability, posterior_probability, prediction_position };
            continue;
        }

        if (posterior_probability == 0) continue; // Skip null predictions

        int type = i % 3;
        int offset_from_variant = ((int) (i / NUM_SCORES)) - CAI_DISTANCE;
        int prediction_position = vposition + offset_from_variant;

        probabilities[prob_idx++] = (SpliceProbability) { type, prior_probability, posterior_probability, prediction_position };
    }

    *splice_probabilities = probabilities;
    return num_splice_probabilities;
}

void resize_splice_sites(SpliceSite **sites_array, size_t *current_size) {
    *current_size *= 2;
    *sites_array = realloc(*sites_array, (*current_size) * sizeof(SpliceSite));
}

int check_proximity(uint32_t predicted_position, uint32_t canonical_site_position) {
    return abs((int32_t)predicted_position - (int32_t)canonical_site_position) <= PROXIMITY_THRESHOLD;
}

void store_splice_sites(SpliceProbability prediction, int *predict_for_pair, int *predict_for_competitor, SpliceSites *splice_sites) {
    if (prediction.type == DONOR & prediction.posterior_probability > 0) {
        if (splice_sites->canonical_donors[0].gff_canonical && prediction.posterior_probability && prediction.position == splice_sites->canonical_donors[0].position) { // If canonical exists and is the only splice site within the window
            splice_sites->canonical_donors[0].predicted = 1;
            splice_sites->canonical_donors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_donors[0].posterior_probability = prediction.posterior_probability;
        } else if (splice_sites->canonical_donors[0].gff_canonical && check_proximity(prediction.position, splice_sites->canonical_donors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical donor site
            int canonical_donor_exists = 0;
            int canonical_donor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_donors; i++) {
                if (splice_sites->canonical_donors[i].position == prediction.position) {
                    canonical_donor_exists = 1;
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
                        splice_sites->canonical_donors[i].predicted = 1;
                        splice_sites->canonical_donors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_donors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_donors >= splice_sites->m_canonical_donors) {
                    resize_splice_sites(&splice_sites->canonical_donors, &splice_sites->m_canonical_donors);
                }
                SpliceSite site = { prediction.position, 1, 0, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_donors[splice_sites->n_canonical_donors++] = site;
            }
        // Store cryptic donor site
        } else {
            // Check if prediction position already exists as a cryptic donor site
            int cryptic_donor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_donors; i++) {
                if (splice_sites->cryptic_donors[i].position == prediction.position) {
                    cryptic_donor_exists = 1;
                }
            }

            if (!cryptic_donor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_donors >= splice_sites->m_cryptic_donors) {
                    resize_splice_sites(&splice_sites->cryptic_donors, &splice_sites->m_cryptic_donors);
                }
                SpliceSite site = { prediction.position, 1, 0, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_donors[splice_sites->n_cryptic_donors++] = site;
            }
        }
    } else if (prediction.type == ACCEPTOR & prediction.posterior_probability > 0) {
        // Store canonical acceptor site
        if (splice_sites->canonical_acceptors[0].gff_canonical && prediction.posterior_probability && prediction.position == splice_sites->canonical_acceptors[0].position) { // If canonical exists and is the only splice site within the window
            splice_sites->canonical_acceptors[0].predicted = 1;
            splice_sites->canonical_acceptors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_acceptors[0].posterior_probability = prediction.posterior_probability;
        } else if (splice_sites->canonical_acceptors[0].gff_canonical & check_proximity(prediction.position, splice_sites->canonical_acceptors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical acceptor site
            int canonical_acceptor_exists = 0;
            int canonical_acceptor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_acceptors; i++) {
                if (splice_sites->canonical_acceptors[i].position == prediction.position) {
                    canonical_acceptor_exists = 1;
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
                        splice_sites->canonical_acceptors[i].predicted = 1;
                        splice_sites->canonical_acceptors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_acceptors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_acceptors >= splice_sites->m_canonical_acceptors) {
                    resize_splice_sites(&splice_sites->canonical_acceptors, &splice_sites->m_canonical_acceptors);
                }
                SpliceSite site = { prediction.position, 1, 0, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_acceptors[splice_sites->n_canonical_acceptors++] = site;
            }
        // Store cryptic acceptor site
        } else {
            // Check if prediction position already exists as a cryptic acceptor site
            int cryptic_acceptor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_acceptors; i++) {
                if (splice_sites->cryptic_acceptors[i].position == prediction.position) {
                    cryptic_acceptor_exists = 1;
                }
            }

            if (!cryptic_acceptor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_acceptors >= splice_sites->m_cryptic_acceptors) {
                    resize_splice_sites(&splice_sites->cryptic_acceptors, &splice_sites->m_cryptic_acceptors);
                }
                SpliceSite site = { prediction.position, 1, 0, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_acceptors[splice_sites->n_cryptic_acceptors++] = site;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    log_set_level(2);
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);

    // log_debug("Loading SpliceAI models...");
    Model *models = load_models(argv[1]);

    // Parse gff file
    const char *gff_path = argv[2];
    log_debug("Loading GFF from: %s", gff_path);
    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);

    // Open fasta file
    const char *fasta_path = argv[3];
    log_debug("Loading FASTA from: %s", fasta_path);
    faidx_t *fai = fai_load(fasta_path);

    // Open VCF
    const char *input_vcf_path = argv[4];
    log_debug("Loading VCF from: %s", input_vcf_path);
    htsFile *fp_in = hts_open(input_vcf_path, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp_in);
    HdrFormat *vep_fmt = parse_csq_hdr_fmt(hdr);

    const char *output_vcf_path = argv[5];
    htsFile *fp_out = hts_open(output_vcf_path, "w");
    add_rescueranger_tags(hdr);
    bcf_hdr_write(fp_out, hdr);

    // const char *raw_score_path = "raw_scores.tsv";
    // FILE *raw_scores = fopen(raw_score_path, "w");

    bcf1_t *v = bcf_init();

    kvec_t(RecordBlock) batches; kv_init(batches);

    Records *batch = malloc(sizeof(Records));
    kv_init(*batch);

    IRegion batch_region = {-1, -1, -1};

    // Read records from vcf and store into coding/non-coding record blocks
    while (bcf_read(fp_in, hdr, v) >= 0) {
        bcf_unpack(v, BCF_UN_STR);

        IRegion var_region = {0};
        int is_non_coding = find_gene_region_for_variant(v, hdr, gff, &var_region);
        if (batch_region.rid == -1) batch_region = var_region;
        int overlaps_batch = is_non_coding ? 0 : iregion_intersects(batch_region, var_region);

        if (batch->n > 0 && !overlaps_batch) { // Variant is either non-coding or falls outside of batch
            RecordBlock block = {kv_A(*batch, 0)->vep == NULL ? 0 : 1, batch};
            kv_push(RecordBlock, batches, block);

            batch_region = var_region;
            batch = malloc(sizeof(Records));
            kv_init(*batch);
        }

        VepCsqs *vep_csqs = is_non_coding ? NULL : parse_vep_transcript_consequences(v, vep_fmt, hdr);
        if (overlaps_batch && var_region.end > batch_region.end) batch_region.end = var_region.end;
        kv_push(Record *, *batch, init_record(v, vep_csqs));
    }

    if (batch->n > 0) {
        RecordBlock block = {kv_A(*batch, 0)->vep == NULL ? 0 : 1, batch};
        kv_push(RecordBlock, batches, block);
    }

    // Convert GFF structures to custom format for ease of use and manipulability
    khash_t(GeneHash) *genes = kh_init_GeneHash();
    khash_t(TranscriptHash) *transcripts = kh_init_TranscriptHash();

    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        // Find matching gff transcripts/genes for records in batch
        for (int i = 0; i < block.records->n; i++) {
            Record *r = kv_A(*block.records, i);
            gf_tscripts *tscripts = NULL;
            find_gf_tscripts_for_record(r, hdr, gff, &tscripts);
            if (tscripts == NULL) continue; // A variant in a coding gene can still have no overlapping coding transcripts

            // Convert tscripts to Transcript
            r->tscripts = malloc(sizeof(Transcripts));
            kv_init(*r->tscripts);
            for (int j = 0; j < tscripts->n; j++) {
                gf_tscript_t *gf_t = kv_A(*tscripts, j);
                if (gf_t == NULL) { kv_push(Transcript *, *r->tscripts, NULL); continue; }

                int key = gf_t->id;
                khint32_t hash_index = kh_get_TranscriptHash(transcripts, key);
                if (hash_index != kh_end(transcripts)) {
                    Transcript *t = kh_val(transcripts, hash_index);
                    kv_push(Transcript *, *r->tscripts, t);
                    continue;
                }

                Transcript *t = malloc(sizeof(Transcript));
                t->id = key;
                t->start = gf_t->beg;
                t->end = gf_t->end;
                t->strand = gf_t->strand;

                t->mcds = gf_t->mcds;
                t->ncds = gf_t->ncds;
                t->cds = malloc(t->mcds * sizeof(Range));
                for (int k = 0; k < t->ncds; k++) {
                    uint32_t start = gf_t->cds[k]->beg;
                    uint32_t end = gf_t->cds[k]->beg + gf_t->cds[k]->len;
                    t->cds[k] = (Range) { start, end };
                }

                const char *chr = bcf_hdr_id2name(hdr, r->v->rid);
                set_transcript_exons(t, chr, gff);

                int is_new;
                hash_index = kh_put_TranscriptHash(transcripts, key, &is_new);
                kh_val(transcripts, hash_index) = t;

                kv_push(Transcript *, *r->tscripts, t);
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
                char *seq = faidx_fetch_seq(fai, bcf_hdr_id2name(hdr, r->v->rid), gene_range.start, gene_range.end, &slen);
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

            r->complete_scores = malloc(r->genes->n * NUM_PHASES * sizeof(Score *));
            for (int i = 0; i < r->genes->n * NUM_PHASES; i++) r->complete_scores[i] = NULL;

            r->cs = malloc(r->genes->n * NUM_PHASES * sizeof(SpliceChange *)); 
            for (int i = 0; i < r->genes->n * NUM_PHASES; i++) r->cs[i] = NULL;
        }

    }

    bcf_destroy(v);
    destroy_hdr_fmt(vep_fmt);
    bcf_close(fp_in);
    fai_destroy(fai);

    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        for (int i = 0; i < block.records->n; i++) {
            Record *r = kv_A(*block.records, i);

            // Check VEP consequences for rescue triggers
            for (int j = 0; j < r->vep->n; j++) {
                TscriptCsqs transcript_consequence = kv_A(*r->vep, j);

                int is_non_coding = 0;
                for (int k = 0; k < transcript_consequence.num_consequences; k++) {
                    const char *consequence = transcript_consequence.consequences[k];
                    if (strcmp(consequence, "non_coding_transcript_variant") == 0 ||
                        strcmp(consequence, "NMD_transcript_variant") == 0) {
                        is_non_coding = 1;
                    }
                }

                if (is_non_coding) continue;

                for (int k = 0; k < transcript_consequence.num_consequences; k++) {
                    const char *consequence = transcript_consequence.consequences[k];

                    if (strcmp(consequence, SPLICE_DONOR_VARIANT) == 0 ||
                        strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0) {
                        r->triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
                    }
                    // if (strcmp(consequence, SPLICE_REGION_VARIANT) == 0) {
                    //     r->triggers[SPLICE_TRIGGER] = 1;
                    // }
                    if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) {
                        r->triggers[FRAMESHIFT_TRIGGER] = 1;
                    }
                    if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) {
                        r->triggers[STOP_TRIGGER] = 1;
                    }
                }
            }

            if (r->triggers[SPLICE_TRIGGER] == 0) r->triggers[SPLICE_TRIGGER] = trigger_cryptic_splice(r->v, hdr);
        }
    }

    // Split records by phase
    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        PhasedRecords *phased_records = split_records_by_phase(block.records, hdr);
        for (int p = 0; p < NUM_PHASES; p++) {
            Records *batch = kv_A(*phased_records, p);
            remove_unworkable_variants(batch, hdr);

            for (int i = 0; i < batch->n; i++) {
                Record *r = kv_A(*batch, i);
                // log_info("%d:%d Initialize splice sites", r->v->rid, r->v->pos);

                for (int ti = 0; ti < r->tscripts->n; ti++) {
                    Transcript *t = kv_A(*r->tscripts, ti);
                    Gene *gene = kv_A(*r->transcript_genes, ti);
                    if (t == NULL || gene == NULL) continue;


                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;

                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    for (int j = 0; j < batch->n; j++) {
                        Record *hap_r = kv_A(*batch, j);
                        Range record_range = variant_to_range(r->v, 0);

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
                        Range canonical_start = (Range) { e->bounds->start, e->bounds->start+1 };
                        Range canonical_end = (Range) { e->bounds->end, e->bounds->end+1 };
                        for (int ri = 0; ri < hap_in.n; ri++) {
                            Record *hap_r = kv_A(hap_in, ri);
                            int is_deletion = allele_offset(hap_r->v, 0, 1, gene->coords) > 0;
                            if (!is_deletion) continue;

                            Range ref = variant_to_range(r->v, 0);

                            Range intersection;
                            if (intersect_ranges(ref, canonical_start, &intersection) && range_len(intersection) > 0) kept_start = 0;
                            if (intersect_ranges(ref, canonical_end, &intersection) && range_len(intersection) > 0) kept_end = 0;
                        }

                        uint32_t exon_hap_start = e->bounds->start - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->bounds->start });
                        uint32_t exon_hap_end = e->bounds->end - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->bounds->end });

                        SpliceSites *sites = malloc(sizeof(SpliceSites));
                        init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
                        init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
                        init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
                        init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

                        if (t->strand == STRAND_FWD) {
                            if (ei > 0) sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { exon_hap_start, 0, kept_start, 1.0, 1.0 };
                            if (ei < t->nexons-1) sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { exon_hap_end, 0, kept_end, 1.0, 1.0 };

                            *e->bounds = (Range) {exon_hap_start, exon_hap_end};
                        } else if (t->strand == STRAND_REV) {
                            if (ei < t->nexons-1) sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { exon_hap_end, 0, kept_end, 1.0, 1.0 };
                            if (ei > 0) sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { exon_hap_start, 0, kept_start, 1.0, 1.0 };

                            *e->bounds = (Range) {exon_hap_end, exon_hap_start};
                        }

                        e->splice_sites = sites;

                        // Convert CDS boundaries to haplotype positions
                        if (e->cds == NULL) continue;
                        uint32_t cds_hap_start = e->cds->start - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->cds->start });
                        uint32_t cds_hap_end = e->cds->end - haplotype_offset(&hap_in, (Range) { gene->coords.start, e->cds->end });
                        if (t->strand == STRAND_FWD) *e->cds = (Range) {cds_hap_start, cds_hap_end};
                        if (t->strand == STRAND_REV) *e->cds = (Range) {cds_hap_end, cds_hap_start};
                    }

                    kv_destroy(hap_in);
                }
            }

            for (int i = 0; i < batch->n; i++) {
                Record *r = kv_A(*batch, i);

                if (r->triggers[SPLICE_TRIGGER] == 0) continue;

                if (strlen(r->v->d.allele[1]) > CAI_DISTANCE) {
                    log_warn("Skipping on %s:%li, insertion is too long", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos);
                    continue;
                }
                if (r->v->rlen > CAI_DISTANCE) {
                    log_warn("Skipping prediction on %s:%li, deletion is too long", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos);
                    continue;
                }

                for (int g = 0; g < r->genes->n; g++) {
                    Gene *gene = kv_A(*r->genes, g);
                    if (gene == NULL) continue;

                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;

                    // Find records for the current allele
                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    Records hap_ex; kv_init(hap_ex); // Variant for current allele excluding current variant

                    for (int j = 0; j < batch->n; j++) {
                        Record *hap_r = kv_A(*batch, j);
                        Range record_range = variant_to_range(r->v, 0);

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

                    log_info("%d:%d", r->v->rid, r->v->pos);

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

                    Range ref_boundary = { start, end };
                    char *ref_padded = pad_sequence(ref_seq, ref_boundary, WIDTH);
                    float *ref_predictions; int num_ref_predictions;
                    int ref_len = strlen(ref_padded);
                    predict_for_seq(models, ref_padded, ref_len, gene->strand, &ref_predictions, &num_ref_predictions);
                    free(ref_padded);

                    // HAP REF
                    int pred_offset = prediction_offset(r->v, &hap_ex);

                    int64_t gene_start = gene->coords.start;
                    int64_t gene_end = gene->coords.end - haplotype_offset(&hap_ex, gene->coords);
                    int64_t vposition = r->v->pos - haplotype_offset(&hap_ex, (Range) {gene->coords.start, r->v->pos});

                    int32_t hap_ref_start = VARIANT_POS + (gene_start - vposition);
                    int32_t hap_ref_end = VARIANT_POS - (gene_end - vposition) - pred_offset;

                    char *hap_ref_seq = hap_ref.s - hap_ref_start;
                    if (hap_ref_start < 0) hap_ref_start = 0;
                    if (hap_ref_end < 0) hap_ref_end = 0;

                    Range hap_ref_boundary = {hap_ref_start, hap_ref_end};
                    char *hap_ref_padded = pad_sequence(hap_ref_seq, hap_ref_boundary, WIDTH - pred_offset);
                    float *hap_ref_predictions; int num_hap_ref_predictions;
                    int hap_ref_len = strlen(hap_ref_padded);
                    predict_for_seq(models, hap_ref_padded, hap_ref_len, gene->strand, &hap_ref_predictions, &num_hap_ref_predictions);
                    float *corrected_hap_ref_predictions = match_predictions_to_reference(hap_ref_predictions, r->v, &hap_ex, num_ref_predictions);

                    free(hap_ref_padded);

                    // HAP ALT
                    pred_offset = prediction_offset(r->v, &hap_in);
                    gene_end -= allele_offset(r->v, 0, 1, gene->coords);
                    int hap_alt_pred_offset = prediction_offset(r->v, &hap_in);

                    int32_t hap_alt_start = VARIANT_POS + (gene_start - vposition);
                    int32_t hap_alt_end = VARIANT_POS - (gene_end - vposition) - pred_offset;

                    char *hap_alt_seq = hap_alt.s - hap_alt_start;
                    if (hap_alt_start < 0) hap_alt_start = 0;
                    if (hap_alt_end < 0) hap_alt_end = 0;

                    Range hap_alt_boundary = {hap_alt_start, hap_alt_end};
                    char *hap_alt_padded = pad_sequence(hap_alt_seq, hap_alt_boundary, WIDTH - pred_offset);
                    float *hap_alt_predictions; int num_hap_alt_predictions;
                    int hap_alt_len = strlen(hap_alt_padded);
                    predict_for_seq(models, hap_alt_padded, hap_alt_len, gene->strand, &hap_alt_predictions, &num_hap_alt_predictions);
                    float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, r->v, &hap_in, num_ref_predictions);
                    free(hap_alt_padded);

                    // Calculate scores between ref and hap_alt
                    Score score = calculate_delta_scores(r->v->d.allele[1], gene->name, ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);

                    r->cs[g * NUM_PHASES + p] = malloc(sizeof(SpliceChange));
                    int ag_position = (score.ag_idx + CAI_DISTANCE) * NUM_SCORES + ACCEPTOR_POS;
                    int al_position = (score.al_idx + CAI_DISTANCE) * NUM_SCORES + ACCEPTOR_POS;
                    int dg_position = (score.dg_idx + CAI_DISTANCE) * NUM_SCORES + DONOR_POS;
                    int dl_position = (score.dl_idx + CAI_DISTANCE) * NUM_SCORES + DONOR_POS;
                    float cs_ag = corrected_hap_alt_predictions[ag_position] - corrected_hap_ref_predictions[ag_position];
                    float cs_al = corrected_hap_alt_predictions[al_position] - corrected_hap_ref_predictions[al_position];
                    float cs_dg = corrected_hap_alt_predictions[dg_position] - corrected_hap_ref_predictions[dg_position];
                    float cs_dl = corrected_hap_alt_predictions[dl_position] - corrected_hap_ref_predictions[dl_position];
                    *r->cs[g * NUM_PHASES + p] = (SpliceChange) {cs_ag, cs_al, cs_dg, cs_dl, fabs(cs_ag) + fabs(cs_al) + fabs(cs_dg) + fabs(cs_dl) };

                    r->complete_scores[g * NUM_PHASES + p] = malloc(sizeof(Score));
                    *r->complete_scores[g * NUM_PHASES + p] = score;

                    // Remove predictions that are most likely to be neutral
                    filter_to_most_likely_prediction(num_hap_alt_predictions, &hap_alt_predictions);

                    // Scale donor/acceptor probabilities by their prior
                    SpliceProbability *splice_probabilities;
                    int num_splice_probabilities = create_splice_probabilities(vposition, r->tscripts, num_hap_alt_predictions, hap_alt_predictions, &splice_probabilities);

                    for (int i = 0; i < r->tscripts->n; i++) {
                        Transcript *t = kv_A(*r->tscripts, i);
                        Gene *tgene = kv_A(*r->transcript_genes, i);

                        if (t == NULL || tgene->gff_id != gene->gff_id) continue;
                        for (int si = 0; si < num_splice_probabilities; si++) {
                            SpliceProbability prediction = splice_probabilities[si];
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
                                    positions = (Range) { e->bounds->start + 1, t->exons[ei + direction]->bounds->start };
                                } else { // ACCEPTOR
                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                    positions = (Range) { t->exons[ei - direction]->bounds->end, e->bounds->end - 1 };
                                }

                                Range intersection;
                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                int predict_for_pair = 0;
                                int predict_for_competitor = 1;
                                store_splice_sites(prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);

                                if (predict_for_pair) {
                                    if (prediction.type == DONOR) {
                                        int pi = ei + direction;
                                        if (pi < 0 || pi >= t->nexons) continue;
                                        SpliceSite pair = t->exons[pi]->splice_sites->canonical_acceptors[0];

                                        if (pair.predicted) continue;

                                        float *pair_predictions = predict_on_sequence(models, hap_alt.s, pair.position, tgene->coords, tgene->strand, hap_in);
                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, &pair_predictions);

                                        SpliceProbability *pair_probabilities;
                                        int num_pair_probabilities = create_splice_probabilities(pair.position, r->tscripts, PRED_WIDTH, pair_predictions, &pair_probabilities);
                                        for (int pi = 0; pi < num_pair_probabilities; pi++) {
                                            SpliceProbability pair_prediction = pair_probabilities[pi];

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->bounds->start + 1, t->exons[ei + direction]->bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->bounds->end, e->bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                store_splice_sites(pair_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    } else { // ACCEPTOR
                                        int pi = ei - direction;
                                        if (pi < 0 || pi >= t->nexons) continue;
                                        SpliceSite pair = t->exons[pi]->splice_sites->canonical_donors[0];

                                        if (pair.predicted) continue;

                                        float *pair_predictions = predict_on_sequence(models, hap_alt.s, pair.position, tgene->coords, tgene->strand, hap_in);
                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, &pair_predictions);

                                        SpliceProbability *pair_probabilities;
                                        int num_pair_probabilities = create_splice_probabilities(pair.position, r->tscripts, PRED_WIDTH, pair_predictions, &pair_probabilities);
                                        for (int pi = 0; pi < num_pair_probabilities; pi++) {
                                            SpliceProbability pair_prediction = pair_probabilities[pi];

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->bounds->start + 1, t->exons[ei + direction]->bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->bounds->end, e->bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

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

                                        if (competitor.predicted) continue;

                                        float *competitor_predictions = predict_on_sequence(models, hap_alt.s, competitor.position, tgene->coords, tgene->strand, hap_in);
                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, &competitor_predictions);

                                        SpliceProbability *competitor_probabilities;
                                        int num_competitor_probabilities = create_splice_probabilities(competitor.position, r->tscripts, PRED_WIDTH, competitor_predictions, &competitor_probabilities);
                                        for (int ci = 0; ci < num_competitor_probabilities; ci++) {
                                            SpliceProbability competitor_prediction = competitor_probabilities[ci];

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->bounds->start + 1, t->exons[ei + direction]->bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->bounds->end, e->bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                store_splice_sites(competitor_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    } else { // ACCEPTOR
                                        int ci = ei + direction;
                                        if (ci < 0 || ci >= t->nexons) continue;
                                        SpliceSite competitor = t->exons[ci]->splice_sites->canonical_acceptors[0];

                                        if (competitor.predicted) continue;

                                        float *competitor_predictions = predict_on_sequence(models, hap_alt.s, competitor.position, tgene->coords, tgene->strand, hap_in);
                                        filter_to_most_likely_prediction(PRED_WIDTH * NUM_SCORES, &competitor_predictions);

                                        SpliceProbability *competitor_probabilities;
                                        int num_competitor_probabilities = create_splice_probabilities(competitor.position, r->tscripts, PRED_WIDTH, competitor_predictions, &competitor_probabilities);
                                        for (int ci = 0; ci < num_competitor_probabilities; ci++) {
                                            SpliceProbability competitor_prediction = competitor_probabilities[ci];

                                            for (int ei = 0; ei < t->nexons; ei++) {
                                                Exon *e = t->exons[ei];

                                                Range positions;
                                                if (prediction.type == DONOR) {
                                                    if (ei == t->nexons-1 && t->strand == STRAND_FWD) continue; // Last exon, cannot form a pair with a downstream acceptor
                                                    if (ei == 0 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { e->bounds->start + 1, t->exons[ei + direction]->bounds->start };
                                                } else { // ACCEPTOR
                                                    if (ei == 0 && t->strand == STRAND_FWD) continue;
                                                    if (ei == t->nexons-1 && t->strand == STRAND_REV) continue;
                                                    positions = (Range) { t->exons[ei - direction]->bounds->end, e->bounds->end - 1 };
                                                }

                                                Range intersection;
                                                if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) continue;

                                                store_splice_sites(competitor_prediction, &predict_for_pair, &predict_for_competitor, e->splice_sites);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    free(ref_predictions);
                    free(hap_ref_predictions);
                    free(corrected_hap_ref_predictions);
                    free(hap_alt_predictions);
                    free(corrected_hap_alt_predictions);

                    free(hap_ref.s);
                    free(hap_alt.s);

                    kv_destroy(hap_in);
                    kv_destroy(hap_ex);
                }
            }

            for (int i = 0; i < batch->n; i++) {
                Record *r = kv_A(*batch, i);

                if (!r->triggers[SPLICE_TRIGGER] && !r->triggers[FRAMESHIFT_TRIGGER] && !r->triggers[STOP_TRIGGER]) continue;

                for (int ti = 0; ti < r->tscripts->n; ti++) {
                    Transcript *t = kv_A(*r->tscripts, ti);
                    Gene *gene = kv_A(*r->transcript_genes, ti);
                    if (t == NULL || gene == NULL) continue;

                    Range intersection;
                    intersect_ranges(gene->coords, variant_to_range(r->v, 0), &intersection);
                    if (range_len(intersection) == 0) continue;


                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    for (int j = 0; j < batch->n; j++) {
                        Record *hap_r = kv_A(*batch, j);
                        Range record_range = variant_to_range(r->v, 0);

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
                    int n_nc_regions = 0, m_nc_regions = t->nexons+5; // TODO: rework this, +5 is messy, and probably not needed in the vast majority of cases
                    Ranges *nc_regions = malloc(m_nc_regions * sizeof(Ranges));
                    for (int i2 = 0; i2 < m_nc_regions; i2++) {
                        kv_init(nc_regions[i2]);
                    }

                    int is_coding = 0;
                    if (t->strand == STRAND_FWD) {
                        for (int ei = 0; ei < t->nexons; ei++) {
                            Exon *e = t->exons[ei];
                            if (e->cds == NULL) continue;
                            is_coding = 1;

                            if (n_nc_regions == 0 && range_len(*e->bounds) != range_len(*e->cds) && e->bounds->start != e->cds->start) { // 5' UTR
                                Range ncr = { gene->coords.start, e->cds->start };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                            }

                            if (n_nc_regions > 0 && range_len(*e->bounds) != range_len(*e->cds) && e->bounds->end != e->cds->end && e->bounds->end+1 != e->cds->end && e->bounds->end+2 != e->cds->end) { // 3' UTR
                                Range ncr = { e->cds->end,  gene->coords.end };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                break; // End of transcript
                            }

                            if (e->splice_sites->canonical_donors[0].gff_canonical == 0 && e->splice_sites->n_cryptic_donors == 0) {
                                Range ncr = { e->cds->start, e->cds->start };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                continue;
                            } // No donors, skip

                            int lost_canonical_acceptor = e->splice_sites->canonical_acceptors[0].gff_canonical == 0 && e->splice_sites->n_cryptic_acceptors == 0;

                            // Get the exon in which the search for the acceptor will be performed
                            Exon *target = NULL;
                            if (lost_canonical_acceptor) {
                                if (ei+2 < t->nexons) target = t->exons[ei+2];
                            } else {
                                if (ei+1 < t->nexons) target = t->exons[ei+1];
                            }

                            // TODO: Add a way to handle missing canonical donor sites, and the possible exon skipped that would be induced

                            for (int di = 0; di < e->splice_sites->n_canonical_donors; di++) { // Always at least 1 canonical donor, and 1 canonical acceptor
                                uint32_t donor_pos = e->splice_sites->canonical_donors[di].position;

                                if (target == NULL) {
                                    Range ncr = { donor_pos, donor_pos }; // No splicing happens for this donor
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                    continue;
                                }

                                for (int ai = 0; ai < target->splice_sites->n_canonical_acceptors; ai++) {
                                    if (ai == 0 && !e->splice_sites->canonical_acceptors[0].gff_canonical) continue; // Canonical acceptor is lost, can't be used
                                    Range ncr = { donor_pos, target->splice_sites->canonical_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }

                                for (int ai = 0; ai < target->splice_sites->n_cryptic_acceptors; ai++) {
                                    Range ncr = { donor_pos, target->splice_sites->cryptic_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            for (int di = 0; di < e->splice_sites->n_cryptic_donors; di++) {
                                uint32_t donor_pos = e->splice_sites->cryptic_donors[di].position;

                                if (target == NULL) {
                                    Range ncr = { donor_pos, donor_pos }; // No splicing happens for this donor
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                    continue;
                                }

                                for (int ai = 0; ai < target->splice_sites->n_canonical_acceptors; ai++) {
                                    if (ai == 0 && !e->splice_sites->canonical_acceptors[0].gff_canonical) continue; // Canonical acceptor is lost, can't be used
                                    Range ncr = { donor_pos, target->splice_sites->canonical_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }

                                for (int ai = 0; ai < target->splice_sites->n_cryptic_acceptors; ai++) {
                                    Range ncr = { donor_pos, target->splice_sites->cryptic_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            n_nc_regions++;
                        }
                    } else if (t->strand == STRAND_REV) {
                        for (int ei = t->nexons-1; ei >= 0; ei--) {
                            Exon *e = t->exons[ei];
                            if (e->cds == NULL) continue;
                            is_coding = 1;

                            if (n_nc_regions == 0 && range_len(*e->bounds) != range_len(*e->cds) && e->bounds->start != e->cds->start) { // 5' UTR
                                Range ncr = { gene->coords.end, e->cds->start };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                            }
                            if (n_nc_regions > 0 && range_len(*e->bounds) != range_len(*e->cds) && e->bounds->end != e->cds->end && e->bounds->end+1 != e->cds->end && e->bounds->end+2 != e->cds->start) { // 3' UTR
                                Range ncr = { e->cds->end,  gene->coords.start };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                break; // End of transcript
                            }

                            if (e->splice_sites->canonical_donors[0].gff_canonical == 0 && e->splice_sites->n_cryptic_donors == 0) {
                                Range ncr = { e->cds->start, e->cds->start };
                                kv_push(Range, nc_regions[n_nc_regions], ncr);
                                n_nc_regions++;
                                continue;
                            } // No donors, intron_retention, check TODO for exon skipping

                            int lost_canonical_acceptor = e->splice_sites->canonical_acceptors[0].gff_canonical == 0 && e->splice_sites->n_cryptic_acceptors == 0;

                            // Get the exon in which the search for the acceptor will be performed
                            Exon *target = NULL;
                            if (lost_canonical_acceptor) {
                                if (ei-2 >= 0) target = t->exons[ei-2];
                            } else {
                                if (ei-1 >= 0) target = t->exons[ei-1];
                            }

                            // TODO: Add a way to handle missing canonical donor sites, and the possible exon skipped that would be induced

                            for (int di = 0; di < e->splice_sites->n_canonical_donors; di++) { // Always at least 1 canonical donor, and 1 canonical acceptor
                                uint32_t donor_pos = e->splice_sites->canonical_donors[di].position;

                                if (target == NULL) {
                                    Range ncr = { donor_pos, donor_pos }; // No splicing happens for this donor
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                    continue;
                                }

                                for (int ai = 0; ai < target->splice_sites->n_canonical_acceptors; ai++) {
                                    if (ai == 0 && !e->splice_sites->canonical_acceptors[0].gff_canonical) continue; // Canonical acceptor is lost, can't be used
                                    Range ncr = { donor_pos, target->splice_sites->canonical_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }

                                for (int ai = 0; ai < target->splice_sites->n_cryptic_acceptors; ai++) {
                                    Range ncr = { donor_pos, target->splice_sites->cryptic_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            for (int di = 0; di < e->splice_sites->n_cryptic_donors; di++) {
                                uint32_t donor_pos = e->splice_sites->cryptic_donors[di].position;

                                if (target == NULL) {
                                    Range ncr = { donor_pos, donor_pos }; // No splicing happens for this donor
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                    continue;
                                }

                                for (int ai = 0; ai < target->splice_sites->n_canonical_acceptors; ai++) {
                                    if (ai == 0 && !e->splice_sites->canonical_acceptors[0].gff_canonical) continue; // Canonical acceptor is lost, can't be used
                                    Range ncr = { donor_pos, target->splice_sites->canonical_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }

                                for (int ai = 0; ai < target->splice_sites->n_cryptic_acceptors; ai++) {
                                    Range ncr = { donor_pos, target->splice_sites->cryptic_acceptors[ai].position };
                                    kv_push(Range, nc_regions[n_nc_regions], ncr);
                                }
                            }

                            n_nc_regions++;
                        }
                    }

                    if (!is_coding) continue; // Transcript is not protein-coding after all!

                    int nisoforms = 0, misoforms = 1;
                    for (int i = 0; i < n_nc_regions; i++) {
                        Ranges ranges = nc_regions[i];
                        if (ranges.n == 0) continue;

                        misoforms *= ranges.n;
                    }

                    Ranges **isoforms = malloc(misoforms * sizeof(Ranges *));
                    for (int i = 0; i < misoforms; i++) {
                        isoforms[i] = malloc(sizeof(Ranges));
                        kv_init(*isoforms[i]);
                    }

                    for (int i = 0; i < n_nc_regions; i++) {
                        Ranges ranges = nc_regions[i];
                        if (ranges.n == 0) continue;

                        if (nisoforms == 0) { // First isoform
                            for (int j = 0; j < ranges.n; j++) {
                                Ranges *first = isoforms[j];
                                Range to_push = kv_A(ranges, j);
                                kv_push(Range, *first, to_push);
                                nisoforms++;
                            }
                        } else { // Extend/duplicate isoforms already created
                            // Copy over isoforms to add new endings
                            int prior_n_isoforms = nisoforms;
                            for (int j = 1; j < ranges.n; j++) {
                                for (int k = 0; k < prior_n_isoforms; k++) {
                                    Ranges *original = isoforms[k];
                                    Ranges *new = isoforms[j * prior_n_isoforms + k];
                                    kv_copy(Range, *new, *original);
                                    nisoforms++;
                                }
                            }

                            for (int j = 0; j < ranges.n; j++) {
                                for (int k = 0; k < prior_n_isoforms; k++) {
                                    Ranges *target = isoforms[j * prior_n_isoforms + k];
                                    Range to_push = kv_A(ranges, j);
                                    kv_push(Range, *target, to_push);
                                }
                            }
                        }
                    }

                    kstring_t hap_seq;
                    get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_seq);

                    log_info("Number of isoforms: %i", nisoforms);
                    for (int i = 0; i < nisoforms; i++) {
                        Ranges *ncrs = isoforms[i];
                        kstring_t spliced_seq = {0};
                        int seq_idx = 0;
                        log_info("Processing isoform: %i", i);
                        log_info("Gene length is: %i", hap_seq.l);

                        int gene_start = gene->coords.start;
                        if (t->strand == STRAND_FWD) {
                            int length = 0;
                            for (int j = 0; j < ncrs->n; j++) {
                                Range ncr = kv_A(*ncrs, j);
                                int start = ncr.start - gene_start;
                                int end = ncr.end - gene_start;
                                if (start < seq_idx) {
                                    continue;
                                }

                                log_info("Adding region %i, %i:%i to isoform %i", j, seq_idx, start, i);
                                length += start - seq_idx;

                                kputsn(hap_seq.s + seq_idx, start - seq_idx, &spliced_seq);
                                seq_idx = end;
                            }
                            log_info("Lengths: %i - %i", spliced_seq.l, length);
                            kputsn(hap_seq.s + seq_idx, hap_seq.l - seq_idx, &spliced_seq);
                        } else if (t->strand == STRAND_REV) {
                            int length = 0;
                            for (int j = ncrs->n-1; j >= 0; j--) {
                                Range ncr = kv_A(*ncrs, j);
                                int start = ncr.end - gene_start;
                                int end = ncr.start - gene_start;
                                if (start < seq_idx) {
                                    continue;
                                }

                                log_info("Adding region %i, %i:%i to isoform %i", j, seq_idx, start, i);
                                length += start - seq_idx;
                                kputsn(hap_seq.s + seq_idx, start - seq_idx, &spliced_seq);
                                seq_idx = end;
                            }
                            log_info("Lengths: %i - %i", spliced_seq.l, length);
                            kputsn(hap_seq.s + seq_idx, hap_seq.l - seq_idx, &spliced_seq);
                        }

                        for (int ri = 0; ri < batch->n; ri++) {
                            
                        }

                        // Frameshift check
                        // Stop gain check



                        // if (spliced_seq.s != NULL) free(spliced_seq.s);
                    }


                    // Loop over isoforms and splice
                    // Generate sequence
                    // Generate complementary sequence if necessary
                    // Walk over sequence while keeping track of records that are responsible for changes
                    // Look for frameshift rescues
                    // Consider stop gains

                    kv_destroy(hap_in);
                    free(hap_seq.s);
                }
            }

            kv_destroy(*batch);
            free(batch);
        }
        kv_destroy(*phased_records);
        free(phased_records);

        for (int i = 0; i < block.records->n; i++) {
            Record *r = kv_A(*block.records, i);
            if (!r->triggers[SPLICE_TRIGGER]) continue;

            kstring_t splice_scores = {0};
            kstring_t splice_change = {0};

            int has_scores = 0;
            for (int g = 0; g < r->genes->n; g++) {
                if (g > 0) {
                    kputc(',', &splice_scores);
                    kputc(',', &splice_change);
                }

                for (int p = 0; p < NUM_PHASES; p++) {

                    if (p > 0) {
                        kputc('&', &splice_scores);
                        kputc('&', &splice_change);
                    }

                    Score *score = r->complete_scores[g * NUM_PHASES + p];
                    kputs(r->v->d.allele[1], &splice_scores);
                    if (score == NULL) {
                        kputc('|', &splice_scores); kputs(kv_A(*r->genes, g)->name, &splice_scores);
                        kputs("|.|.|.|.|.|.|.|.", &splice_scores);
                        kputs(".|.|.|.|.", &splice_change);
                        continue;
                    }
                    has_scores = 1;

                    kputc('|', &splice_scores); kputs(score->gene, &splice_scores);

                    // Format and append DS scores
                    char ds_scores[20];
                    snprintf(ds_scores, 20, "%.2f|%.2f|%.2f|%.2f", score->ag, score->al, score->dg, score->dl);
                    kputc('|', &splice_scores); kputsn(ds_scores, 19, &splice_scores);

                    kputc('|', &splice_scores); kputw(score->ag_idx, &splice_scores);
                    kputc('|', &splice_scores); kputw(score->al_idx, &splice_scores);
                    kputc('|', &splice_scores); kputw(score->dg_idx, &splice_scores);
                    kputc('|', &splice_scores); kputw(score->dl_idx, &splice_scores);

                    // Splice change
                    SpliceChange *change = r->cs[g * NUM_PHASES + p];

                    char change_scores[29];
                    snprintf(change_scores, 29, "%.2f|%.2f|%.2f|%.2f|%.2f", change->AG, change->AL, change->DG, change->DL, change->sum);
                    kputsn(change_scores, strlen(change_scores), &splice_change);
                }
            }

            if (has_scores) log_info(splice_scores.s);

            bcf_update_info_string(hdr, r->v, SPLICEAI_HAP_TAG, splice_scores.s);
            bcf_update_info_string(hdr, r->v, CS_TAG, splice_change.s);

            free(splice_scores.s);
            free(splice_change.s);
        }
    }

    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);

        for (int i = 0; i < block.records->n; i++) {
            Record *r = kv_A(*block.records, i);
            bcf_write(fp_out, hdr, r->v);
        }
    }


    // Free genes hashmap
    for (int k = kh_begin(genes); k != kh_end(genes); ++k) {
        if (kh_exist(genes, k)) {
            Gene *gene = kh_value(genes, k);
            free(gene->seq.s);
            free(gene);
        }
    }
    kh_destroy_GeneHash(genes);

    // Free CDS ranges
    // for (int k = kh_begin(cds_ranges); k != kh_end(cds_ranges); ++k) {
    //     if (!kh_exist(cds_ranges, k)) break;
    //
    //     Range *cd = kh_value(cds_ranges, k);
    //     free(cd);
    // }

    // Free batches
    for (int i = 0; i < batches.n; i++) {
        RecordBlock block = kv_A(batches, i);

        for (int j = 0; j < block.records->n; j++) {
            Record *r = kv_A(*block.records, j);
            destroy_record(r);
        }
        kv_destroy(*block.records);
        free(block.records);
    }
    kv_destroy(batches);

    // Free raw scores file
    // fclose(raw_scores);

    // Free remaining allocs
    gff_destroy(gff);
    bcf_hdr_destroy(hdr);
    bcf_close(fp_out);


    // Destroy SpliceAI models
    destroy_models(models);
}

