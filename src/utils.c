#include "utils.h"
#include "region.h"
#include <cpliceai/utils.h>


Range variant_to_range(const bcf1_t *v, const int allele) {
    return (Range) { v->pos, v->pos + strlen(v->d.allele[allele]) };
}

int allele_offset(const bcf1_t *v, const int ref_index, const int alt_index, const Range bounds) {
    Range ref = variant_to_range(v, ref_index);
    Range ref_range;
    if (intersect_ranges(ref, bounds, &ref_range) == 0 || range_len(ref_range) == 0) {
        return 0;
    }

    Range alt = variant_to_range(v, alt_index);
    Range intersected_alt;
    int alt_len = 0;
    if (intersect_ranges((Range) {v->pos, v->pos+1}, bounds, &intersected_alt) != 0 && range_len(intersected_alt) > 0) {
        alt_len = range_len(alt);
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

int prediction_offset(const bcf1_t *v, const Records *rs) {
    return haplotype_offset(rs, (Range) { v->pos - 50, v->pos + 51 }); // TODO: Change to a const/defined variable
}

void get_best_insertion_score(const float *insert_start, const int length, float *best) {
    for (int i = 0; i < length; i+=3) {
        for (int j = 0; j < 3; j++) {
            if (insert_start[i + j] > best[j]) best[j] = insert_start[i+j];
        }
    }
}

float *match_predictions_to_reference(const float *predictions, const bcf1_t *v, const Records *rs, const size_t num_reference_predictions) {
    float *corrected_predictions = malloc(num_reference_predictions * sizeof(float));
    Range pred_range = {v->pos - (50), v->pos + 51}; // +1 for open ended coordinates

    int ref_idx = 0, hap_idx = 0;
    for (int l = 0; l < rs->n; l++) {
        Record *r = kv_A(*rs, l);

        Range ref;
        if (intersect_ranges(variant_to_range(r->v, 0), pred_range, &ref) == 0) continue;
        shift_range(&ref, -pred_range.start);
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt = variant_to_range(r->v, 1);
        Range intersected_alt;
        int alt_segment = 0;
        if (intersect_ranges(alt, pred_range, &intersected_alt) > 0) {
            // shift_range(&alt, -pred_range.start);
            alt_segment = range_len(alt);
        }

        int skipped_alt_bases = (intersected_alt.start - alt.start);
        if (alt_segment > 0) alt_segment -= skipped_alt_bases;
        Range ref_segment = scale_range(ref, NUM_SCORES);

        // Copy segments between variants/boundaries
        int segment_len = (ref_segment.start - ref_idx);
        if (segment_len > 0) memcpy(corrected_predictions + ref_idx, predictions + hap_idx, segment_len * sizeof(float));

        int variant_offset = haplotype_offset(rs, (Range) { pred_range.start, pred_range.start + ref.start }) * NUM_SCORES;
        int offset = ref_len - alt_segment; // Naive offset in end position vs ref and alt
        float *corrected_start = corrected_predictions + ref_segment.start;
        if (offset == 0) {
            memcpy(corrected_start, predictions + ref_segment.start - variant_offset, range_len(ref_segment) * sizeof(float));
        } else if (offset < 0) {
            float best_scores[NUM_SCORES] = {0};
            get_best_insertion_score(predictions + ref_segment.start - variant_offset, alt_segment * NUM_SCORES, best_scores);
            memcpy(corrected_start, best_scores, NUM_SCORES * sizeof(float));
        } else { // offset > 0
            memcpy(corrected_start, predictions + ref_segment.start - variant_offset, NUM_SCORES * sizeof(float));
            memset(corrected_start + NUM_SCORES, 0, (offset - (skipped_alt_bases > 0)) * NUM_SCORES * sizeof(float));
        }


        if (ref_segment.end > ref_idx) ref_idx = ref_segment.end;
    }

    int variant_offset = haplotype_offset(rs, pred_range) * NUM_SCORES;
    memcpy(corrected_predictions + ref_idx, predictions + ref_idx - variant_offset, (num_reference_predictions - ref_idx) * sizeof(float));

    return corrected_predictions;
}
