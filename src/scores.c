#include "scores.h"
#include "klib/kvec.h"

void parse_prediction_pairs(const float *ref_priors, const float *ref_posteriors, const float *alt_priors, const float *alt_posteriors, const Records *rs, const Range pred_range, PairPrediction **ref_alt_maps, size_t *ref_alt_map_size) {
    kvec_t(PairPrediction) map;
    kv_init(map);

    uint64_t curr_position = pred_range.start;
    int refs_idx = 0, alts_idx = 0;
    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        Range ref = variant_to_range(r->v, 0);
        Range intersected_ref;
        if (intersect_ranges(ref, pred_range, &intersected_ref) == 0) continue;
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt;
        Range intersected_alt;
        int alt_segment = 0;
        if (intersect_ranges(alt, pred_range, &intersected_alt) > 0) {
            const int skipped_bases = intersected_alt.start - alt.start;
            alt_segment = range_len(alt) - skipped_bases;
        }

        while (curr_position < r->v->pos) {
            Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        int offset = allele_offset(r->v, 0, 1, pred_range);
        if (offset == 0) continue; // Doesn't need special handling, can be copied over in while loops

        // Copy over preceding reference base, if it hasn't already (e.g. no overlapping reference base should be copied over in the case of touching indels)
        if (curr_position == r->v->pos) {
            Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }


        Prediction empty_pred = { true, 0, 0.0 };
        if (offset > 0) { // Deletion
            Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            PairPrediction pp = { ref_pred, empty_pred };
            kv_push(PairPrediction, map, pp);
            curr_position++;
        } else { // Offset < 0, insertion
            Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            PairPrediction pp = { empty_pred, alt_pred };
            kv_push(PairPrediction, map, pp);
        }

    }

    while (curr_position < pred_range.end) {
        Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
        refs_idx++;
        Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
        alts_idx++;
        PairPrediction pp = {ref_pred, alt_pred};
        kv_push(PairPrediction, map, pp);
        curr_position++;
    }

    *ref_alt_maps = map.a;
    *ref_alt_map_size = map.n;
}

