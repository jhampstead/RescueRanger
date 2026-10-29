#include "scores.h"
#include "klib/kvec.h"

void parse_prediction_pairs(const float *ref_priors, const float *ref_posteriors, const float *alt_priors, const float *alt_posteriors, const Records *rs, const Range pred_range, PairPrediction **ref_alt_maps, size_t *ref_alt_map_size) {
    kvec_t(PairPrediction) map;
    kv_init(map);

    uint64_t curr_position = pred_range.start;
    int refs_idx = 0, alts_idx = 0;
    for (int i = 0; i < rs->n; i++) {
        const Record *r = kv_A(*rs, i);

        // The variant starts at the base it is attached to (the base before an indel), which is where an indel's pairs start
        const Range footprint = variant_footprint(&r->alleles);
        Range touched;
        if (intersect_ranges(footprint, pred_range, &touched) == 0) continue;
        const int64_t attached_to = footprint.start;

        while (curr_position < attached_to) {
            const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            const PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        const int offset = allele_offset(&r->alleles, pred_range);
        if (offset == 0) continue; // Doesn't need special handling, can be copied over in while loops

        // The base the variant is attached to keeps an ordinary pair, unless it was already paired (touching indels)
        if (curr_position == attached_to) {
            const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            const PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        // The REF and ALT bases are paired one to one from the left: the paired bases of a complex indel are ordinary pairs, and
        // only the bases left over after them are deleted or inserted. For a pure insertion or deletion nothing is paired.
        const int paired = r->alleles.ref_len < r->alleles.alt_len ? r->alleles.ref_len : r->alleles.alt_len;
        while ((int64_t) curr_position < r->alleles.pos + paired && (int64_t) curr_position < pred_range.end) {
            const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            const PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        const Prediction empty_pred = { true, 0, 0.0, 0.0 };
        if (offset > 0) { // Deletion
            for (int i = 0; i < offset; i++) {
                const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
                refs_idx++;
                const PairPrediction pp = { ref_pred, empty_pred };
                kv_push(PairPrediction, map, pp);
                curr_position++;
            }
        } else { // Offset < 0, insertion
            for (int i = 0; i < -offset; i++) {
                const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
                alts_idx++;
                const PairPrediction pp = { empty_pred, alt_pred };
                kv_push(PairPrediction, map, pp);
            }
        }

    }

    while (curr_position < pred_range.end) {
        const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
        refs_idx++;
        const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
        alts_idx++;
        const PairPrediction pp = {ref_pred, alt_pred};
        kv_push(PairPrediction, map, pp);
        curr_position++;
    }

    *ref_alt_maps = map.a;
    *ref_alt_map_size = map.n;
}

