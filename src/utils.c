#include "utils.h"
#include "region.h"
#include <cpliceai/utils.h>


Range variant_to_range(const bcf1_t *v, const int allele) {
    return (Range) { v->pos, v->pos + strlen(v->d.allele[allele]) };
}

bool variant_fits_window(const Alleles *a, const int prediction_radius) {
    return a->alt_len <= prediction_radius && a->ref_len <= prediction_radius;
}

// The base the ALT of a variant is attached to: the base before an indel, the first base of a substitution
static int64_t attachment_base(const Alleles *a) {
    return a->pos - (a->prefix_len > 0);
}

int64_t variant_position(const Alleles *a) {
    return a->ref_len > 0 ? a->pos : attachment_base(a);
}

Range variant_footprint(const Alleles *a) {
    return (Range) { attachment_base(a), a->pos + a->ref_len };
}

bool variant_overlaps(const Alleles *a, const Range range) {
    Range intersection;
    return intersect_ranges(variant_footprint(a), range, &intersection) != 0 && range_len(intersection) > 0;
}

bool variant_in_region(const Alleles *a, const Range region) {
    // An insertion changes no REF base. It lies between the bases pos - 1 and pos, and is in the region when both are
    if (a->ref_len == 0) return a->pos - 1 >= region.start && a->pos < region.end;

    return a->pos < region.end && a->pos + a->ref_len > region.start;
}

int64_t variant_shift(const Alleles *a, const int64_t b, const bool insertions_at_b) {
    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
    int64_t shift = 0;

    // The REF bases left over after pairing them with the ALT bases are deleted: the ones before b move what comes after them
    const int64_t deleted_start = a->pos + paired;
    const int64_t deleted_end = a->pos + a->ref_len < b ? a->pos + a->ref_len : b;
    if (deleted_end > deleted_start) shift += deleted_end - deleted_start;

    // The ALT bases left over are inserted at the point after the last REF base, between the bases point - 1 and point
    if (a->alt_len > a->ref_len) {
        const int64_t point = a->pos + a->ref_len;
        if (point < b || (point == b && insertions_at_b)) shift -= a->alt_len - a->ref_len;
    }

    return shift;
}

int allele_offset(const Alleles *a, const Range bounds) {
    if (bounds.end <= bounds.start) return 0;

    // The deleted bases in bounds, and the inserted bases whose attached base (the base before their point) is in bounds, which
    // is when the point is in (start, end]
    return (int) (variant_shift(a, bounds.end, true) - variant_shift(a, bounds.start, true));
}

int haplotype_offset(const Records *rs, const Range bounds) {
    int offset = 0;

    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);
        offset += allele_offset(&r->alleles, bounds);
    }

    return offset;
}

int64_t haplotype_shift(const Records *rs, const Range gene, const int64_t b, const bool insertions_at_b) {
    if (b <= gene.start) return 0;

    int64_t shift = 0;

    for (size_t i = 0; i < kv_size(*rs); i++) {
        const Alleles *a = &kv_A(*rs, i)->alleles;
        if (!variant_in_region(a, gene)) continue;

        // Only what lies in the gene counts: the shift at b minus the shift at the start of the gene
        shift += variant_shift(a, b, insertions_at_b) - variant_shift(a, gene.start, insertions_at_b);
    }

    return shift;
}

Range region_to_haplotype(const Records *rs, const Range gene, const Range region, const bool is_utr) {
    const int64_t start = region.start - haplotype_shift(rs, gene, region.start, !is_utr);

    // An empty region has no inside: an insertion at its position cannot be before its start and after its end, which would
    // turn it around. It stays empty.
    if (region.end <= region.start) return (Range) { start, start };

    return (Range) { start, region.end - haplotype_shift(rs, gene, region.end, is_utr) };
}

Range reference_range_to_haplotype(const Records *rs, const Range gene, const Range range, const uint32_t strand) {
    // An insertion at b lies between the bases b - 1 and b. Forward it is attached to b - 1, so it counts as before b at both
    // borders. Reverse it is attached to b, and does not
    const bool insertions_at_b = strand != STRAND_REV;

    const int64_t start = range.start - haplotype_shift(rs, gene, range.start, insertions_at_b);
    if (range.end <= range.start) return (Range) { start, start };

    return (Range) { start, range.end - haplotype_shift(rs, gene, range.end, insertions_at_b) };
}

int prediction_offset(const Alleles *a, const Records *rs, const int prediction_radius) {
    return haplotype_offset(rs, range_around(variant_footprint(a).start, prediction_radius));
}

void get_best_insertion_score(const float *insert_start, const int length, float *best) {
    for (int i = 0; i < length; i+=3) {
        for (int j = 0; j < 3; j++) {
            if (insert_start[i + j] > best[j]) best[j] = insert_start[i+j];
        }
    }
}

float *match_predictions_to_reference(const float *predictions, const Range pred_range, const Records *rs, const size_t num_reference_predictions) {
    float *corrected = malloc(num_reference_predictions * sizeof(float));
    const size_t score_size = NUM_SCORES * sizeof(float);

    int64_t idx = 0; // First position of the window (relative to its start) that has not been written yet
    for (int l = 0; l < rs->n; l++) {
        const Alleles *a = &kv_A(*rs, l)->alleles;
        if (!variant_overlaps(a, pred_range)) continue;

        // The positions of the window the variant touches, relative to the window start
        Range touched;
        intersect_ranges(variant_footprint(a), pred_range, &touched);
        shift_range(&touched, -pred_range.start);

        // Up to the variant, the haplotype predictions are shifted by the variants before it
        const int64_t shift = haplotype_offset(rs, (Range) { pred_range.start, pred_range.start + touched.start });
        const float *hap = predictions - shift * NUM_SCORES; // hap + q * NUM_SCORES is the prediction for window position q

        if (touched.start > idx) memcpy(corrected + idx * NUM_SCORES, hap + idx * NUM_SCORES, (touched.start - idx) * score_size);

        if (a->ref_len == a->alt_len) { // A substitution replaces its bases one by one
            for (int64_t q = touched.start; q < touched.end; q++) memcpy(corrected + q * NUM_SCORES, hap + q * NUM_SCORES, score_size);
        } else {
            // The REF and ALT bases are paired one to one from the left. The junction is the base the indel is attached to, if
            // that base lies in the window; it keeps its own prediction. The paired bases keep the prediction of the base that
            // takes their place, and the REF bases left over are removed and get 0, also the first one of a deletion whose
            // attachment base lies before the window.
            const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
            const bool junction_in_window = a->prefix_len > 0 && variant_footprint(a).start >= pred_range.start;
            for (int64_t q = touched.start; q < touched.end; q++) {
                const bool kept = (junction_in_window && q == touched.start) || pred_range.start + q < a->pos + paired;
                if (kept) memcpy(corrected + q * NUM_SCORES, hap + q * NUM_SCORES, score_size);
                else memset(corrected + q * NUM_SCORES, 0, score_size);
            }

            // The ALT bases left over are inserted after the last paired base (the base the insertion is attached to when
            // nothing is paired). If it is in the window it takes the best of itself and the inserted bases.
            if (a->alt_len > a->ref_len) {
                const int64_t last_paired = a->pos + paired - 1;
                if (last_paired >= pred_range.start && last_paired < pred_range.end) {
                    const int64_t q = last_paired - pred_range.start;
                    float best_scores[NUM_SCORES] = {0};
                    get_best_insertion_score(hap + q * NUM_SCORES, (1 + a->alt_len - a->ref_len) * NUM_SCORES, best_scores);
                    memcpy(corrected + q * NUM_SCORES, best_scores, score_size);
                }
            }
        }

        if (touched.end > idx) idx = touched.end;
    }

    const int64_t shift = haplotype_offset(rs, pred_range);
    memcpy(corrected + idx * NUM_SCORES, predictions + (idx - shift) * NUM_SCORES, (num_reference_predictions - idx * NUM_SCORES) * sizeof(float));

    return corrected;
}
