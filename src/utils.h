#ifndef RESCUERANGER_UTILS_H
#define RESCUERANGER_UTILS_H

#include <limits.h>
#include <stdint.h>
#include <htslib/vcf.h>

#include <klib/kstring.h>
#include <klib/kvec.h>

#include "record.h"

// Bases of sequence context the model needs on each side of the prediction window.
#ifndef SPLICEAI_CONTEXT
#define SPLICEAI_CONTEXT 5000
#endif

// Largest prediction radius for which pred_array_len (the first derived size to overflow) fits in an int.
#define MAX_PREDICTION_RADIUS (((INT_MAX / NUM_SCORES) - 1) / 2)

// Half-open range [pos - offset, pos + offset + 1). Not clamped, so it can start below 0.
static inline Range range_around(const int64_t pos, const int offset) {
    return (Range) { pos - offset, pos + offset + 1 };
}

// Half-open range [pos - offset, pos), the left half of range_around without pos itself.
static inline Range range_before(const int64_t pos, const int offset) {
    return (Range) { pos - offset, pos };
}

// Number of positions in a prediction window of radius `prediction_radius`.
static inline int pred_width(const int prediction_radius) {
    return 2 * prediction_radius + 1;
}

// Length of the model input sequence for a prediction window of radius `prediction_radius`.
static inline int model_width(const int prediction_radius) {
    return pred_width(prediction_radius) + 2 * SPLICEAI_CONTEXT;
}

// Position within the model input sequence that corresponds to the variant's genomic position.
static inline int variant_position_in_sequence(const int prediction_radius) {
    return SPLICEAI_CONTEXT + prediction_radius;
}

// Number of floats in a prediction array (NUM_SCORES per position).
static inline int pred_array_len(const int prediction_radius) {
    return pred_width(prediction_radius) * NUM_SCORES;
}

Range variant_to_range(const bcf1_t *v, const int allele_index);

// True when neither the alt allele nor the reference span is longer than the prediction radius.
bool variant_fits_window(const bcf1_t *v, const int prediction_radius);

int allele_offset(const bcf1_t *v, const int ref_index, const int alt_index, const Range bounds);

int haplotype_offset(const Records *rs, const Range bounds);

int prediction_offset(const bcf1_t *v, const Records *rs, const int prediction_radius);

float *match_predictions_to_reference(const float *predictions, const Range pred_range, const Records *rs, const size_t num_reference_predictions);

#endif
