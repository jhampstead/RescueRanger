#ifndef SCORE_H
#define SCORE_H

#include <cpliceai/range.h>
#include <stdbool.h>
#include <stdint.h>

#include "record.h"
#include "utils.h"


// A change of a splice site probability from this size on counts as a change: a site that weakens by this much is ablated, a
// score that differs this much between the haplotype and the variant on its own is haplotype specific, and a position with a
// probability above it is written to the raw scores
#define SPLICEAI_RESCUE_THRESHOLD 0.2

typedef struct {
    bool missing;
    uint64_t position;
    float prior;
    float posterior;
} Prediction;

typedef struct {
    Prediction ref;
    Prediction alt;
} PairPrediction;


void parse_prediction_pairs(const float *ref_priors, const float *ref_posteriors, const float *alt_priors, const float *alt_posteriors, const Records *rs, const Range pred_range, PairPrediction **ref_alt_maps, size_t *ref_alt_map_size);


#endif
