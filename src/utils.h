#include <stdint.h>
#include <htslib/vcf.h>

#include <klib/kstring.h>
#include <klib/kvec.h>

#include "record.h"


Range variant_to_range(const bcf1_t *v, const int allele_index);

int allele_offset(const bcf1_t *v, const int ref_index, const int alt_index, const Range bounds);

int haplotype_offset(const Records *rs, const Range bounds);

int prediction_offset(const bcf1_t *v, const Records *rs);

float *match_predictions_to_reference(const float *predictions, const Range pred_range, const Records *rs, const size_t num_reference_predictions);

