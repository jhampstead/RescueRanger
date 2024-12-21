#ifndef HAPLOTYPE_HEADER
#define HAPLOTYPE_HEADER

#include "klib/kstring.h"
#include "cpliceai/range.h"
#include "record.h"
#include "utils.h"

void get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap);



#endif
