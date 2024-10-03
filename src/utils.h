#include <stdint.h>
#include <htslib/vcf.h>

#include <klib/kstring.h>
#include <klib/kvec.h>

#include "record.h"

int get_allele_offset(int64_t pos, int rlen, int alen, int boundary_start, int boundary_end);

int get_allele_offset_variant(bcf1_t *v, int boundary_start, int boundary_end);

int64_t get_position_in_haplotype(int64_t pos, Records *rs, int64_t start, int64_t end);

int get_haplotype_sequence(kstring_t ref, Records *rs, int seq_offset, int pos_offset, kstring_t *out_hap);


