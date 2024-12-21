#ifndef MOCKS_GUARD
#define MOCKS_GUARD

#include <stdint.h>
#include <htslib/vcf.h>

bcf1_t * mock_bcf1_t(const uint64_t position, char *ref, char *alt);

#endif
