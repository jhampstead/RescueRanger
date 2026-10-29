#ifndef MOCKS_GUARD
#define MOCKS_GUARD

#include <check.h>
#include <stdbool.h>
#include <stdint.h>
#include <htslib/vcf.h>

// A bare record (positions and alleles only) for ref and alt. The alleles are stored, not copied, so they have to outlive the
// record. The record belongs to the mocks: mocks_free_all frees it, so do not bcf_destroy it. Create the test cases of a suite
// that uses it with mocked_tcase, which calls mocks_free_all after every test.
bcf1_t * mock_bcf1_t(const uint64_t position, const char *ref, const char *alt);

// Frees every record mock_bcf1_t made since the last call. A copy made from one (Record_init does bcf_dup) is not affected.
void mocks_free_all(void);

// The number of records mock_bcf1_t made since the last call of mocks_free_all
size_t mocks_live_count(void);

// tcase_create, plus mocks_free_all as the checked teardown, so it runs in the process of the test
TCase *mocked_tcase(const char *name);

// A VCF header with the contig chr1, INFO/AF, FORMAT/GT and a single sample. Free it with bcf_hdr_destroy.
bcf_hdr_t *mock_hdr(void);

// A real record on chr1 for the header of mock_hdr. alleles is "REF,ALT" (for example "A,G"), the genotype is gt0 and gt1
// (allele indices), phased or not, and af is INFO/AF, or a negative value for a record without INFO/AF. Free it with bcf_destroy.
bcf1_t *mock_genotyped_bcf1_t(bcf_hdr_t *hdr, const int64_t position, const char *alleles, const int gt0, const int gt1, const bool phased, const float af);

#endif
