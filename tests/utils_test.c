#include <htslib/vcf.h>
#include <unity.h>
#include <utils.h>

#include "logging/log.h"

bcf_hdr_t *init_hdr(void) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");
    bcf_hdr_append(hdr, "##fileformat=VCFv4.2\n");
    bcf_hdr_append(hdr, "##contig=<ID=1,length=249250621>\n");
    bcf_hdr_add_sample(hdr, "sample1");
    bcf_hdr_sync(hdr);

    return hdr;
}

bcf1_t *init_variant(int rid, int pos, char *ref, char *alt) {
    bcf_hdr_t *hdr = init_hdr();
    bcf1_t *v = bcf_init();

    v->rid = rid;
    v->pos = pos;
    v->n_allele = 2;
    v->d.allele = malloc(v->n_allele * sizeof(char *));
    v->rlen = strlen(ref);
    v->d.allele[0] = ref;
    v->d.allele[1] = alt;

    bcf_hdr_destroy(hdr);

    return v;
}

Variants as_kvec(bcf1_t ** vs, size_t n) {
    return (Variants) {n, n, vs};
}

void setUp(void) {}
void tearDown(void) {}

void GetHaplotypeSequence_SingleSNV_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 2, "C", "T");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; get_haplotype_sequence(ref, vs, 0, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATTGT", hap.s);

    bcf_destroy(v);
    free(hap.s);
}

void GetHaplotypeSequence_MultipleSNVs_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v[3];
    v[0] = init_variant(0, 2, "C", "A");
    v[1] = init_variant(0, 3, "G", "A");
    v[2] = init_variant(0, 4, "T", "A");
    Variants vs = as_kvec(v, 3);

    kstring_t hap; get_haplotype_sequence(ref, vs, 0, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATAAA", hap.s);

    bcf_destroy(v[0]); bcf_destroy(v[1]); bcf_destroy(v[2]);
    free(hap.s);
}

void GetHaplotypeSequence_SingleInsertion_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 2, "C", "CC");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; get_haplotype_sequence(ref, vs, 1, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATCCGT", hap.s);

    bcf_destroy(v);
    free(hap.s);
}

void GetHaplotypeSequence_SingleDeletion_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 2, "CG", "C");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; get_haplotype_sequence(ref, vs, -1, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATCT", hap.s);

    bcf_destroy(v);
    free(hap.s);
}

void GetHaplotypeSequence_PositionLessThanZero_ShouldFail(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, -1, "C", "G");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 0, 0, &hap);

    TEST_ASSERT_EQUAL(1, res);

    bcf_destroy(v);
}

void GetHaplotypeSequence_PositionMoreThanLength_ShouldFail(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 7, "C", "G");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 0, 0, &hap);

    TEST_ASSERT_EQUAL(1, res);

    bcf_destroy(v);
}

void GetHaplotypeSequence_PositionOffset_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 7, "C", "G");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 0, 5, &hap);

    TEST_ASSERT_EQUAL(0, res);
    TEST_ASSERT_EQUAL_STRING("ATGGT", hap.s);

    bcf_destroy(v);
    free(hap.s);
}


void GetHaplotypeSequence_TouchingInsertions_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCAT" };
    bcf1_t *v[2];
    v[0] = init_variant(0, 2, "C", "CACGT");
    v[1] = init_variant(0, 3, "A", "ACGTT");
    Variants vs = as_kvec(v, 2);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 8, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATCACGTACGTTT", hap.s);

    bcf_destroy(v[0]); bcf_destroy(v[1]);
    free(hap.s);
}

void GetHaplotypeSequence_TouchingInsertionSNV_ShouldSucceed(void) {
    kstring_t ref = {5, 6, "ATCAT" };
    bcf1_t *v[2];
    v[0] = init_variant(0, 2, "C", "CC");
    v[1] = init_variant(0, 3, "A", "T");
    Variants vs = as_kvec(v, 2);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 1, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATCCTT", hap.s);
    bcf_destroy(v[0]); bcf_destroy(v[1]);
    free(hap.s);
}

void GetHaplotypeSequence_TouchingDeletionSNV_ShouldSucceed(void) {
    kstring_t ref = {5, 6, "ATCAT" };
    bcf1_t *v[2];
    v[0] = init_variant(0, 2, "CA", "C");
    v[1] = init_variant(0, 4, "T", "G");
    Variants vs = as_kvec(v, 2);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, -1, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATCG", hap.s);

    bcf_destroy(v[0]); bcf_destroy(v[1]);
    free(hap.s);
}

// void GetHaplotypeSequence_MismatchedReference_ShouldFail(void) {
//     kstring_t ref = {5, 6, "ATCGT" };
//     bcf1_t *v = init_variant(0, 3, "A", "T");
//     Variants vs = as_kvec(&v, 1);
//
//     kstring_t hap; int res = get_haplotype_sequence(ref, vs, 0, 0, &hap);
//
//     TEST_ASSERT_EQUAL(1, res);
//
//     bcf_destroy(v);
// }

void GetHaplotypeSequence_HasOverlap_ShouldSuccedWithWarning(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v[2];
    v[0] = init_variant(0, 2, "CGT", "C");
    v[1] = init_variant(0, 3, "G", "T");
    Variants vs = as_kvec(v, 2);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, -2, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATC", hap.s);

    bcf_destroy(v[0]); bcf_destroy(v[1]);
    free(hap.s);
}

void GetHaplotypeSequence_BeyondStart_ShouldSucceed(void) {
    kstring_t ref = { 5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 0, "AAAT", "A");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, -2, 2, &hap);

    TEST_ASSERT_EQUAL_STRING("CGT", hap.s);

    bcf_destroy(v);
    free(hap.s);
}

void GetHaplotypeSequence_BeyondEnd_ShouldSucceed(void) {
    kstring_t ref = {5, 6, "ATCGT" };
    bcf1_t *v = init_variant(0, 0, "CGTAAAA", "C");
    Variants vs = as_kvec(&v, 1);

    kstring_t hap; int res = get_haplotype_sequence(ref, vs, 6, 0, &hap);

    TEST_ASSERT_EQUAL_STRING("ATC", hap.s);

    bcf_destroy(v);
    free(hap.s);
}


int main(void) {
    log_set_level(2);
    UNITY_BEGIN();

    RUN_TEST(GetHaplotypeSequence_SingleSNV_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_MultipleSNVs_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_SingleInsertion_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_SingleDeletion_ShouldSucceed);
    // RUN_TEST(GetHaplotypeSequence_PositionLessThanZero_ShouldFail);
    // RUN_TEST(GetHaplotypeSequence_PositionMoreThanLength_ShouldFail);
    RUN_TEST(GetHaplotypeSequence_PositionOffset_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_TouchingInsertions_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_TouchingInsertionSNV_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_TouchingDeletionSNV_ShouldSucceed);
    // RUN_TEST(GetHaplotypeSequence_MismatchedReference_ShouldFail);
    RUN_TEST(GetHaplotypeSequence_HasOverlap_ShouldSuccedWithWarning);
    RUN_TEST(GetHaplotypeSequence_BeyondStart_ShouldSucceed);
    RUN_TEST(GetHaplotypeSequence_BeyondEnd_ShouldSucceed);

    return UNITY_END();
}
