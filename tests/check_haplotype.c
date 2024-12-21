#include "check_haplotype.h"

#include <check.h>
#include <klib/kstring.h>

#include "haplotype.h"
#include "klib/kvec.h"
#include "mocks.h"
#include "record.h"

kstring_t seq;
kstring_t hap;
Records rs;
static Range range = {0, 20};

void checked_setup() {
    kv_init(rs);

    seq = (kstring_t) {0};
    kputs("AAAAAAAAAAAAAAAAAAAA", &seq);

    hap = (kstring_t) {0};
}

void checked_teardown() {
    kv_destroy(rs);

    free(seq.s);
    if (hap.s != NULL) free(hap.s);
}


/*
 * SNVs
 */
START_TEST(check_get_haplotype_sequence__snvs_within_bounds) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(6, "A", "T"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "G"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(19, "A", "C"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAATAAAGAAAAAAAAC");
} END_TEST

START_TEST(check_get_haplotype_sequence__snv_outside_bounds) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(20, "A", "T"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, seq.s);
} END_TEST

/*
 * Insertions
 */
START_TEST(check_get_haplotype_sequence__insertion_within_bounds) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "AT"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAATAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__insertion_overlaps_bounds_start) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(0, "A", "AT"), NULL));

    seq.s[19] = '\0';
    get_haplotype_sequence(seq, &rs, (Range) {1, 20}, &hap);

    ck_assert_str_eq(hap.s, "TAAAAAAAAAAAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__insertion_overlaps_bound_end) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(19, "A", "ATT"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAAATT");
} END_TEST

/*
* Deletions
*/
START_TEST(check_get_haplotype_sequence__deletion_within_bounds) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "AA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__deletion_overlaps_bound_start) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(0, "AAA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, (Range) {1, 20}, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__deletion_overlaps_bound_end) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(18, "AAA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAA");
} END_TEST

/*
* Complex cases
*/
START_TEST(check_get_haplotype_sequence__consecutive_indels) {
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(10, "A", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "");
} END_TEST

Suite *haplotype_suite() {
    Suite *s = suite_create("Haplotype");

    TCase *tc_ghs = tcase_create("get_haplotype_sequence");

    suite_add_tcase(s, tc_ghs);
    tcase_add_checked_fixture(tc_ghs, checked_setup, checked_teardown);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__snvs_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__snv_outside_bounds);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_overlaps_bounds_start);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_overlaps_bound_end);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_overlaps_bound_start);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_overlaps_bound_end);

    return s;
}

