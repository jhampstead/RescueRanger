#include <check.h>
#include <cpliceai/range.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <string.h>

#include "klib/kvec.h"
#include "record.h"
#include "utils.h"

#include "check_utils.h"
#include "mocks.h"


/*
* VARIANT_TO_RANGE
*/
START_TEST(check_variant_to_range__snv) {
    bcf1_t *v = mock_bcf1_t(10, "A", "C");

    Range ref = variant_to_range(v, 0);
    Range alt = variant_to_range(v, 1);

    ck_assert_int_eq(ref.start, 10);
    ck_assert_int_eq(ref.end, 11);
    ck_assert_int_eq(range_len(ref), 1);

    ck_assert_int_eq(alt.start, 10);
    ck_assert_int_eq(alt.end, 11);
    ck_assert_int_eq(range_len(alt), 1);
} END_TEST

START_TEST(check_variant_to_range__insertion) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ACCC");

    Range alt = variant_to_range(v, 1);

    ck_assert_int_eq(alt.start, 10);
    ck_assert_int_eq(alt.end, 14);
    ck_assert_int_eq(range_len(alt), 4);
} END_TEST

START_TEST(check_variant_to_range__deletion) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    Range ref = variant_to_range(v, 0);

    ck_assert_int_eq(ref.start, 10);
    ck_assert_int_eq(ref.end, 17);
    ck_assert_int_eq(range_len(ref), 7);
} END_TEST


/*
* ALLELE_OFFSET
*/
START_TEST(check_allele_offset__snv_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "C");

    int offset = allele_offset(v, 0, 1, (Range) {0, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__snv_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "C");

    int offset = allele_offset(v, 0, 1, (Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__snv_outside_touches_bounds) {
    bcf1_t *v_start = mock_bcf1_t(19, "A", "C");

    int offset_start = allele_offset(v_start, 0, 1, (Range) {20, 30});
    ck_assert_int_eq(offset_start, 0);

    bcf1_t *v_end = mock_bcf1_t(29, "A", "C");

    int offset_end = allele_offset(v_end, 0, 1, (Range) {20, 30});
    ck_assert_int_eq(offset_end, 0);
} END_TEST

START_TEST(check_allele_offset__insertion_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ACCC");

    int offset = allele_offset(v, 0, 1, (Range) {0, 30});
    ck_assert_int_eq(offset, -3);
} END_TEST

START_TEST(check_allele_offset__insertion_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ACCC");

    int offset = allele_offset(v, 0, 1, (Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__deletion_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = allele_offset(v, 0, 1, (Range) {0, 30});
    ck_assert_int_eq(offset, 6);
} END_TEST

START_TEST(check_allele_offset__deletion_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = allele_offset(v, 0, 1, (Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__deletion_overlapping_bounds_start) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = allele_offset(v, 0, 1, (Range) {15, 30});
    ck_assert_int_eq(offset, 2);
} END_TEST

START_TEST(check_allele_offset__deletion_overlapping_bounds_end) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = allele_offset(v, 0, 1, (Range) {0, 15});
    ck_assert_int_eq(offset, 4);
} END_TEST

/*
* HAPLOTYPE OFFSET
*/
START_TEST(check_haplotype_offset) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(0, "A", "AT"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(4, "ATTTT", "A"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(9, "GC", "G"), NULL));

    int offset = haplotype_offset(&rs, (Range) {5, 10});

    ck_assert_int_eq(offset, 4);
} END_TEST

START_TEST(check_haplotype_offset__touching_indels) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(0, "A", "AT"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(1, "ATT", "A"), NULL));

    int offset = haplotype_offset(&rs, (Range) {0, 10} );

    ck_assert_int_eq(offset, 1);
} END_TEST


/*
* PREDICTION OFFSET
*/
START_TEST(check_prediction_offset) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(50, "A", "AT"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(54, "ATTTT", "A"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(59, "GCC", "G"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs);

    ck_assert_int_eq(offset, 5);
} END_TEST

START_TEST(check_prediction_offset__long_deletion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(59, "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "G"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs);

    ck_assert_int_eq(offset, 41);
} END_TEST

START_TEST(check_prediction_offset__long_insertion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, init_record(mock_bcf1_t(59, "G", "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs);

    ck_assert_int_eq(offset, -96);
} END_TEST


/*
* Match predictions to reference
*/
START_TEST(check_match_predictions_to_reference__insertion_within_bounds) {
    float predictions[306];
    for (int i = 0; i < 306; i++) {
        predictions[i] = 0.0;
    }

    predictions[150] = 0.7;
    predictions[151] = 0.15;
    predictions[152] = 0.15;
    predictions[153] = 0;
    predictions[154] = 0.9;
    predictions[155] = 0.1;

    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(50, "A", "AT"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, kv_A(rs, 0)->v, &rs, 303);

    ck_assert_float_eq(ref_preds[150], 0.7);
    ck_assert_float_eq(ref_preds[151], 0.9);
    ck_assert_float_eq(ref_preds[152], 0.15);
} END_TEST

START_TEST(check_match_predictions_to_reference__deletion_within_bounds) {
    float predictions[300];
    for (int i = 0; i < 300; i++) {
        predictions[i] = 0.1;
    }

    predictions[150] = 0.7;
    predictions[151] = 0.15;
    predictions[152] = 0.15;

    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(50, "AT", "A"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, kv_A(rs, 0)->v, &rs, 303);

    ck_assert_float_eq(ref_preds[150], 0.7);
    ck_assert_float_eq(ref_preds[151], 0.15);
    ck_assert_float_eq(ref_preds[152], 0.15);
    ck_assert_float_eq(ref_preds[153], 0);
    ck_assert_float_eq(ref_preds[154], 0);
    ck_assert_float_eq(ref_preds[155], 0);
} END_TEST

START_TEST(check_match_predictions_to_reference__deletion_overlaps_bounds_start) {
    float predictions[297];
    for (int i = 0; i < 297; i++) {
        predictions[i] = 0.1;
    }

    predictions[0] = 0.7;
    predictions[1] = 0.15;
    predictions[2] = 0.15;
    predictions[3] = 0.05;
    predictions[4] = 0.70;
    predictions[5] = 0.25;

    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(0, "ATT", "A"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, mock_bcf1_t(51, "A", "T"), &rs, 303);

    ck_assert_float_eq(ref_preds[0], 0.7);
    ck_assert_float_eq(ref_preds[1], 0.15);
    ck_assert_float_eq(ref_preds[2], 0.15);
    ck_assert_float_eq(ref_preds[3], 0);
    ck_assert_float_eq(ref_preds[4], 0);
    ck_assert_float_eq(ref_preds[5], 0);
} END_TEST

START_TEST(check_match_predictions_to_reference__deletion_overlaps_bounds_end) {
    float predictions[309];
    for (int i = 0; i < 309; i++) {
        predictions[i] = 0.0;
    }

    predictions[300] = 0.7;
    predictions[301] = 0.15;
    predictions[302] = 0.15;
    predictions[303] = 0;
    predictions[304] = 0.9;
    predictions[305] = 0.1;
    predictions[306] = 0.45;
    predictions[306] = 0.25;
    predictions[306] = 0.3;

    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, init_record(mock_bcf1_t(99, "A", "ATT"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, mock_bcf1_t(50, "A", "T"), &rs, 303);

    ck_assert_float_eq(ref_preds[0], 0.7);
    ck_assert_float_eq(ref_preds[1], 0.9);
    ck_assert_float_eq(ref_preds[2], 0.15);

} END_TEST

Suite * utils_suite(void) {
    Suite *s = suite_create("Utils");

    TCase *tc_vtr = tcase_create("variant_to_range");
    suite_add_tcase(s, tc_vtr);

    tcase_add_test(tc_vtr, check_variant_to_range__snv);
    tcase_add_test(tc_vtr, check_variant_to_range__insertion);
    tcase_add_test(tc_vtr, check_variant_to_range__deletion);

    TCase *tc_ao = tcase_create("allele_offset");
    suite_add_tcase(s, tc_ao);

    tcase_add_test(tc_ao, check_allele_offset__snv_within_bounds);
    tcase_add_test(tc_ao, check_allele_offset__snv_outside_bounds);
    tcase_add_test(tc_ao, check_allele_offset__snv_outside_touches_bounds);

    tcase_add_test(tc_ao, check_allele_offset__insertion_within_bounds);
    tcase_add_test(tc_ao, check_allele_offset__insertion_outside_bounds);

    tcase_add_test(tc_ao, check_allele_offset__deletion_within_bounds);
    tcase_add_test(tc_ao, check_allele_offset__deletion_outside_bounds);
    tcase_add_test(tc_ao, check_allele_offset__deletion_overlapping_bounds_start);
    tcase_add_test(tc_ao, check_allele_offset__deletion_overlapping_bounds_end);

    TCase *tc_ho = tcase_create("haplotype_offset");
    suite_add_tcase(s, tc_ho);

    tcase_add_test(tc_ho, check_haplotype_offset);
    tcase_add_test(tc_ho, check_haplotype_offset__touching_indels);

    TCase *tc_po = tcase_create("prediction_offset");
    suite_add_tcase(s, tc_po);

    tcase_add_test(tc_po, check_prediction_offset);
    tcase_add_test(tc_po, check_prediction_offset__long_deletion);
    tcase_add_test(tc_po, check_prediction_offset__long_insertion);

    TCase *tc_cmptf = tcase_create("match_predictions_to_reference");
    suite_add_tcase(s, tc_cmptf);

    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_overlaps_bounds_start);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_overlaps_bounds_end);


    return s;
}

