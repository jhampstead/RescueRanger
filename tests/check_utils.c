#include <check.h>
#include <limits.h>
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
* RANGE_AROUND / RANGE_BEFORE
*/
START_TEST(check_range_around__default_distance) {
    Range r = range_around(50, 50);

    ck_assert_int_eq(r.start, 0);
    ck_assert_int_eq(r.end, 101);
    ck_assert_int_eq(range_len(r), 101);
} END_TEST

START_TEST(check_range_around__zero_offset) {
    Range r = range_around(50, 0);

    ck_assert_int_eq(r.start, 50);
    ck_assert_int_eq(r.end, 51);
    ck_assert_int_eq(range_len(r), 1);
} END_TEST

START_TEST(check_range_around__near_contig_start_is_not_clamped) {
    Range r = range_around(10, 50);

    ck_assert_int_eq(r.start, -40);
    ck_assert_int_eq(r.end, 61);
} END_TEST

START_TEST(check_range_before__default_distance) {
    Range r = range_before(50, 50);

    ck_assert_int_eq(r.start, 0);
    ck_assert_int_eq(r.end, 50);
    ck_assert_int_eq(range_len(r), 50);
} END_TEST

START_TEST(check_range_before__ends_at_pos_and_starts_like_range_around) {
    Range before = range_before(200, 30);
    Range around = range_around(200, 30);

    ck_assert_int_eq(before.end, 200);
    ck_assert_int_eq(before.start, around.start);
} END_TEST

START_TEST(check_range_before__near_contig_start_is_not_clamped) {
    Range r = range_before(10, 50);

    ck_assert_int_eq(r.start, -40);
    ck_assert_int_eq(r.end, 10);
} END_TEST


/*
* WINDOW WIDTHS (derived from distance)
*/
START_TEST(check_window_widths__default_distance) {
    ck_assert_int_eq(pred_width(50), 101);
    ck_assert_int_eq(model_width(50), 10101);
    ck_assert_int_eq(variant_position_in_sequence(50), 5050);
    ck_assert_int_eq(pred_array_len(50), 303);
} END_TEST

START_TEST(check_window_widths__other_distances) {
    ck_assert_int_eq(pred_width(1), 3);
    ck_assert_int_eq(model_width(1), 10003);
    ck_assert_int_eq(variant_position_in_sequence(1), 5001);
    ck_assert_int_eq(pred_array_len(1), 9);

    ck_assert_int_eq(pred_width(200), 401);
    ck_assert_int_eq(model_width(200), 10401);
    ck_assert_int_eq(variant_position_in_sequence(200), 5200);
    ck_assert_int_eq(pred_array_len(200), 1203);
} END_TEST

START_TEST(check_max_prediction_radius__derived_sizes_fit_in_int) {
    ck_assert_int_gt(MAX_PREDICTION_RADIUS, 0);
    ck_assert_int_gt(pred_array_len(MAX_PREDICTION_RADIUS), 0);
    ck_assert_int_gt(model_width(MAX_PREDICTION_RADIUS), 0);

    ck_assert((2LL * MAX_PREDICTION_RADIUS + 1) * NUM_SCORES <= INT_MAX);
} END_TEST

START_TEST(check_max_prediction_radius__is_tight) {
    const long long one_more = MAX_PREDICTION_RADIUS + 1LL;

    ck_assert((2LL * one_more + 1) * NUM_SCORES > INT_MAX);
} END_TEST

START_TEST(check_window_widths__consistent_with_range_around) {
    for (int d = 1; d <= 100; d += 7) {
        ck_assert_int_eq(range_len(range_around(1000, d)), pred_width(d));
    }
} END_TEST


/*
* VARIANT_FITS_WINDOW
*/
static void fill_bases(char *dst, const int n) {
    memset(dst, 'A', n);
    dst[n] = '\0';
}

START_TEST(check_variant_fits_window__snv) {
    ck_assert(variant_fits_window(mock_bcf1_t(10, "A", "C"), 50));
} END_TEST

START_TEST(check_variant_fits_window__alt_length_at_limit) {
    char at_limit[64], over_limit[64];
    fill_bases(at_limit, 50);
    fill_bases(over_limit, 51);

    ck_assert(variant_fits_window(mock_bcf1_t(10, "A", at_limit), 50));
    ck_assert(!variant_fits_window(mock_bcf1_t(10, "A", over_limit), 50));
} END_TEST

START_TEST(check_variant_fits_window__ref_length_at_limit) {
    char at_limit[64], over_limit[64];
    fill_bases(at_limit, 50);
    fill_bases(over_limit, 51);

    ck_assert(variant_fits_window(mock_bcf1_t(10, at_limit, "A"), 50));
    ck_assert(!variant_fits_window(mock_bcf1_t(10, over_limit, "A"), 50));
} END_TEST

START_TEST(check_variant_fits_window__depends_on_prediction_radius) {
    char alt[64];
    fill_bases(alt, 30);

    ck_assert(variant_fits_window(mock_bcf1_t(10, "A", alt), 50));
    ck_assert(!variant_fits_window(mock_bcf1_t(10, "A", alt), 20));
} END_TEST


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

    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "A", "AT"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(4, "ATTTT", "A"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(9, "GC", "G"), NULL));

    int offset = haplotype_offset(&rs, (Range) {5, 10});

    ck_assert_int_eq(offset, 4);
} END_TEST

START_TEST(check_haplotype_offset__touching_indels) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "A", "AT"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(1, "ATT", "A"), NULL));

    int offset = haplotype_offset(&rs, (Range) {0, 10} );

    ck_assert_int_eq(offset, 1);
} END_TEST


/*
* PREDICTION OFFSET
*/
START_TEST(check_prediction_offset) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "AT"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(54, "ATTTT", "A"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(59, "GCC", "G"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs, 50);

    ck_assert_int_eq(offset, 5);
} END_TEST

START_TEST(check_prediction_offset__long_deletion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(59, "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "G"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs, 50);

    ck_assert_int_eq(offset, 41);
} END_TEST

START_TEST(check_prediction_offset__long_insertion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(59, "G", "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), NULL));

    int offset = prediction_offset(kv_A(rs, 1)->v, &rs, 50);

    ck_assert_int_eq(offset, -96);
} END_TEST

START_TEST(check_prediction_offset__depends_on_prediction_radius) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(80, "AT", "A"), NULL)); // 30 bases from the variant at 50

    // Window {0, 101} contains the deletion at 80, which shortens the sequence by 1
    ck_assert_int_eq(prediction_offset(kv_A(rs, 0)->v, &rs, 50), 1);
    ck_assert_int_eq(prediction_offset(kv_A(rs, 0)->v, &rs, 100), 1);

    // Window {30, 71} does not
    ck_assert_int_eq(prediction_offset(kv_A(rs, 0)->v, &rs, 20), 0);
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

    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "AT"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, range_around(50, 50), &rs, 303);

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

    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "AT", "A"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, range_around(50, 50), &rs, 303); // NOTE: predictions has 300 floats, 303 are read

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

    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "ATT", "A"), NULL));

    float *ref_preds = match_predictions_to_reference(predictions, range_around(51, 50), &rs, 303); // NOTE: suspect, window centre 51 vs record at pos 0

    ck_assert_float_eq(ref_preds[0], 0.7);
    ck_assert_float_eq(ref_preds[1], 0.15);
    ck_assert_float_eq(ref_preds[2], 0.15);
    ck_assert_float_eq(ref_preds[3], 0);
    ck_assert_float_eq(ref_preds[4], 0);
    ck_assert_float_eq(ref_preds[5], 0);
} END_TEST

Suite * utils_suite(void) {
    Suite *s = suite_create("Utils");

    TCase *tc_rng = tcase_create("range_helpers");
    suite_add_tcase(s, tc_rng);

    tcase_add_test(tc_rng, check_range_around__default_distance);
    tcase_add_test(tc_rng, check_range_around__zero_offset);
    tcase_add_test(tc_rng, check_range_around__near_contig_start_is_not_clamped);
    tcase_add_test(tc_rng, check_range_before__default_distance);
    tcase_add_test(tc_rng, check_range_before__ends_at_pos_and_starts_like_range_around);
    tcase_add_test(tc_rng, check_range_before__near_contig_start_is_not_clamped);

    TCase *tc_ww = tcase_create("window_widths");
    suite_add_tcase(s, tc_ww);

    tcase_add_test(tc_ww, check_window_widths__default_distance);
    tcase_add_test(tc_ww, check_window_widths__other_distances);
    tcase_add_test(tc_ww, check_window_widths__consistent_with_range_around);
    tcase_add_test(tc_ww, check_max_prediction_radius__derived_sizes_fit_in_int);
    tcase_add_test(tc_ww, check_max_prediction_radius__is_tight);

    TCase *tc_vfw = tcase_create("variant_fits_window");
    suite_add_tcase(s, tc_vfw);

    tcase_add_test(tc_vfw, check_variant_fits_window__snv);
    tcase_add_test(tc_vfw, check_variant_fits_window__alt_length_at_limit);
    tcase_add_test(tc_vfw, check_variant_fits_window__ref_length_at_limit);
    tcase_add_test(tc_vfw, check_variant_fits_window__depends_on_prediction_radius);

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
    tcase_add_test(tc_po, check_prediction_offset__depends_on_prediction_radius);

    TCase *tc_cmptf = tcase_create("match_predictions_to_reference");
    suite_add_tcase(s, tc_cmptf);

    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_overlaps_bounds_start);


    return s;
}

