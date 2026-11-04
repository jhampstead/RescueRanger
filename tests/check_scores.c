#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <cpliceai/range.h>

#include "klib/kvec.h"
#include "record.h"
#include "scores.h"
#include "utils.h"

#include "check_scores.h"
#include "mocks.h"

// parse_prediction_pairs as it was before the anchor-free view, kept to compare against
static void old_parse_prediction_pairs(const float *ref_priors, const float *ref_posteriors, const float *alt_priors, const float *alt_posteriors, const Records *rs, const Range pred_range, PairPrediction **ref_alt_maps, size_t *ref_alt_map_size) {
    kvec_t(PairPrediction) map;
    kv_init(map);

    uint64_t curr_position = pred_range.start;
    int refs_idx = 0, alts_idx = 0;
    for (int i = 0; i < rs->n; i++) {
        const Record *r = kv_A(*rs, i);

        const Range ref = variant_to_range(r->v, 0);
        Range intersected_ref;
        if (intersect_ranges(ref, pred_range, &intersected_ref) == 0) continue;
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        while (curr_position < r->v->pos) {
            const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            const PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        const int offset = allele_offset(&r->alleles, pred_range);
        if (offset == 0) continue;

        if (curr_position == r->v->pos) {
            const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
            refs_idx++;
            const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
            alts_idx++;
            const PairPrediction pp = {ref_pred, alt_pred};
            kv_push(PairPrediction, map, pp);
            curr_position++;
        }

        const Prediction empty_pred = { true, 0, 0.0, 0.0 };
        if (offset > 0) {
            for (int i = 0; i < offset; i++) {
                const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
                refs_idx++;
                const PairPrediction pp = { ref_pred, empty_pred };
                kv_push(PairPrediction, map, pp);
                curr_position++;
            }
        } else {
            for (int i = 0; i < -offset; i++) {
                const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
                alts_idx++;
                const PairPrediction pp = { empty_pred, alt_pred };
                kv_push(PairPrediction, map, pp);
            }
        }
    }

    while (curr_position < pred_range.end) {
        const Prediction ref_pred = { false, curr_position, ref_priors[refs_idx], ref_posteriors[refs_idx] };
        refs_idx++;
        const Prediction alt_pred = { false, curr_position, alt_priors[alts_idx], alt_posteriors[alts_idx] };
        alts_idx++;
        const PairPrediction pp = {ref_pred, alt_pred};
        kv_push(PairPrediction, map, pp);
        curr_position++;
    }

    *ref_alt_maps = map.a;
    *ref_alt_map_size = map.n;
}

// Distinct values per array, so a value taken from the wrong index shows
#define SCORE_VALUES 300
typedef struct { float ref_priors[SCORE_VALUES], ref_posteriors[SCORE_VALUES], alt_priors[SCORE_VALUES], alt_posteriors[SCORE_VALUES]; } ScoreArrays;

static void fill_scores(ScoreArrays *s) {
    for (int i = 0; i < SCORE_VALUES; i++) {
        s->ref_priors[i] = 1000 + i;
        s->ref_posteriors[i] = 2000 + i;
        s->alt_priors[i] = 3000 + i;
        s->alt_posteriors[i] = 4000 + i;
    }
}

static void assert_same_prediction(const Prediction new, const Prediction old, const char *what, const size_t k) {
    ck_assert_msg(new.missing == old.missing && new.position == old.position && new.prior == old.prior && new.posterior == old.posterior, "%s, pair %zu", what, k);
}

static const Range pairs_window = { 100, 141 };

#define N_PAIR_VARIANTS 8
static const char *pair_variants[N_PAIR_VARIANTS][2] = {
    { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "N", "Nxyz" }, { "NN", "N" }, { "NNNN", "N" }, { "NNN", "Nx" }, { "NN", "Nxyz" }
};

// Every normalised variant kind at every position around both window edges, alone and before a second variant. Complex indels
// are left out: the old version paired them as a net insertion or deletion straight after the attached base, where they are
// paired base by base from the left now (the tests below).
START_TEST(check_parse_prediction_pairs__gives_the_old_pairs) {
    ScoreArrays s;
    fill_scores(&s);

    for (int i = 0; i < N_PAIR_VARIANTS; i++) {
        for (int pos = 95; pos <= 145; pos++) {
            for (int with_second = 0; with_second <= 1; with_second++) {
                Records rs; kv_init(rs);
                kv_push(Record *, rs, Record_init(mock_bcf1_t(pos, (char *) pair_variants[i][0], (char *) pair_variants[i][1]), NULL));
                const Alleles *first = &kv_A(rs, 0)->alleles;
                if (first->ref_len > 0 && first->alt_len > 0 && first->ref_len != first->alt_len) {
                    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
                    continue;
                }
                if (with_second && pos < 125) kv_push(Record *, rs, Record_init(mock_bcf1_t(130, "NN", "N"), NULL));

                PairPrediction *new = NULL, *old = NULL;
                size_t n_new = 0, n_old = 0;
                parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &new, &n_new);
                old_parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &old, &n_old);

                char what[64];
                snprintf(what, sizeof(what), "%s>%s at %d%s", pair_variants[i][0], pair_variants[i][1], pos, with_second ? " with a second variant" : "");
                ck_assert_msg(n_new == n_old, "%s: %zu pairs, old %zu", what, n_new, n_old);
                for (size_t k = 0; k < n_new; k++) {
                    assert_same_prediction(new[k].ref, old[k].ref, what, k);
                    assert_same_prediction(new[k].alt, old[k].alt, what, k);
                }

                free(new); free(old);
                for (size_t r = 0; r < kv_size(rs); r++) Record_destroy(kv_A(rs, r));
                kv_destroy(rs);
            }
        }
    }
} END_TEST

// An insertion attached to 110: the pair for 110 itself, then one pair per inserted base without a reference prediction
START_TEST(check_parse_prediction_pairs__insertion) {
    ScoreArrays s;
    fill_scores(&s);

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(110, "N", "Nxy"), NULL));

    PairPrediction *map = NULL;
    size_t n = 0;
    parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &map, &n);

    ck_assert_uint_eq(n, 41 + 2);
    ck_assert_uint_eq(map[10].ref.position, 110);
    ck_assert(!map[10].alt.missing);
    for (int k = 11; k <= 12; k++) {
        ck_assert(map[k].ref.missing);
        ck_assert_uint_eq(map[k].alt.position, 111);
        ck_assert_float_eq(map[k].alt.prior, s.alt_priors[k]);
    }
    ck_assert_uint_eq(map[13].ref.position, 111);
    ck_assert_float_eq(map[13].ref.prior, s.ref_priors[11]);
    ck_assert_float_eq(map[13].alt.prior, s.alt_priors[13]);

    free(map);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// A deletion attached to 110 removing 111 and 112: the pair for 110, then pairs without an ALT prediction for the removed bases
START_TEST(check_parse_prediction_pairs__deletion) {
    ScoreArrays s;
    fill_scores(&s);

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(110, "NNN", "N"), NULL));

    PairPrediction *map = NULL;
    size_t n = 0;
    parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &map, &n);

    ck_assert_uint_eq(n, 41);
    ck_assert(!map[10].alt.missing);
    for (int k = 11; k <= 12; k++) {
        ck_assert_uint_eq(map[k].ref.position, 100 + k);
        ck_assert(map[k].alt.missing);
    }
    ck_assert_uint_eq(map[13].ref.position, 113);
    ck_assert_float_eq(map[13].alt.prior, s.alt_priors[11]);

    free(map);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// A complex indel is paired base by base from the left: the paired bases are ordinary pairs, and the leftover bases come after
// them. NNN>Nx at 110: x pairs with the base at 111 (a substitution), and the base at 112 is deleted.
START_TEST(check_parse_prediction_pairs__complex_deletion) {
    ScoreArrays s;
    fill_scores(&s);

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(110, "NNN", "Nx"), NULL));

    PairPrediction *map = NULL;
    size_t n = 0;
    parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &map, &n);

    ck_assert_uint_eq(n, 41);

    // The attached base, then the substituted base: both ordinary pairs
    ck_assert_uint_eq(map[10].ref.position, 110);
    ck_assert_float_eq(map[10].alt.prior, s.alt_priors[10]);
    ck_assert_uint_eq(map[11].ref.position, 111);
    ck_assert(!map[11].alt.missing);
    ck_assert_float_eq(map[11].ref.prior, s.ref_priors[11]);
    ck_assert_float_eq(map[11].alt.prior, s.alt_priors[11]);

    // The deleted base has no ALT prediction
    ck_assert_uint_eq(map[12].ref.position, 112);
    ck_assert(map[12].alt.missing);
    ck_assert_float_eq(map[12].ref.prior, s.ref_priors[12]);

    // The next base is one ALT prediction behind
    ck_assert_uint_eq(map[13].ref.position, 113);
    ck_assert_float_eq(map[13].ref.prior, s.ref_priors[13]);
    ck_assert_float_eq(map[13].alt.prior, s.alt_priors[12]);

    free(map);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// NN>Nxyz at 110: x pairs with the base at 111, and y and z are inserted after it, without a reference prediction
START_TEST(check_parse_prediction_pairs__complex_insertion) {
    ScoreArrays s;
    fill_scores(&s);

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(110, "NN", "Nxyz"), NULL));

    PairPrediction *map = NULL;
    size_t n = 0;
    parse_prediction_pairs(s.ref_priors, s.ref_posteriors, s.alt_priors, s.alt_posteriors, &rs, pairs_window, &map, &n);

    ck_assert_uint_eq(n, 41 + 2);

    ck_assert_uint_eq(map[10].ref.position, 110);
    ck_assert_uint_eq(map[11].ref.position, 111);                       // The substituted base
    ck_assert(!map[11].alt.missing);
    ck_assert_float_eq(map[11].alt.prior, s.alt_priors[11]);

    for (int k = 12; k <= 13; k++) {                                    // The inserted bases
        ck_assert(map[k].ref.missing);
        ck_assert_uint_eq(map[k].alt.position, 112);
        ck_assert_float_eq(map[k].alt.prior, s.alt_priors[k]);
    }

    ck_assert_uint_eq(map[14].ref.position, 112);                       // The next base is two ALT predictions ahead
    ck_assert_float_eq(map[14].ref.prior, s.ref_priors[12]);
    ck_assert_float_eq(map[14].alt.prior, s.alt_priors[14]);

    free(map);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

Suite *scores_suite() {
    Suite *s = suite_create("Scores");

    TCase *tc_ppp = mocked_tcase("parse_prediction_pairs");
    suite_add_tcase(s, tc_ppp);

    tcase_add_test(tc_ppp, check_parse_prediction_pairs__gives_the_old_pairs);
    tcase_add_test(tc_ppp, check_parse_prediction_pairs__insertion);
    tcase_add_test(tc_ppp, check_parse_prediction_pairs__deletion);
    tcase_add_test(tc_ppp, check_parse_prediction_pairs__complex_deletion);
    tcase_add_test(tc_ppp, check_parse_prediction_pairs__complex_insertion);

    return s;
}
