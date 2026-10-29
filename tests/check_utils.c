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

// Destroys the records and the vector that holds them. The bcf1_t a record was made from belongs to the mocks, not to the record
// (Record_init copies it).
static void destroy_records(Records *rs) {
    for (size_t i = 0; i < kv_size(*rs); i++) Record_destroy(kv_A(*rs, i));
    kv_destroy(*rs);
}

/*
* MOCKS
*/
START_TEST(check_mocks_free_all__frees_what_mock_bcf1_t_made) {
    ck_assert_uint_eq(mocks_live_count(), 0);

    mock_bcf1_t(10, "A", "C");
    mock_bcf1_t(20, "AT", "A");
    mock_bcf1_t(30, "A", "AT");
    ck_assert_uint_eq(mocks_live_count(), 3);

    mocks_free_all();
    ck_assert_uint_eq(mocks_live_count(), 0);

    // Nothing left to free, and the registry works again
    mocks_free_all();
    mock_bcf1_t(10, "A", "C");
    ck_assert_uint_eq(mocks_live_count(), 1);
} END_TEST

START_TEST(check_mocks_free_all__keeps_a_record_made_from_a_mock) {
    Record *r = Record_init(mock_bcf1_t(10, "A", "C"), NULL);

    mocks_free_all();

    ck_assert_int_eq(r->alleles.pos, 10);
    ck_assert_int_eq(r->alleles.ref_len, 1);
    Record_destroy(r);
} END_TEST


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

// variant_fits_window works on the anchor-free view: only the bases that differ count, not the shared anchor base
static bool fits_window(const int pos, const char *ref, const char *alt, const int prediction_radius) {
    const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
    return variant_fits_window(&a, prediction_radius);
}

START_TEST(check_variant_fits_window__snv) {
    ck_assert(fits_window(10, "A", "C", 50));
} END_TEST

START_TEST(check_variant_fits_window__alt_length_at_limit) {
    // A substitution has no shared prefix, so its ALT length is the VCF one
    char at_limit[64], over_limit[64];
    fill_bases(at_limit, 50);
    fill_bases(over_limit, 51);
    at_limit[0] = over_limit[0] = 'C';

    ck_assert(fits_window(10, "A", at_limit, 50));
    ck_assert(!fits_window(10, "A", over_limit, 50));
} END_TEST

START_TEST(check_variant_fits_window__ref_length_at_limit) {
    char at_limit[64], over_limit[64];
    fill_bases(at_limit, 50);
    fill_bases(over_limit, 51);
    at_limit[0] = over_limit[0] = 'C';

    ck_assert(fits_window(10, at_limit, "A", 50));
    ck_assert(!fits_window(10, over_limit, "A", 50));
} END_TEST

START_TEST(check_variant_fits_window__insertion_does_not_count_anchor) {
    // Insertion of 50 bases after the anchor A: VCF ALT is 51 bases long, the view's is 50
    char inserted[64], over_limit[64];
    fill_bases(inserted, 51);
    fill_bases(over_limit, 52);

    ck_assert(fits_window(10, "A", inserted, 50));
    ck_assert(!fits_window(10, "A", over_limit, 50));
} END_TEST

START_TEST(check_variant_fits_window__deletion_does_not_count_anchor) {
    // Deletion of 50 bases after the anchor A: VCF REF is 51 bases long, the view's is 50
    char deleted[64], over_limit[64];
    fill_bases(deleted, 51);
    fill_bases(over_limit, 52);

    ck_assert(fits_window(10, deleted, "A", 50));
    ck_assert(!fits_window(10, over_limit, "A", 50));
} END_TEST

START_TEST(check_variant_fits_window__depends_on_prediction_radius) {
    char alt[64];
    fill_bases(alt, 30);
    alt[0] = 'C';

    ck_assert(fits_window(10, "A", alt, 50));
    ck_assert(!fits_window(10, "A", alt, 20));
} END_TEST


/*
* VARIANT_IN_REGION
*
* Whether a variant belongs to a region such as a gene: some base it changes lies in it. The base an indel is attached to is
* not changed, so a deletion anchored on the last base of the region, which only deletes bases after it, is not in it.
*/
static bool in_region(const int64_t pos, const char *ref, const char *alt, const Range region) {
    const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
    return variant_in_region(&a, region);
}

static const Range region_100_200 = { 100, 200 };

START_TEST(check_variant_in_region__substitutions) {
    ck_assert(!in_region(99, "N", "x", region_100_200));
    ck_assert(in_region(100, "N", "x", region_100_200));
    ck_assert(in_region(199, "N", "x", region_100_200));
    ck_assert(!in_region(200, "N", "x", region_100_200));

    // An MNV is in the region when any of its bases is
    ck_assert(!in_region(98, "NN", "xy", region_100_200));
    ck_assert(in_region(99, "NN", "xy", region_100_200));
    ck_assert(in_region(199, "NN", "xy", region_100_200));
    ck_assert(!in_region(200, "NN", "xy", region_100_200));
} END_TEST

START_TEST(check_variant_in_region__deletions_count_by_the_deleted_bases) {
    // NN>N at p deletes p + 1
    ck_assert(!in_region(98, "NN", "N", region_100_200)); // Deletes 99
    ck_assert(in_region(99, "NN", "N", region_100_200));  // Deletes 100, the anchor is outside
    ck_assert(in_region(198, "NN", "N", region_100_200)); // Deletes 199
    ck_assert(!in_region(199, "NN", "N", region_100_200)); // Deletes 200, the anchor is inside

    // NNNN>N at p deletes p + 1 to p + 3
    ck_assert(in_region(97, "NNNN", "N", region_100_200));  // Deletes 98, 99, 100
    ck_assert(in_region(197, "NNNN", "N", region_100_200)); // Deletes 198, 199, 200
    ck_assert(!in_region(199, "NNNN", "N", region_100_200)); // Deletes 200, 201, 202
} END_TEST

START_TEST(check_variant_in_region__complex_indels_count_by_the_replaced_bases) {
    ck_assert(!in_region(199, "NN", "Nxyz", region_100_200)); // Replaces 200
    ck_assert(in_region(99, "NN", "Nxyz", region_100_200));   // Replaces 100
    ck_assert(!in_region(199, "NNN", "Nx", region_100_200));  // Replaces 200, 201
    ck_assert(in_region(198, "NNN", "Nx", region_100_200));   // Replaces 199, 200
} END_TEST

// An insertion changes no base. It lies between two bases and belongs to the region only when both are in it, so an insertion at
// either edge of the region (before the first base or after the last) is not in it. It used to follow the base it is attached
// to, which put the insertion after the last base in and the one before the first base out.
START_TEST(check_variant_in_region__an_insertion_needs_both_neighbours_in_the_region) {
    ck_assert(!in_region(98, "N", "Nx", region_100_200));  // Between 98 and 99
    ck_assert(!in_region(99, "N", "Nx", region_100_200));  // Between 99 and 100: before the first base
    ck_assert(in_region(100, "N", "Nx", region_100_200));  // Between 100 and 101
    ck_assert(in_region(198, "N", "Nx", region_100_200));  // Between 198 and 199
    ck_assert(!in_region(199, "N", "Nx", region_100_200)); // Between 199 and 200: after the last base
    ck_assert(!in_region(200, "N", "Nx", region_100_200));

    ck_assert(in_region(150, "N", "Nxyz", region_100_200));
    ck_assert(!in_region(199, "N", "Nxyz", region_100_200));
} END_TEST

START_TEST(check_variant_in_region__differs_from_the_footprint_only_by_the_anchor) {
    // The footprint includes the anchor, so it overlaps the region where the changed bases do not: an indel attached to the
    // last base of the region whose change lies after it
    const Alleles del = Alleles_from_variant(mock_bcf1_t(199, "NN", "N"));
    ck_assert(variant_overlaps(&del, region_100_200));
    ck_assert(!variant_in_region(&del, region_100_200));

    // Everywhere else they agree
    for (int pos = 95; pos <= 205; pos++) {
        char *kinds[][2] = { { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "NN", "N" }, { "NNNN", "N" } };
        for (int k = 0; k < 5; k++) {
            const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, kinds[k][0], kinds[k][1]));
            const int anchored_only = variant_overlaps(&a, region_100_200) && !variant_in_region(&a, region_100_200);
            const int attached_to_the_last_base = a.prefix_len > 0 && a.pos == region_100_200.end;
            ck_assert_msg(anchored_only == attached_to_the_last_base, "%s>%s at %d", kinds[k][0], kinds[k][1], pos);
            ck_assert_msg(!(variant_in_region(&a, region_100_200) && !variant_overlaps(&a, region_100_200)), "%s>%s at %d", kinds[k][0], kinds[k][1], pos);
        }
    }
} END_TEST

/*
* HAPLOTYPE_SHIFT and REGION_TO_HAPLOTYPE
*
* How far the haplotype moves a coordinate b of a gene [100, 200): the bases deleted before it minus the bases inserted before
* it. The pieces of a variant are its paired bases (substitutions, which move nothing), the REF bases left over (deleted) and
* the ALT bases left over (inserted at the point after the last REF base). An insertion exactly at b counts as before b only
* when insertions_at_b is set, because which side of a region border it belongs to depends on the region.
*/
static int64_t shift_at(const int pos, const char *ref, const char *alt, const int64_t b, const bool insertions_at_b) {
    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(pos, ref, alt), NULL));

    const int64_t shift = haplotype_shift(&rs, region_100_200, b, insertions_at_b);

    Record_destroy(kv_A(rs, 0));
    kv_destroy(rs);
    return shift;
}

START_TEST(check_haplotype_shift__substitutions_move_nothing) {
    for (int b = 99; b <= 152; b++) {
        ck_assert_int_eq(shift_at(150, "N", "x", b, true), 0);
        ck_assert_int_eq(shift_at(109, "NN", "xy", b, true), 0); // An MNV that straddles any bound moves nothing
        ck_assert_int_eq(shift_at(149, "NN", "xy", b, false), 0);
    }
} END_TEST

START_TEST(check_haplotype_shift__a_deletion_moves_what_lies_after_its_deleted_bases) {
    // NN>N at 120 deletes 121
    ck_assert_int_eq(shift_at(120, "NN", "N", 121, true), 0);  // The deleted base is not before 121
    ck_assert_int_eq(shift_at(120, "NN", "N", 122, true), 1);
    ck_assert_int_eq(shift_at(120, "NN", "N", 150, false), 1);

    // NNNN>N at 97 deletes 98, 99, 100: only 100 is in the gene
    ck_assert_int_eq(shift_at(97, "NNNN", "N", 100, true), 0);
    ck_assert_int_eq(shift_at(97, "NNNN", "N", 101, true), 1);
    ck_assert_int_eq(shift_at(97, "NNNN", "N", 150, true), 1);

    // NNNN>N at 107 deletes 108, 109, 110
    ck_assert_int_eq(shift_at(107, "NNNN", "N", 110, true), 2);
    ck_assert_int_eq(shift_at(107, "NNNN", "N", 111, true), 3);
} END_TEST

START_TEST(check_haplotype_shift__an_insertion_at_b_counts_only_when_asked) {
    // N>Nx at 120 inserts between 120 and 121, the point 121
    ck_assert_int_eq(shift_at(120, "N", "Nx", 120, true), 0);
    ck_assert_int_eq(shift_at(120, "N", "Nx", 121, true), -1);
    ck_assert_int_eq(shift_at(120, "N", "Nx", 121, false), 0);
    ck_assert_int_eq(shift_at(120, "N", "Nx", 122, false), -1);
    ck_assert_int_eq(shift_at(120, "N", "Nxyz", 122, false), -3);

    // An insertion at the gene end is not in the gene and moves nothing
    ck_assert_int_eq(shift_at(199, "N", "Nx", 200, true), 0);
} END_TEST

START_TEST(check_haplotype_shift__complex_indels_are_a_substitution_and_a_leftover) {
    // NNN>Nx at 109: x pairs with 110, 111 is deleted
    ck_assert_int_eq(shift_at(109, "NNN", "Nx", 111, true), 0);
    ck_assert_int_eq(shift_at(109, "NNN", "Nx", 112, true), 1);

    // NN>Nxyz at 109: x pairs with 110, y and z are inserted at the point 111
    ck_assert_int_eq(shift_at(109, "NN", "Nxyz", 111, true), -2);
    ck_assert_int_eq(shift_at(109, "NN", "Nxyz", 111, false), 0);
    ck_assert_int_eq(shift_at(109, "NN", "Nxyz", 112, false), -2);
} END_TEST

START_TEST(check_haplotype_shift__the_records_add_up) {
    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(120, "NN", "N"), NULL));   // Deletes 121
    kv_push(Record *, rs, Record_init(mock_bcf1_t(130, "N", "Nxy"), NULL));  // Inserts 2 at 131

    ck_assert_int_eq(haplotype_shift(&rs, region_100_200, 140, true), -1);

    for (int i = 0; i < 2; i++) Record_destroy(kv_A(rs, i));
    kv_destroy(rs);
} END_TEST

static Range region_in_haplotype(const int pos, const char *ref, const char *alt, const Range region, const bool is_utr) {
    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(pos, ref, alt), NULL));

    const Range hap = region_to_haplotype(&rs, region_100_200, region, is_utr);

    Record_destroy(kv_A(rs, 0));
    kv_destroy(rs);
    return hap;
}

START_TEST(check_region_to_haplotype__an_insertion_at_a_border_goes_to_the_side_without_the_feature) {
    const Range intron = { 130, 150 }, utr5 = { 100, 110 }, utr3 = { 180, 200 };
    Range hap;

    // Between the intron and exon 2 (after 149): the insertion is in exon 2, so the intron does not grow
    hap = region_in_haplotype(149, "N", "Nx", intron, false);
    ck_assert_int_eq(hap.start, 130); ck_assert_int_eq(hap.end, 150);

    // Between exon 1 and the intron (after 129): the insertion is in exon 1, so the intron moves up
    hap = region_in_haplotype(129, "N", "Nx", intron, false);
    ck_assert_int_eq(hap.start, 131); ck_assert_int_eq(hap.end, 151);

    // Between the 5' UTR and the CDS (after 109): the insertion is in the UTR, so the UTR grows
    hap = region_in_haplotype(109, "N", "Nx", utr5, true);
    ck_assert_int_eq(hap.start, 100); ck_assert_int_eq(hap.end, 111);

    // Between the CDS and the 3' UTR (after 179): the insertion is in the UTR, so the UTR starts where it did and grows
    hap = region_in_haplotype(179, "N", "Nx", utr3, true);
    ck_assert_int_eq(hap.start, 180); ck_assert_int_eq(hap.end, 201);
} END_TEST

// An intron of length 0 (touching exons) has no inside: an insertion at its position cannot be before its start and after its
// end, which used to turn the region around ([131, 130)) and made the spliced sequence copy a base twice
START_TEST(check_region_to_haplotype__an_empty_region_stays_empty) {
    const Range touching = { 130, 130 };

    // Between the two exons that touch (after 129)
    Range hap = region_in_haplotype(129, "N", "Nx", touching, false);
    ck_assert_int_eq(hap.start, hap.end);

    // Before and after it
    hap = region_in_haplotype(120, "N", "Nxyz", touching, false);
    ck_assert_int_eq(hap.start, 133); ck_assert_int_eq(hap.end, 133);

    hap = region_in_haplotype(140, "N", "Nxyz", touching, false);
    ck_assert_int_eq(hap.start, 130); ck_assert_int_eq(hap.end, 130);

    // A deletion before it moves both ends the same
    hap = region_in_haplotype(120, "NN", "N", touching, false);
    ck_assert_int_eq(hap.start, 129); ck_assert_int_eq(hap.end, 129);

    // An empty UTR the same way
    hap = region_in_haplotype(110, "N", "Nx", (Range) { 100, 100 }, true);
    ck_assert_int_eq(hap.start, 100); ck_assert_int_eq(hap.end, 100);
} END_TEST

START_TEST(check_region_to_haplotype__a_substitution_across_a_border_changes_nothing) {
    Range hap = region_in_haplotype(109, "NN", "xy", (Range) { 100, 110 }, true);
    ck_assert_int_eq(hap.start, 100); ck_assert_int_eq(hap.end, 110);

    hap = region_in_haplotype(149, "NN", "xy", (Range) { 130, 150 }, false);
    ck_assert_int_eq(hap.start, 130); ck_assert_int_eq(hap.end, 150);
} END_TEST

/*
* VARIANT_POSITION
*/
static int64_t position_of(const int64_t pos, const char *ref, const char *alt) {
    const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
    return variant_position(&a);
}

START_TEST(check_variant_position__substitutions_count_at_their_first_base) {
    ck_assert_int_eq(position_of(10, "A", "C"), 10);
    ck_assert_int_eq(position_of(10, "AT", "CG"), 10); // MNV
} END_TEST

START_TEST(check_variant_position__a_deletion_counts_at_its_first_deleted_base) {
    ck_assert_int_eq(position_of(10, "AT", "A"), 11);
    ck_assert_int_eq(position_of(10, "ATGC", "A"), 11);
} END_TEST

START_TEST(check_variant_position__an_insertion_counts_at_the_base_it_is_attached_to) {
    // No REF bases are changed, so it counts at the base before it
    ck_assert_int_eq(position_of(10, "A", "AT"), 10);
    ck_assert_int_eq(position_of(10, "A", "ATGC"), 10);
} END_TEST

START_TEST(check_variant_position__a_complex_indel_counts_at_its_first_replaced_base) {
    ck_assert_int_eq(position_of(10, "ATG", "AC"), 11);
    ck_assert_int_eq(position_of(10, "AT", "ACG"), 11);
} END_TEST

START_TEST(check_variant_position__window_centre_is_still_the_attached_base) {
    // The footprint, used for window centres and prediction matching, keeps starting at the attached base
    const Alleles del = Alleles_from_variant(mock_bcf1_t(10, "AT", "A"));
    ck_assert_int_eq(variant_footprint(&del).start, 10);
    ck_assert_int_eq(variant_position(&del), 11);
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
*
* The expected values are the ones the version on bcf1_t gave, before the anchor-free view.
*/
static int offset_of(const bcf1_t *v, const Range bounds) {
    const Alleles a = Alleles_from_variant(v);
    return allele_offset(&a, bounds);
}

START_TEST(check_allele_offset__snv_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "C");

    int offset = offset_of(v,(Range) {0, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__snv_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "C");

    int offset = offset_of(v,(Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__snv_outside_touches_bounds) {
    bcf1_t *v_start = mock_bcf1_t(19, "A", "C");

    int offset_start = offset_of(v_start, (Range) {20, 30});
    ck_assert_int_eq(offset_start, 0);

    bcf1_t *v_end = mock_bcf1_t(29, "A", "C");

    int offset_end = offset_of(v_end, (Range) {20, 30});
    ck_assert_int_eq(offset_end, 0);
} END_TEST

START_TEST(check_allele_offset__insertion_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ACCC");

    int offset = offset_of(v,(Range) {0, 30});
    ck_assert_int_eq(offset, -3);
} END_TEST

START_TEST(check_allele_offset__insertion_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ACCC");

    int offset = offset_of(v,(Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__deletion_within_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = offset_of(v,(Range) {0, 30});
    ck_assert_int_eq(offset, 6);
} END_TEST

START_TEST(check_allele_offset__deletion_outside_bounds) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = offset_of(v,(Range) {20, 30});
    ck_assert_int_eq(offset, 0);
} END_TEST

START_TEST(check_allele_offset__deletion_overlapping_bounds_start) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = offset_of(v,(Range) {15, 30});
    ck_assert_int_eq(offset, 2);
} END_TEST

START_TEST(check_allele_offset__deletion_overlapping_bounds_end) {
    bcf1_t *v = mock_bcf1_t(10, "ATCGTCG", "A");

    int offset = offset_of(v,(Range) {0, 15});
    ck_assert_int_eq(offset, 4);
} END_TEST

/*
* The anchor-free view gives the old answers
*/

// allele_offset as it was on bcf1_t, with the anchor, kept to compare against
static int old_allele_offset(const bcf1_t *v, const Range bounds) {
    Range ref = variant_to_range(v, 0);
    Range ref_range;
    if (intersect_ranges(ref, bounds, &ref_range) == 0 || range_len(ref_range) == 0) return 0;

    Range alt = variant_to_range(v, 1);
    Range intersected_alt;
    int alt_len = 0;
    if (intersect_ranges((Range) {v->pos, v->pos+1}, bounds, &intersected_alt) != 0 && range_len(intersected_alt) > 0) {
        alt_len = range_len(alt);
    }

    return range_len(ref_range) - alt_len;
}

// The old overlap test: the record's REF range, anchor included, shares a base with the range
static bool old_overlaps(const bcf1_t *v, const Range range) {
    Range intersection = {0};
    intersect_ranges(range, variant_to_range(v, 0), &intersection);
    return range_len(intersection) != 0;
}

// Normalised variants (at most one shared base), all at position 10
#define N_GRID_VARIANTS 9
static const char *grid_variants[N_GRID_VARIANTS][2] = {
    { "A", "C" }, { "AT", "GC" }, { "A", "AT" }, { "A", "ATTT" }, { "AT", "A" }, { "ATTT", "A" }, { "ATG", "AC" }, { "A", "A" },
    { "A", "<DEL>" }
};

START_TEST(check_variant_footprint__is_the_old_ref_range_for_normalised_variants) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        bcf1_t *v = mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]);
        const Alleles a = Alleles_from_variant(v);

        const Range footprint = variant_footprint(&a);
        const Range old = variant_to_range(v, 0);
        ck_assert_msg(footprint.start == old.start && footprint.end == old.end, "%s>%s", grid_variants[i][0], grid_variants[i][1]);
    }
} END_TEST

START_TEST(check_variant_overlaps__gives_the_old_answer_for_every_window) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        bcf1_t *v = mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]);
        const Alleles a = Alleles_from_variant(v);

        for (int s = 4; s <= 17; s++) {
            for (int e = s - 2; e <= 18; e++) {
                const Range window = { s, e };
                ck_assert_msg(variant_overlaps(&a, window) == old_overlaps(v, window), "%s>%s in [%d, %d)", grid_variants[i][0], grid_variants[i][1], s, e);
            }
        }
    }
} END_TEST

// The offset of a variant in a window, counted base by base from the pieces of the variant: the REF bases left over after pairing
// them with the ALT bases from the left are deleted, and count when they lie in the window; the ALT bases left over are inserted
// after the last paired base, and count when the base they are attached to lies in the window. Paired bases move nothing.
static int piecewise_offset(const Alleles *a, const Range window) {
    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
    int offset = 0;

    for (int64_t b = a->pos + paired; b < a->pos + a->ref_len; b++) {
        if (b >= window.start && b < window.end) offset++;
    }

    if (a->alt_len > a->ref_len) {
        const int64_t attached = a->pos + paired - 1;
        if (attached >= window.start && attached < window.end) offset -= a->alt_len - a->ref_len;
    }

    return offset;
}

// A window edge can cut through an MNV or a complex indel (or a symbolic allele, which reads as one), and the old formula
// counted the REF bases inside the window as removed whatever the ALT was. Those get the piecewise value instead. An empty or
// inverted window holds nothing, where the old code swapped the ends.
START_TEST(check_allele_offset__gives_the_old_answer_for_every_window) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        bcf1_t *v = mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]);
        const Alleles a = Alleles_from_variant(v);
        const bool can_be_cut = a.ref_len > 0 && a.alt_len > 0 && !(a.ref_len == 1 && a.alt_len == 1);

        for (int s = 4; s <= 17; s++) {
            for (int e = s - 2; e <= 18; e++) {
                const Range window = { s, e };
                const int offset = allele_offset(&a, window);

                if (e <= s) {
                    ck_assert_msg(offset == 0, "%s>%s in [%d, %d): %d, expected 0 for an empty window", grid_variants[i][0], grid_variants[i][1], s, e, offset);
                } else if (can_be_cut) {
                    ck_assert_msg(offset == piecewise_offset(&a, window), "%s>%s in [%d, %d): %d, expected %d", grid_variants[i][0], grid_variants[i][1], s, e, offset, piecewise_offset(&a, window));
                } else {
                    ck_assert_msg(offset == old_allele_offset(v, window), "%s>%s in [%d, %d): %d, old %d", grid_variants[i][0], grid_variants[i][1], s, e, offset, old_allele_offset(v, window));
                }
            }
        }
    }
} END_TEST

// A window edge through a variant: only the part of it in the window counts, and a substitution never moves anything
START_TEST(check_allele_offset__a_window_edge_through_a_variant) {
    // AT>GC at 10 replaces 10 and 11 one for one
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "GC"), (Range) { 4, 11 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "GC"), (Range) { 11, 18 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "GC"), (Range) { 10, 12 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "GC"), (Range) { 4, 18 }), 0);

    // ATG>AC at 10: C replaces the T at 11 and the G at 12 is deleted
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "AC"), (Range) { 4, 12 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "AC"), (Range) { 4, 13 }), 1);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "AC"), (Range) { 12, 18 }), 1);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "AC"), (Range) { 11, 18 }), 1);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "AC"), (Range) { 13, 18 }), 0);

    // AT>ACCC at 10: C replaces the T at 11 and CC is inserted after it
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "ACCC"), (Range) { 4, 11 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "ACCC"), (Range) { 4, 12 }), -2);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "ACCC"), (Range) { 11, 12 }), -2);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "AT", "ACCC"), (Range) { 12, 18 }), 0);
} END_TEST

// variant_shift(a, b, insertions_at_b): what one variant moves the coordinate b by, the REF bases deleted before b minus the ALT
// bases inserted at points before b (at b only with insertions_at_b). There is no gene here: nothing is clipped.
static int64_t shift_of(const char *ref, const char *alt, const int pos, const int64_t b, const bool at_b) {
    const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
    return variant_shift(&a, b, at_b);
}

START_TEST(check_variant_shift__substitutions_move_nothing) {
    for (int b = 0; b <= 30; b++) {
        ck_assert_int_eq(shift_of("A", "C", 10, b, true), 0);
        ck_assert_int_eq(shift_of("AT", "GC", 10, b, false), 0); // An MNV, whatever b cuts through
    }
} END_TEST

START_TEST(check_variant_shift__deleted_bases_move_what_is_after_them) {
    // AT>A at 10 deletes 11; ATTT>A at 10 deletes 11 to 13
    ck_assert_int_eq(shift_of("AT", "A", 10, 11, true), 0);   // The deleted base is not before 11
    ck_assert_int_eq(shift_of("AT", "A", 10, 12, true), 1);
    ck_assert_int_eq(shift_of("AT", "A", 10, 500, false), 1);
    ck_assert_int_eq(shift_of("ATTT", "A", 10, 12, true), 1);
    ck_assert_int_eq(shift_of("ATTT", "A", 10, 13, true), 2);
    ck_assert_int_eq(shift_of("ATTT", "A", 10, 14, false), 3);
    ck_assert_int_eq(shift_of("ATTT", "A", 10, 5, false), 0);  // Before the variant
} END_TEST

START_TEST(check_variant_shift__an_insertion_at_b_counts_only_when_asked) {
    // A>ACC at 10 inserts two bases between 10 and 11, the point 11
    ck_assert_int_eq(shift_of("A", "ACC", 10, 10, true), 0);
    ck_assert_int_eq(shift_of("A", "ACC", 10, 11, true), -2);
    ck_assert_int_eq(shift_of("A", "ACC", 10, 11, false), 0);
    ck_assert_int_eq(shift_of("A", "ACC", 10, 12, false), -2);
} END_TEST

START_TEST(check_variant_shift__complex_indels_are_a_substitution_and_a_leftover) {
    // ATG>AC at 10: C pairs with the T at 11, the G at 12 is deleted
    ck_assert_int_eq(shift_of("ATG", "AC", 10, 12, true), 0);
    ck_assert_int_eq(shift_of("ATG", "AC", 10, 13, true), 1);

    // AT>ACCC at 10: C pairs with the T at 11 and CC is inserted at the point 12
    ck_assert_int_eq(shift_of("AT", "ACCC", 10, 12, true), -2);
    ck_assert_int_eq(shift_of("AT", "ACCC", 10, 12, false), 0);
    ck_assert_int_eq(shift_of("AT", "ACCC", 10, 13, false), -2);
} END_TEST

// The offset in a window [s, e) is the difference of the shifts at its ends: the deleted bases in it, minus the inserted bases
// whose attached base is in it (their point is then in (s, e], which is what the shift at e counts and the shift at s does not).
// piecewise_offset counts that base by base, so this does not depend on how allele_offset is written.
START_TEST(check_variant_shift__the_difference_over_a_window_is_the_offset) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        const Alleles a = Alleles_from_variant(mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]));

        for (int s = 4; s <= 17; s++) {
            for (int e = s + 1; e <= 18; e++) {
                const Range window = { s, e };
                ck_assert_msg(variant_shift(&a, e, true) - variant_shift(&a, s, true) == piecewise_offset(&a, window),
                    "%s>%s in [%d, %d)", grid_variants[i][0], grid_variants[i][1], s, e);
            }
        }
    }
} END_TEST

// The prediction code adds and subtracts the offsets of neighbouring windows (pred_offset - start_offset is the offset of the
// part of the window from the centre on), so the offset has to be additive over adjacent windows, whatever cuts through a
// variant. And a window that holds the whole footprint of a variant gives the offset of the variant, which is what the offset of
// the whole gene is compared with.
START_TEST(check_allele_offset__is_additive_over_adjacent_windows) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        const Alleles a = Alleles_from_variant(mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]));

        for (int s = 4; s <= 17; s++) {
            for (int m = s; m <= 18; m++) {
                for (int e = m; e <= 18; e++) {
                    ck_assert_msg(allele_offset(&a, (Range) { s, e }) == allele_offset(&a, (Range) { s, m }) + allele_offset(&a, (Range) { m, e }),
                        "%s>%s over [%d, %d) split at %d", grid_variants[i][0], grid_variants[i][1], s, e, m);
                }
            }
        }
    }
} END_TEST

START_TEST(check_allele_offset__a_window_around_the_variant_gives_its_whole_offset) {
    for (int i = 0; i < N_GRID_VARIANTS; i++) {
        const Alleles a = Alleles_from_variant(mock_bcf1_t(10, (char *) grid_variants[i][0], (char *) grid_variants[i][1]));
        const Range footprint = variant_footprint(&a);
        const int whole = allele_offset(&a, (Range) { -1000, 1000 });

        for (int before = 0; before <= 5; before++) {
            for (int after = 0; after <= 5; after++) {
                ck_assert_msg(allele_offset(&a, (Range) { footprint.start - before, footprint.end + after }) == whole,
                    "%s>%s with %d before and %d after", grid_variants[i][0], grid_variants[i][1], before, after);
            }
        }
    }
} END_TEST

START_TEST(check_allele_offset__an_empty_or_inverted_window_holds_nothing) {
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "A"), (Range) { 11, 11 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "ATG", "A"), (Range) { 12, 11 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "A", "ACC"), (Range) { 10, 9 }), 0);
    ck_assert_int_eq(offset_of(mock_bcf1_t(10, "A", "ACC"), (Range) { 10, 10 }), 0);
} END_TEST

// The one case where the view is not just "minus the anchor": a longer shared prefix. The ALT is attached to the base
// before the first differing one (here 11), where the old code attached it to the first base of the record (10) and
// counted the shared base 11 as removed when the window started there.
START_TEST(check_allele_offset__longer_prefix_attaches_to_the_base_before_the_change) {
    bcf1_t *v = mock_bcf1_t(10, "AT", "ATCC");
    const Alleles a = Alleles_from_variant(v);

    ck_assert_int_eq(allele_offset(&a, (Range) { 11, 20 }), -2);
    ck_assert_int_eq(old_allele_offset(v, (Range) { 11, 20 }), 1);

    ck_assert_int_eq(allele_offset(&a, (Range) { 0, 20 }), old_allele_offset(v, (Range) { 0, 20 }));
} END_TEST

START_TEST(check_variant_footprint__insertion_is_the_base_before_it) {
    bcf1_t *v = mock_bcf1_t(10, "A", "ATT");
    const Alleles a = Alleles_from_variant(v);

    const Range footprint = variant_footprint(&a);
    ck_assert_int_eq(footprint.start, 10);
    ck_assert_int_eq(footprint.end, 11);

    ck_assert(variant_overlaps(&a, (Range) { 10, 11 }));
    ck_assert(!variant_overlaps(&a, (Range) { 11, 20 })); // An insertion at 11 counts for [s, e) when s < 11 <= e
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
    destroy_records(&rs);
} END_TEST

START_TEST(check_haplotype_offset__touching_indels) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "A", "AT"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(1, "ATT", "A"), NULL));

    int offset = haplotype_offset(&rs, (Range) {0, 10} );

    ck_assert_int_eq(offset, 1);
    destroy_records(&rs);
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

    int offset = prediction_offset(&kv_A(rs, 1)->alleles,&rs, 50);

    ck_assert_int_eq(offset, 5);
    destroy_records(&rs);
} END_TEST

START_TEST(check_prediction_offset__long_deletion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(59, "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "G"), NULL));

    int offset = prediction_offset(&kv_A(rs, 1)->alleles,&rs, 50);

    ck_assert_int_eq(offset, 41);
    destroy_records(&rs);
} END_TEST

START_TEST(check_prediction_offset__long_insertion) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(46, "A", "T"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(59, "G", "GACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), NULL));

    int offset = prediction_offset(&kv_A(rs, 1)->alleles,&rs, 50);

    ck_assert_int_eq(offset, -96);
    destroy_records(&rs);
} END_TEST

START_TEST(check_prediction_offset__depends_on_prediction_radius) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(80, "AT", "A"), NULL)); // 30 bases from the variant at 50

    // Window {0, 101} contains the deletion at 80, which shortens the sequence by 1
    ck_assert_int_eq(prediction_offset(&kv_A(rs, 0)->alleles,&rs, 50), 1);
    ck_assert_int_eq(prediction_offset(&kv_A(rs, 0)->alleles,&rs, 100), 1);

    // Window {30, 71} does not
    ck_assert_int_eq(prediction_offset(&kv_A(rs, 0)->alleles,&rs, 20), 0);
    destroy_records(&rs);
} END_TEST


// The window is centred on the base an indel is attached to (its VCF position), not on its first changed base: centred on
// 50, [0, 101) holds the insertion attached to 0; centred on 51 it would not
START_TEST(check_prediction_offset__indel_window_is_centred_on_the_attached_base) {
    Records rs;
    kv_init(rs);

    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "N", "Nx"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(50, "N", "Nxy"), NULL));

    ck_assert_int_eq(prediction_offset(&kv_A(rs, 1)->alleles, &rs, 50), -3);
    destroy_records(&rs);
} END_TEST

/*
* Region to haplotype, a plain region (the border cases are with haplotype_shift)
*/
START_TEST(check_region_to_haplotype__without_variants_is_unchanged) {
    Records rs;
    kv_init(rs);

    Range r = region_to_haplotype(&rs, (Range) { 0, 1000 }, (Range) { 200, 300 }, false);

    ck_assert_int_eq(r.start, 200);
    ck_assert_int_eq(r.end, 300);
} END_TEST

START_TEST(check_region_to_haplotype__deletion_before_range_shifts_both_ends) {
    Records rs;
    kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(80, "AT", "A"), NULL)); // Removes 1 base

    Range r = region_to_haplotype(&rs, (Range) { 0, 1000 }, (Range) { 200, 300 }, false);

    ck_assert_int_eq(r.start, 199);
    ck_assert_int_eq(r.end, 299);
    destroy_records(&rs);
} END_TEST

START_TEST(check_region_to_haplotype__insertion_before_range_shifts_both_ends) {
    Records rs;
    kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(80, "A", "AGG"), NULL)); // Adds 2 bases

    Range r = region_to_haplotype(&rs, (Range) { 0, 1000 }, (Range) { 200, 300 }, false);

    ck_assert_int_eq(r.start, 202);
    ck_assert_int_eq(r.end, 302);
    destroy_records(&rs);
} END_TEST

START_TEST(check_region_to_haplotype__variant_after_range_is_ignored) {
    Records rs;
    kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(500, "AT", "A"), NULL));

    Range r = region_to_haplotype(&rs, (Range) { 0, 1000 }, (Range) { 200, 300 }, false);

    ck_assert_int_eq(r.start, 200);
    ck_assert_int_eq(r.end, 300);
    destroy_records(&rs);
} END_TEST

START_TEST(check_region_to_haplotype__variant_inside_range_only_shifts_the_end) {
    Records rs;
    kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(250, "AT", "A"), NULL));

    Range r = region_to_haplotype(&rs, (Range) { 0, 1000 }, (Range) { 200, 300 }, false);

    ck_assert_int_eq(r.start, 200);
    ck_assert_int_eq(r.end, 299);
    destroy_records(&rs);
} END_TEST

/*
* reference_range_to_haplotype: a stretch of reference bases of the gene [100, 200) in haplotype coordinates. N>Nxy at p inserts two
* bases between p and p + 1; NN>N at p deletes p + 1.
*/
static Range reference_range_in_haplotype(const int pos, const char *ref, const char *alt, const Range range, const uint32_t strand) {
    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(pos, ref, alt), NULL));

    const Range hap = reference_range_to_haplotype(&rs, region_100_200, range, strand);

    destroy_records(&rs);
    return hap;
}

#define CHECK_HAP_RANGE_ON(strand, pos, ref, alt, range_start, range_end, expected_start, expected_end) do { \
    const Range got = reference_range_in_haplotype((pos), (ref), (alt), (Range) { (range_start), (range_end) }, (strand)); \
    ck_assert_msg(got.start == (expected_start) && got.end == (expected_end), "strand %d, %s>%s at %d, [%d, %d): got [%ld, %ld), expected [%d, %d)", \
                  (int) (strand), (ref), (alt), (int) (pos), (int) (range_start), (int) (range_end), (long) got.start, (long) got.end, \
                  (int) (expected_start), (int) (expected_end)); \
} while (0)

// The forward strand, where an insertion is attached to the base before it
#define CHECK_HAP_RANGE(pos, ref, alt, range_start, range_end, expected_start, expected_end) \
    CHECK_HAP_RANGE_ON(STRAND_FWD, pos, ref, alt, range_start, range_end, expected_start, expected_end)

// The reverse strand, where it is attached to the base after it: in the direction of the transcript that is the one before
#define CHECK_HAP_RANGE_REV(pos, ref, alt, range_start, range_end, expected_start, expected_end) \
    CHECK_HAP_RANGE_ON(STRAND_REV, pos, ref, alt, range_start, range_end, expected_start, expected_end)

START_TEST(check_reference_range_to_haplotype__without_records_is_unchanged) {
    Records rs; kv_init(rs);

    const uint32_t strands[3] = { STRAND_FWD, STRAND_REV, STRAND_UNK };
    for (int si = 0; si < 3; si++) {
        Range hap = reference_range_to_haplotype(&rs, region_100_200, (Range) { 150, 160 }, strands[si]);
        ck_assert_int_eq(hap.start, 150);
        ck_assert_int_eq(hap.end, 160);

        hap = reference_range_to_haplotype(&rs, region_100_200, (Range) { 150, 150 }, strands[si]);
        ck_assert_int_eq(hap.start, 150);
        ck_assert_int_eq(hap.end, 150);
    }

    kv_destroy(rs);
} END_TEST

START_TEST(check_reference_range_to_haplotype__substitutions_move_nothing) {
    CHECK_HAP_RANGE(150, "N", "x", 150, 160, 150, 160);
    CHECK_HAP_RANGE(159, "N", "x", 150, 160, 150, 160);
    CHECK_HAP_RANGE(155, "NN", "xy", 150, 160, 150, 160);
    CHECK_HAP_RANGE(149, "NN", "xy", 150, 160, 150, 160);   // An MNV across the start
} END_TEST

START_TEST(check_reference_range_to_haplotype__a_deletion_moves_what_lies_after_its_deleted_bases) {
    CHECK_HAP_RANGE(120, "NN", "N", 150, 160, 149, 159);   // Deletes 121, before
    CHECK_HAP_RANGE(152, "NN", "N", 150, 160, 150, 159);   // Deletes 153, inside
    CHECK_HAP_RANGE(170, "NN", "N", 150, 160, 150, 160);   // Deletes 171, after
    CHECK_HAP_RANGE(159, "NN", "N", 150, 160, 150, 160);   // Deletes 160, the first base after the range: the range is untouched
} END_TEST

START_TEST(check_reference_range_to_haplotype__a_deletion_across_a_border_removes_only_the_deleted_bases) {
    CHECK_HAP_RANGE(148, "NNNN", "N", 150, 160, 149, 157);   // Deletes 149 to 151: 150 and 151 were in the range
    CHECK_HAP_RANGE(157, "NNNN", "N", 150, 160, 150, 158);   // Deletes 158 to 160: 158 and 159 were in the range
    CHECK_HAP_RANGE(148, "NNNNNNNNNNNNNN", "N", 150, 160, 149, 149);   // Deletes 149 to 161: all of it, so it is empty
} END_TEST

START_TEST(check_reference_range_to_haplotype__an_insertion_inside_lengthens_it) {
    CHECK_HAP_RANGE(155, "N", "Nxy", 150, 160, 150, 162);
    CHECK_HAP_RANGE(150, "N", "Nxy", 150, 160, 150, 162);   // Between 150 and 151: attached to the first base
    CHECK_HAP_RANGE(158, "N", "Nxyzw", 150, 160, 150, 164);
} END_TEST

// An insertion at a border belongs to the side of the base it is attached to, the one before it
START_TEST(check_reference_range_to_haplotype__an_insertion_at_the_start_is_outside_and_at_the_end_inside) {
    CHECK_HAP_RANGE(149, "N", "Nxy", 150, 160, 152, 162);   // Between 149 and 150: outside, so the range starts after it
    CHECK_HAP_RANGE(159, "N", "Nxy", 150, 160, 150, 162);   // Between 159 and 160: attached to 159, inside
    CHECK_HAP_RANGE(160, "N", "Nxy", 150, 160, 150, 160);   // Between 160 and 161: outside
    CHECK_HAP_RANGE(120, "N", "Nxy", 150, 160, 152, 162);   // Before: everything moves
    CHECK_HAP_RANGE(170, "N", "Nxy", 150, 160, 150, 160);   // After
} END_TEST

START_TEST(check_reference_range_to_haplotype__an_empty_range_stays_empty) {
    CHECK_HAP_RANGE(149, "N", "Nxy", 150, 150, 152, 152);
    CHECK_HAP_RANGE(120, "NN", "N", 150, 150, 149, 149);
    CHECK_HAP_RANGE(155, "N", "Nxy", 150, 150, 150, 150);
    CHECK_HAP_RANGE(120, "NN", "N", 160, 150, 159, 159);   // And so does an inverted one
} END_TEST

START_TEST(check_reference_range_to_haplotype__at_the_gene_start) {
    CHECK_HAP_RANGE(100, "NN", "N", 100, 110, 100, 109);   // Deletes 101
    CHECK_HAP_RANGE(99, "N", "Nxy", 100, 110, 100, 110);    // An insertion before the gene start is not in the gene
    CHECK_HAP_RANGE(100, "N", "Nxy", 100, 110, 100, 112);   // Between 100 and 101
} END_TEST

START_TEST(check_reference_range_to_haplotype__records_add_up) {
    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(120, "NN", "N"), NULL));    // Deletes 121
    kv_push(Record *, rs, Record_init(mock_bcf1_t(130, "N", "Nxy"), NULL));   // Inserts 2 after 130
    kv_push(Record *, rs, Record_init(mock_bcf1_t(155, "N", "Nxyz"), NULL));  // Inserts 3 after 155

    // None of the insertions is at a border, so the strand does not matter
    const uint32_t strands[2] = { STRAND_FWD, STRAND_REV };
    for (int si = 0; si < 2; si++) {
        const Range hap = reference_range_to_haplotype(&rs, region_100_200, (Range) { 150, 160 }, strands[si]);
        ck_assert_int_eq(hap.start, 151);   // One deleted base before it and two inserted ones: 150 - 1 + 2
        ck_assert_int_eq(hap.end, 164);     // One deleted, 2 + 3 inserted before 160: 160 - 1 + 5
    }

    destroy_records(&rs);
} END_TEST

// An unknown strand is read as forward
START_TEST(check_reference_range_to_haplotype__an_unknown_strand_is_read_as_forward) {
    CHECK_HAP_RANGE_ON(STRAND_UNK, 149, "N", "Nxy", 150, 160, 152, 162);
    CHECK_HAP_RANGE_ON(STRAND_UNK, 159, "N", "Nxy", 150, 160, 150, 162);
} END_TEST

// On the reverse strand an insertion is attached to the base after it, which is the one before it in the direction of the
// transcript: at the start border (the base after it is the first of the range) it is inside, and at the end border it is outside
START_TEST(check_reference_range_to_haplotype__reverse_an_insertion_at_the_start_is_inside_and_at_the_end_outside) {
    CHECK_HAP_RANGE_REV(149, "N", "Nxy", 150, 160, 150, 162);   // Between 149 and 150: attached to 150, inside
    CHECK_HAP_RANGE_REV(159, "N", "Nxy", 150, 160, 150, 160);   // Between 159 and 160: attached to 160, outside
    CHECK_HAP_RANGE_REV(160, "N", "Nxy", 150, 160, 150, 160);   // Between 160 and 161: outside
    CHECK_HAP_RANGE_REV(148, "N", "Nxy", 150, 160, 152, 162);   // Between 148 and 149: before the range, everything moves
    CHECK_HAP_RANGE_REV(155, "N", "Nxy", 150, 160, 150, 162);   // Inside
    CHECK_HAP_RANGE_REV(150, "N", "Nxy", 150, 160, 150, 162);   // Between 150 and 151: inside
    CHECK_HAP_RANGE_REV(120, "N", "Nxy", 150, 160, 152, 162);   // Before
    CHECK_HAP_RANGE_REV(170, "N", "Nxy", 150, 160, 150, 160);   // After
} END_TEST

// Everything but a border insertion is the same on both strands
START_TEST(check_reference_range_to_haplotype__reverse_deletions_and_substitutions_are_as_on_the_forward_strand) {
    CHECK_HAP_RANGE_REV(150, "N", "x", 150, 160, 150, 160);
    CHECK_HAP_RANGE_REV(155, "NN", "xy", 150, 160, 150, 160);
    CHECK_HAP_RANGE_REV(120, "NN", "N", 150, 160, 149, 159);
    CHECK_HAP_RANGE_REV(152, "NN", "N", 150, 160, 150, 159);
    CHECK_HAP_RANGE_REV(159, "NN", "N", 150, 160, 150, 160);
    CHECK_HAP_RANGE_REV(148, "NNNN", "N", 150, 160, 149, 157);
    CHECK_HAP_RANGE_REV(157, "NNNN", "N", 150, 160, 150, 158);
    CHECK_HAP_RANGE_REV(148, "NNNNNNNNNNNNNN", "N", 150, 160, 149, 149);
} END_TEST

START_TEST(check_reference_range_to_haplotype__reverse_an_empty_range_stays_empty) {
    CHECK_HAP_RANGE_REV(149, "N", "Nxy", 150, 150, 150, 150);   // The insertion is attached to the base after it, 150: not before an empty start
    CHECK_HAP_RANGE_REV(148, "N", "Nxy", 150, 150, 152, 152);
    CHECK_HAP_RANGE_REV(120, "NN", "N", 160, 150, 159, 159);
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
    free(ref_preds);
    destroy_records(&rs);
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
    free(ref_preds);
    destroy_records(&rs);
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

    // Attached to 0, before the window [1, 102): the deleted bases 1 and 2 are both deleted
    ck_assert_float_eq(ref_preds[0], 0);
    ck_assert_float_eq(ref_preds[1], 0);
    ck_assert_float_eq(ref_preds[2], 0);
    ck_assert_float_eq(ref_preds[3], 0);
    ck_assert_float_eq(ref_preds[4], 0);
    ck_assert_float_eq(ref_preds[5], 0);
    free(ref_preds);
    destroy_records(&rs);
} END_TEST

// match_predictions_to_reference as it was before the anchor-free view, kept to compare against. It leaves positions
// unwritten for complex indels and reads an uninitialised range for a deletion anchored two or more bases before the window,
// so it is only compared where it was well defined.
static float *old_match_predictions_to_reference(const float *predictions, const Range pred_range, const Records *rs, const size_t num_reference_predictions) {
    float *corrected_predictions = malloc(num_reference_predictions * sizeof(float));

    int idx = 0;
    for (int l = 0; l < rs->n; l++) {
        Record *r = kv_A(*rs, l);

        Range ref;
        if (intersect_ranges(variant_to_range(r->v, 0), pred_range, &ref) == 0) continue;
        shift_range(&ref, -pred_range.start);
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt = variant_to_range(r->v, 1);
        Range intersected_alt;
        int alt_segment = 0;
        if (intersect_ranges(alt, pred_range, &intersected_alt) > 0) {
            alt_segment = range_len(alt);
        }

        int skipped_alt_bases = (intersected_alt.start - alt.start);
        if (alt_segment > 0) alt_segment -= skipped_alt_bases;
        Range ref_segment = scale_range(ref, NUM_SCORES);

        int variant_offset = haplotype_offset(rs, (Range) { pred_range.start, pred_range.start + ref.start }) * NUM_SCORES;

        int segment_len = (ref_segment.start - idx);
        if (segment_len > 0) memcpy(corrected_predictions + idx, predictions + idx - variant_offset, segment_len * sizeof(float));

        int offset = ref_len - alt_segment;
        float *corrected_start = corrected_predictions + ref_segment.start;
        if (offset == 0) {
            memcpy(corrected_start, predictions + ref_segment.start - variant_offset, range_len(ref_segment) * sizeof(float));
        } else if (offset < 0) {
            float best_scores[NUM_SCORES] = {0};
            get_best_insertion_score(predictions + ref_segment.start - variant_offset, alt_segment * NUM_SCORES, best_scores);
            memcpy(corrected_start, best_scores, NUM_SCORES * sizeof(float));
        } else {
            memcpy(corrected_start, predictions + ref_segment.start - variant_offset, NUM_SCORES * sizeof(float));
            memset(corrected_start + NUM_SCORES, 0, (offset - (skipped_alt_bases > 0)) * NUM_SCORES * sizeof(float));
        }

        if (ref_segment.end > idx) idx = ref_segment.end;
    }

    int variant_offset = haplotype_offset(rs, pred_range) * NUM_SCORES;
    memcpy(corrected_predictions + idx, predictions + idx - variant_offset, (num_reference_predictions - idx) * sizeof(float));

    return corrected_predictions;
}

// Distinct values, so a prediction taken from the wrong position shows; long enough for any window below
#define GRID_PREDICTIONS 1200
static float *grid_predictions(void) {
    float *predictions = malloc(GRID_PREDICTIONS * sizeof(float));
    for (int i = 0; i < GRID_PREDICTIONS; i++) predictions[i] = (float) (i + 1) / GRID_PREDICTIONS;
    return predictions;
}

#define N_PRED_VARIANTS 6
static const char *pred_variants[N_PRED_VARIANTS][2] = {
    { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "N", "Nxyz" }, { "NN", "N" }, { "NNNN", "N" }
};

// Every normalised SNV, MNV, insertion and deletion around both window edges, where the old version was well defined: not
// a deletion attached outside the window (pinned separately below) and not an MNV that runs past the window end.
START_TEST(check_match_predictions_to_reference__single_variants_give_the_old_predictions) {
    const Range window = range_around(150, 50); // [100, 201)
    const size_t n = pred_array_len(50);
    float *predictions = grid_predictions();

    for (int i = 0; i < N_PRED_VARIANTS; i++) {
        for (int pos = 95; pos <= 205; pos++) {
            Records rs; kv_init(rs);
            kv_push(Record *, rs, Record_init(mock_bcf1_t(pos, (char *) pred_variants[i][0], (char *) pred_variants[i][1]), NULL));
            const Alleles *a = &kv_A(rs, 0)->alleles;

            const bool deletion_attached_outside = a->ref_len > a->alt_len && a->pos - 1 < window.start;
            const bool mnv_past_the_end = a->prefix_len == 0 && a->ref_len > 1 && a->pos + a->ref_len > window.end;
            if (!deletion_attached_outside && !mnv_past_the_end) {
                float *new = match_predictions_to_reference(predictions, window, &rs, n);
                float *old = old_match_predictions_to_reference(predictions, window, &rs, n);

                for (size_t k = 0; k < n; k++) {
                    ck_assert_msg(new[k] == old[k], "%s>%s at %d, float %zu: %f, old %f", pred_variants[i][0], pred_variants[i][1], pos, k, new[k], old[k]);
                }

                free(new); free(old);
            }

            Record_destroy(kv_A(rs, 0));
            kv_destroy(rs);
        }
    }

    free(predictions);
} END_TEST

// A deletion attached to the base before the window: the attachment base is not in the window, so no position stands in for
// the junction and every deleted base in the window gets 0. (The old version gave the first one the prediction of the first
// haplotype base.)
START_TEST(check_match_predictions_to_reference__deletion_attached_before_the_window) {
    const Range window = range_around(60, 50); // [10, 111)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(9, "NNN", "N"), NULL)); // Attached to 9, removes 10 and 11

    float *new = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(new[k], 0); // Position 10
        ck_assert_float_eq(new[NUM_SCORES + k], 0); // Position 11
        ck_assert_float_eq(new[2 * NUM_SCORES + k], predictions[k]); // Position 12 is the first haplotype base
    }

    free(new); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// Further before the window the old version read an uninitialised range; the rule is the same as just before it
START_TEST(check_match_predictions_to_reference__deletion_attached_further_before_the_window) {
    const Range window = range_around(60, 50); // [10, 111)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(7, "NNNNN", "N"), NULL)); // Attached to 7, removes 8 to 11

    float *new = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(new[k], 0); // Position 10
        ck_assert_float_eq(new[NUM_SCORES + k], 0); // Position 11
        ck_assert_float_eq(new[2 * NUM_SCORES + k], predictions[k]); // Position 12
    }

    free(new); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// Every attachment point before the window, for a deletion of 1 to 3 bases: the deleted bases in the window get 0 and the
// bases after the deletion take the predictions of the haplotype, shifted by the deleted bases in the window
START_TEST(check_match_predictions_to_reference__every_deleted_base_in_the_window_gets_0) {
    const Range window = range_around(60, 50); // [10, 111)
    float *predictions = grid_predictions();
    static const char *refs[3] = { "NN", "NNN", "NNNN" };

    for (int deleted = 1; deleted <= 3; deleted++) {
        for (int attached = 4; attached < window.start; attached++) {
            Records rs; kv_init(rs);
            kv_push(Record *, rs, Record_init(mock_bcf1_t(attached, (char *) refs[deleted - 1], "N"), NULL));

            const int64_t first_deleted = attached + 1, after = attached + 1 + deleted;
            if (after <= window.start) { // Nothing deleted in the window
                Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
                continue;
            }

            float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));
            const int64_t in_window = after - (first_deleted < window.start ? window.start : first_deleted);

            for (int64_t x = window.start; x < after; x++) {
                for (int k = 0; k < NUM_SCORES; k++) {
                    ck_assert_msg(corrected[(x - window.start) * NUM_SCORES + k] == 0, "%d deleted after %d, position %ld", deleted, attached, x);
                }
            }
            for (int64_t x = after; x < window.start + 10; x++) {
                for (int k = 0; k < NUM_SCORES; k++) {
                    ck_assert_msg(corrected[(x - window.start) * NUM_SCORES + k] == predictions[(x - window.start - in_window) * NUM_SCORES + k], "%d deleted after %d, position %ld", deleted, attached, x);
                }
            }

            free(corrected);
            Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
        }
    }

    free(predictions);
} END_TEST

// The prediction array of a haplotype has exactly the window width minus the net deletions in the window (the model returns the
// sequence length minus its context, and the caller gives it model_width - haplotype_offset bases). Allocated to that size,
// so a read past it is an out-of-bounds read under AddressSanitizer. Values grow with the index, so the best of several
// positions is the last one.
static float *exact_predictions(const Records *rs, const Range window, int *positions) {
    *positions = pred_width(50) - haplotype_offset(rs, window);
    float *predictions = malloc(*positions * NUM_SCORES * sizeof(float));
    for (int i = 0; i < *positions * NUM_SCORES; i++) predictions[i] = (float) (i + 1) / GRID_PREDICTIONS;
    return predictions;
}

// The insertion fold reads the junction and the inserted bases: 1 + alt_len - ref_len positions from the junction. For an
// insertion attached to the last base of the window these are the last positions of the array, and no more.
START_TEST(check_match_predictions_to_reference__insertion_fold_at_the_last_base_stays_in_the_array) {
    const Range window = range_around(150, 50); // [100, 201)

    for (int inserted = 1; inserted <= 3; inserted++) {
        static const char *alts[3] = { "Nx", "Nxy", "Nxyz" };
        Records rs; kv_init(rs);
        kv_push(Record *, rs, Record_init(mock_bcf1_t(200, "N", (char *) alts[inserted - 1]), NULL)); // Attached to 200

        int positions;
        float *predictions = exact_predictions(&rs, window, &positions);
        ck_assert_int_eq(positions, 101 + inserted);

        float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

        for (int k = 0; k < NUM_SCORES; k++) {
            ck_assert_float_eq(corrected[99 * NUM_SCORES + k], predictions[99 * NUM_SCORES + k]);
            ck_assert_float_eq(corrected[100 * NUM_SCORES + k], predictions[(100 + inserted) * NUM_SCORES + k]); // The last inserted base
        }

        free(corrected); free(predictions);
        Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
    }
} END_TEST

// A deletion before the insertion moves the junction one position up in the array and shortens the array by one
START_TEST(check_match_predictions_to_reference__insertion_fold_after_a_deletion_stays_in_the_array) {
    const Range window = range_around(150, 50); // [100, 201)

    for (int inserted = 1; inserted <= 3; inserted++) {
        static const char *alts[3] = { "Nx", "Nxy", "Nxyz" };
        Records rs; kv_init(rs);
        kv_push(Record *, rs, Record_init(mock_bcf1_t(150, "NN", "N"), NULL));   // Deletes 151
        kv_push(Record *, rs, Record_init(mock_bcf1_t(200, "N", (char *) alts[inserted - 1]), NULL)); // Attached to 200

        int positions;
        float *predictions = exact_predictions(&rs, window, &positions);
        ck_assert_int_eq(positions, 100 + inserted);

        float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

        for (int k = 0; k < NUM_SCORES; k++) {
            ck_assert_float_eq(corrected[51 * NUM_SCORES + k], 0); // The deleted base
            ck_assert_float_eq(corrected[100 * NUM_SCORES + k], predictions[(99 + inserted) * NUM_SCORES + k]);
        }

        free(corrected); free(predictions);
        Record_destroy(kv_A(rs, 0)); Record_destroy(kv_A(rs, 1)); kv_destroy(rs);
    }
} END_TEST

// A deletion of the last base of the window after the insertion: the array is one shorter again
START_TEST(check_match_predictions_to_reference__insertion_fold_before_a_deletion_at_the_end_stays_in_the_array) {
    const Range window = range_around(150, 50); // [100, 201)

    for (int inserted = 1; inserted <= 3; inserted++) {
        static const char *alts[3] = { "Nx", "Nxy", "Nxyz" };
        Records rs; kv_init(rs);
        kv_push(Record *, rs, Record_init(mock_bcf1_t(198, "N", (char *) alts[inserted - 1]), NULL)); // Attached to 198
        kv_push(Record *, rs, Record_init(mock_bcf1_t(199, "NN", "N"), NULL));   // Attached to 199, deletes 200

        int positions;
        float *predictions = exact_predictions(&rs, window, &positions);
        ck_assert_int_eq(positions, 100 + inserted);

        float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

        for (int k = 0; k < NUM_SCORES; k++) {
            ck_assert_float_eq(corrected[98 * NUM_SCORES + k], predictions[(98 + inserted) * NUM_SCORES + k]); // The fold
            ck_assert_float_eq(corrected[99 * NUM_SCORES + k], predictions[(99 + inserted) * NUM_SCORES + k]); // Attachment base of the deletion
            ck_assert_float_eq(corrected[100 * NUM_SCORES + k], 0); // The deleted base
        }

        free(corrected); free(predictions);
        Record_destroy(kv_A(rs, 0)); Record_destroy(kv_A(rs, 1)); kv_destroy(rs);
    }
} END_TEST

// Replaces the removed end-of-window test (Step 2): a deletion attached to the last base but one of the window, removing
// the last base and one beyond. The attachment base keeps its prediction, the removed last base gets 0.
START_TEST(check_match_predictions_to_reference__deletion_overlaps_bounds_end) {
    const Range window = range_around(50, 50); // [0, 101)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(99, "NNN", "N"), NULL)); // Attached to 99, removes 100 and 101

    float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(corrected[98 * NUM_SCORES + k], predictions[98 * NUM_SCORES + k]);
        ck_assert_float_eq(corrected[99 * NUM_SCORES + k], predictions[99 * NUM_SCORES + k]);
        ck_assert_float_eq(corrected[100 * NUM_SCORES + k], 0);
    }

    free(corrected); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// An insertion attached to the last base of the window: that base takes the best of itself and the inserted bases, which
// lie past the window end in the haplotype predictions
START_TEST(check_match_predictions_to_reference__insertion_overlaps_bounds_end) {
    const Range window = range_around(50, 50); // [0, 101)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(100, "N", "Nxy"), NULL)); // Attached to 100, inserts two bases

    float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(corrected[99 * NUM_SCORES + k], predictions[99 * NUM_SCORES + k]);
        ck_assert_float_eq(corrected[100 * NUM_SCORES + k], predictions[102 * NUM_SCORES + k]); // The values grow with the index
    }

    free(corrected); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// An MNV running past the window end is copied position by position; the old version took the best over the whole ALT,
// including the prediction past the window end
START_TEST(check_match_predictions_to_reference__mnv_past_the_window_end_is_copied) {
    const Range window = range_around(50, 50); // [0, 101)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(100, "NN", "xy"), NULL));

    float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) ck_assert_float_eq(corrected[100 * NUM_SCORES + k], predictions[100 * NUM_SCORES + k]);

    free(corrected); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// A complex indel is paired base by base from the left: the paired bases keep the prediction of the base that replaces them, the
// REF bases left over are removed and get 0. NNN>Nx at 40: x replaces the base at 41, and the base at 42 is deleted. (The old
// version gave 41 and 42 both 0, and left the third REF position of ATG>AC unwritten.)
START_TEST(check_match_predictions_to_reference__complex_deletion_keeps_the_replaced_base) {
    const Range window = range_around(50, 50); // [0, 101)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(40, "NNN", "Nx"), NULL)); // Attached to 40, x replaces 41, 42 is deleted

    float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(corrected[40 * NUM_SCORES + k], predictions[40 * NUM_SCORES + k]); // The attached base
        ck_assert_float_eq(corrected[41 * NUM_SCORES + k], predictions[41 * NUM_SCORES + k]); // Replaced by x: its prediction
        ck_assert_float_eq(corrected[42 * NUM_SCORES + k], 0);                                  // Deleted
        ck_assert_float_eq(corrected[43 * NUM_SCORES + k], predictions[42 * NUM_SCORES + k]); // One base shorter
    }

    free(corrected); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

// NN>Nxyz at 40: x replaces the base at 41 and y and z are inserted after it. The last paired base, 41, takes the best of itself
// and the inserted bases (the values grow with the index, so that is z), and what follows is two haplotype bases ahead.
START_TEST(check_match_predictions_to_reference__complex_insertion_gives_the_best_to_the_last_paired_base) {
    const Range window = range_around(50, 50); // [0, 101)
    float *predictions = grid_predictions();

    Records rs; kv_init(rs);
    kv_push(Record *, rs, Record_init(mock_bcf1_t(40, "NN", "Nxyz"), NULL));

    float *corrected = match_predictions_to_reference(predictions, window, &rs, pred_array_len(50));

    for (int k = 0; k < NUM_SCORES; k++) {
        ck_assert_float_eq(corrected[40 * NUM_SCORES + k], predictions[40 * NUM_SCORES + k]); // The attached base
        ck_assert_float_eq(corrected[41 * NUM_SCORES + k], predictions[43 * NUM_SCORES + k]); // Best of x, y and z
        ck_assert_float_eq(corrected[42 * NUM_SCORES + k], predictions[44 * NUM_SCORES + k]); // Two bases longer
    }

    free(corrected); free(predictions);
    Record_destroy(kv_A(rs, 0)); kv_destroy(rs);
} END_TEST

Suite * utils_suite(void) {
    Suite *s = suite_create("Utils");

    TCase *tc_mocks = mocked_tcase("mocks");
    suite_add_tcase(s, tc_mocks);

    tcase_add_test(tc_mocks, check_mocks_free_all__frees_what_mock_bcf1_t_made);
    tcase_add_test(tc_mocks, check_mocks_free_all__keeps_a_record_made_from_a_mock);

    TCase *tc_rng = mocked_tcase("range_helpers");
    suite_add_tcase(s, tc_rng);

    tcase_add_test(tc_rng, check_range_around__default_distance);
    tcase_add_test(tc_rng, check_range_around__zero_offset);
    tcase_add_test(tc_rng, check_range_around__near_contig_start_is_not_clamped);
    tcase_add_test(tc_rng, check_range_before__default_distance);
    tcase_add_test(tc_rng, check_range_before__ends_at_pos_and_starts_like_range_around);
    tcase_add_test(tc_rng, check_range_before__near_contig_start_is_not_clamped);

    TCase *tc_ww = mocked_tcase("window_widths");
    suite_add_tcase(s, tc_ww);

    tcase_add_test(tc_ww, check_window_widths__default_distance);
    tcase_add_test(tc_ww, check_window_widths__other_distances);
    tcase_add_test(tc_ww, check_window_widths__consistent_with_range_around);
    tcase_add_test(tc_ww, check_max_prediction_radius__derived_sizes_fit_in_int);
    tcase_add_test(tc_ww, check_max_prediction_radius__is_tight);

    TCase *tc_vfw = mocked_tcase("variant_fits_window");
    suite_add_tcase(s, tc_vfw);

    tcase_add_test(tc_vfw, check_variant_fits_window__snv);
    tcase_add_test(tc_vfw, check_variant_fits_window__alt_length_at_limit);
    tcase_add_test(tc_vfw, check_variant_fits_window__ref_length_at_limit);
    tcase_add_test(tc_vfw, check_variant_fits_window__insertion_does_not_count_anchor);
    tcase_add_test(tc_vfw, check_variant_fits_window__deletion_does_not_count_anchor);
    tcase_add_test(tc_vfw, check_variant_fits_window__depends_on_prediction_radius);

    TCase *tc_vir = mocked_tcase("variant_in_region");
    suite_add_tcase(s, tc_vir);

    tcase_add_test(tc_vir, check_variant_in_region__substitutions);
    tcase_add_test(tc_vir, check_variant_in_region__deletions_count_by_the_deleted_bases);
    tcase_add_test(tc_vir, check_variant_in_region__complex_indels_count_by_the_replaced_bases);
    tcase_add_test(tc_vir, check_variant_in_region__an_insertion_needs_both_neighbours_in_the_region);
    tcase_add_test(tc_vir, check_variant_in_region__differs_from_the_footprint_only_by_the_anchor);

    TCase *tc_hs = mocked_tcase("haplotype_shift");
    suite_add_tcase(s, tc_hs);

    tcase_add_test(tc_hs, check_haplotype_shift__substitutions_move_nothing);
    tcase_add_test(tc_hs, check_haplotype_shift__a_deletion_moves_what_lies_after_its_deleted_bases);
    tcase_add_test(tc_hs, check_haplotype_shift__an_insertion_at_b_counts_only_when_asked);
    tcase_add_test(tc_hs, check_haplotype_shift__complex_indels_are_a_substitution_and_a_leftover);
    tcase_add_test(tc_hs, check_haplotype_shift__the_records_add_up);
    tcase_add_test(tc_hs, check_region_to_haplotype__an_insertion_at_a_border_goes_to_the_side_without_the_feature);
    tcase_add_test(tc_hs, check_region_to_haplotype__an_empty_region_stays_empty);
    tcase_add_test(tc_hs, check_region_to_haplotype__a_substitution_across_a_border_changes_nothing);

    TCase *tc_rrh = mocked_tcase("reference_range_to_haplotype");
    suite_add_tcase(s, tc_rrh);

    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__without_records_is_unchanged);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__substitutions_move_nothing);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__a_deletion_moves_what_lies_after_its_deleted_bases);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__a_deletion_across_a_border_removes_only_the_deleted_bases);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__an_insertion_inside_lengthens_it);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__an_insertion_at_the_start_is_outside_and_at_the_end_inside);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__an_empty_range_stays_empty);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__at_the_gene_start);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__records_add_up);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__an_unknown_strand_is_read_as_forward);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__reverse_an_insertion_at_the_start_is_inside_and_at_the_end_outside);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__reverse_deletions_and_substitutions_are_as_on_the_forward_strand);
    tcase_add_test(tc_rrh, check_reference_range_to_haplotype__reverse_an_empty_range_stays_empty);

    TCase *tc_vp = mocked_tcase("variant_position");
    suite_add_tcase(s, tc_vp);

    tcase_add_test(tc_vp, check_variant_position__substitutions_count_at_their_first_base);
    tcase_add_test(tc_vp, check_variant_position__a_deletion_counts_at_its_first_deleted_base);
    tcase_add_test(tc_vp, check_variant_position__an_insertion_counts_at_the_base_it_is_attached_to);
    tcase_add_test(tc_vp, check_variant_position__a_complex_indel_counts_at_its_first_replaced_base);
    tcase_add_test(tc_vp, check_variant_position__window_centre_is_still_the_attached_base);

    TCase *tc_vtr = mocked_tcase("variant_to_range");
    suite_add_tcase(s, tc_vtr);

    tcase_add_test(tc_vtr, check_variant_to_range__snv);
    tcase_add_test(tc_vtr, check_variant_to_range__insertion);
    tcase_add_test(tc_vtr, check_variant_to_range__deletion);

    TCase *tc_ao = mocked_tcase("allele_offset");
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

    TCase *tc_view = mocked_tcase("anchor_free_view");
    suite_add_tcase(s, tc_view);

    tcase_add_test(tc_view, check_variant_footprint__is_the_old_ref_range_for_normalised_variants);
    tcase_add_test(tc_view, check_variant_overlaps__gives_the_old_answer_for_every_window);
    tcase_add_test(tc_view, check_allele_offset__gives_the_old_answer_for_every_window);
    tcase_add_test(tc_view, check_variant_shift__substitutions_move_nothing);
    tcase_add_test(tc_view, check_variant_shift__deleted_bases_move_what_is_after_them);
    tcase_add_test(tc_view, check_variant_shift__an_insertion_at_b_counts_only_when_asked);
    tcase_add_test(tc_view, check_variant_shift__complex_indels_are_a_substitution_and_a_leftover);
    tcase_add_test(tc_view, check_variant_shift__the_difference_over_a_window_is_the_offset);
    tcase_add_test(tc_view, check_allele_offset__is_additive_over_adjacent_windows);
    tcase_add_test(tc_view, check_allele_offset__a_window_around_the_variant_gives_its_whole_offset);
    tcase_add_test(tc_view, check_allele_offset__a_window_edge_through_a_variant);
    tcase_add_test(tc_view, check_allele_offset__an_empty_or_inverted_window_holds_nothing);
    tcase_add_test(tc_view, check_allele_offset__longer_prefix_attaches_to_the_base_before_the_change);
    tcase_add_test(tc_view, check_variant_footprint__insertion_is_the_base_before_it);

    TCase *tc_ho = mocked_tcase("haplotype_offset");
    suite_add_tcase(s, tc_ho);

    tcase_add_test(tc_ho, check_haplotype_offset);
    tcase_add_test(tc_ho, check_haplotype_offset__touching_indels);

    TCase *tc_po = mocked_tcase("prediction_offset");
    suite_add_tcase(s, tc_po);

    tcase_add_test(tc_po, check_prediction_offset);
    tcase_add_test(tc_po, check_prediction_offset__long_deletion);
    tcase_add_test(tc_po, check_prediction_offset__long_insertion);
    tcase_add_test(tc_po, check_prediction_offset__depends_on_prediction_radius);
    tcase_add_test(tc_po, check_prediction_offset__indel_window_is_centred_on_the_attached_base);

    TCase *tc_rth = mocked_tcase("region_to_haplotype");
    suite_add_tcase(s, tc_rth);

    tcase_add_test(tc_rth, check_region_to_haplotype__without_variants_is_unchanged);
    tcase_add_test(tc_rth, check_region_to_haplotype__deletion_before_range_shifts_both_ends);
    tcase_add_test(tc_rth, check_region_to_haplotype__insertion_before_range_shifts_both_ends);
    tcase_add_test(tc_rth, check_region_to_haplotype__variant_after_range_is_ignored);
    tcase_add_test(tc_rth, check_region_to_haplotype__variant_inside_range_only_shifts_the_end);

    TCase *tc_cmptf = mocked_tcase("match_predictions_to_reference");
    suite_add_tcase(s, tc_cmptf);

    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_within_bounds);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_overlaps_bounds_start);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__single_variants_give_the_old_predictions);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_attached_before_the_window);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_attached_further_before_the_window);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__every_deleted_base_in_the_window_gets_0);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_fold_at_the_last_base_stays_in_the_array);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_fold_after_a_deletion_stays_in_the_array);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_fold_before_a_deletion_at_the_end_stays_in_the_array);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__deletion_overlaps_bounds_end);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__insertion_overlaps_bounds_end);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__mnv_past_the_window_end_is_copied);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__complex_deletion_keeps_the_replaced_base);
    tcase_add_test(tc_cmptf, check_match_predictions_to_reference__complex_insertion_gives_the_best_to_the_last_paired_base);


    return s;
}

