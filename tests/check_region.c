#include <check.h>
#include <htslib/vcf.h>

#include "mocks.h"
#include "region.h"

#include "check_region.h"

// An IRegion is half-open: the bases start to end - 1, so it can be compared with the others. The GFF and the region index give
// the first and the last base, both inclusive.

START_TEST(check_iregion_from_bases__is_half_open) {
    const IRegion r = iregion_from_bases(2, 100, 200); // The bases 100 to 200

    ck_assert_int_eq(r.rid, 2);
    ck_assert_int_eq(r.start, 100);
    ck_assert_int_eq(r.end, 201);
} END_TEST

START_TEST(check_iregion_last_base__is_the_inclusive_end_for_index_queries) {
    ck_assert_int_eq(iregion_last_base((IRegion) { 0, 100, 201 }), 200);
    ck_assert_int_eq(iregion_last_base((IRegion) { 0, 100, 101 }), 100); // One base
    ck_assert_int_eq(iregion_last_base((IRegion) { 0, 100, 100 }), 100); // Empty: not before its start
} END_TEST

START_TEST(check_iregion_intersects__shares_at_least_one_base) {
    const IRegion a = { 0, 100, 200 }; // The bases 100 to 199

    ck_assert(iregion_intersects(a, (IRegion) { 0, 199, 300 }));   // Shares 199
    ck_assert(iregion_intersects((IRegion) { 0, 199, 300 }, a));
    ck_assert(!iregion_intersects(a, (IRegion) { 0, 200, 300 }));  // Next to it
    ck_assert(!iregion_intersects((IRegion) { 0, 200, 300 }, a));
    ck_assert(iregion_intersects(a, (IRegion) { 0, 0, 101 }));     // Shares 100
    ck_assert(!iregion_intersects(a, (IRegion) { 0, 0, 100 }));
    ck_assert(!iregion_intersects(a, (IRegion) { 1, 100, 200 }));  // Another sequence
} END_TEST

// Two genes that share their last and first base are in one batch. With the inclusive gene end stored as the end of the region
// (the old code) the shared base was lost.
START_TEST(check_iregion_intersects__genes_sharing_one_base_batch_together) {
    const IRegion first = iregion_from_bases(0, 100, 200);  // The bases 100 to 200
    const IRegion second = iregion_from_bases(0, 200, 300); // The bases 200 to 300
    const IRegion apart = iregion_from_bases(0, 201, 300);  // Starts after the first ends

    ck_assert(iregion_intersects(first, second));
    ck_assert(!iregion_intersects(first, apart));

    ck_assert(!iregion_intersects((IRegion) { 0, 100, 200 }, (IRegion) { 0, 200, 300 })); // The old representation
} END_TEST

START_TEST(check_variant_to_iregion__covers_the_ref_bases) {
    bcf1_t *snv = mock_bcf1_t(10, "N", "x");
    bcf1_t *deletion = mock_bcf1_t(10, "NNN", "N");
    bcf1_t *insertion = mock_bcf1_t(10, "N", "Nxy");

    IRegion r = variant_to_iregion(snv);
    ck_assert_int_eq(r.start, 10); ck_assert_int_eq(r.end, 11);

    r = variant_to_iregion(deletion);
    ck_assert_int_eq(r.start, 10); ck_assert_int_eq(r.end, 13);

    r = variant_to_iregion(insertion); // Only the base it is attached to
    ck_assert_int_eq(r.start, 10); ck_assert_int_eq(r.end, 11);
} END_TEST

// The region index takes an inclusive end: a query for a variant reaches its last REF base and not the one after it
START_TEST(check_variant_to_iregion__the_index_query_ends_on_the_last_ref_base) {
    bcf1_t *snv = mock_bcf1_t(10, "N", "x");
    bcf1_t *deletion = mock_bcf1_t(10, "NNN", "N");

    ck_assert_int_eq(iregion_last_base(variant_to_iregion(snv)), 10);
    ck_assert_int_eq(iregion_last_base(variant_to_iregion(deletion)), 12);
} END_TEST

Suite *region_suite() {
    Suite *s = suite_create("Region");

    TCase *tc = mocked_tcase("iregion");
    tcase_add_test(tc, check_iregion_from_bases__is_half_open);
    tcase_add_test(tc, check_iregion_last_base__is_the_inclusive_end_for_index_queries);
    tcase_add_test(tc, check_iregion_intersects__shares_at_least_one_base);
    tcase_add_test(tc, check_iregion_intersects__genes_sharing_one_base_batch_together);
    tcase_add_test(tc, check_variant_to_iregion__covers_the_ref_bases);
    tcase_add_test(tc, check_variant_to_iregion__the_index_query_ends_on_the_last_ref_base);
    suite_add_tcase(s, tc);

    return s;
}
