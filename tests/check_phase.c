#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>
#include <cpliceai/range.h>

#include "klib/kvec.h"
#include "phase.h"
#include "record.h"
#include "utils.h"

#include "check_phase.h"
#include "mocks.h"

// Present on the haplotype of a phase for a genotype, with and without unphased genotypes included
static int present(const int gt0, const int gt1, const bool phased, const int phase, const bool include_unphased) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,G", gt0, gt1, phased, -1);

    int result = has_matching_phase(v, hdr, phase, include_unphased);

    bcf_destroy(v);
    bcf_hdr_destroy(hdr);

    return result;
}

START_TEST(check_has_matching_phase__homozygous_is_on_both_haplotypes) {
    for (int phase = 0; phase < 2; phase++) {
        ck_assert_int_eq(present(1, 1, true, phase, false), 1);
        ck_assert_int_eq(present(1, 1, false, phase, false), 1);
    }
}
END_TEST

START_TEST(check_has_matching_phase__homozygous_reference_also_matches) {
    // Not what the name suggests, but it is what the code does
    for (int phase = 0; phase < 2; phase++) {
        ck_assert_int_eq(present(0, 0, true, phase, false), 1);
    }
}
END_TEST

START_TEST(check_has_matching_phase__phased_heterozygous_is_on_the_haplotype_with_the_alt) {
    // 1|0 carries the alt on the first haplotype
    ck_assert_int_ne(present(1, 0, true, 0, false), 0);
    ck_assert_int_eq(present(1, 0, true, 1, false), 0);

    // 0|1 carries it on the second
    ck_assert_int_eq(present(0, 1, true, 0, false), 0);
    ck_assert_int_ne(present(0, 1, true, 1, false), 0);
}
END_TEST

START_TEST(check_has_matching_phase__unphased_heterozygous_follows_include_unphased) {
    for (int phase = 0; phase < 2; phase++) {
        ck_assert_int_eq(present(0, 1, false, phase, true), 1);
        ck_assert_int_eq(present(0, 1, false, phase, false), 0);
    }
}
END_TEST

// A record with the given genotype. The loader gives a record with coding transcripts a tscripts vector, and leaves it NULL otherwise.
static Record *make_record(bcf_hdr_t *hdr, const int64_t pos, const int gt0, const int gt1, const bool phased, const bool with_transcripts) {
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, pos, "A,G", gt0, gt1, phased, -1);
    Record *r = Record_init(v, NULL);
    bcf_destroy(v);

    if (with_transcripts) {
        r->tscripts = malloc(sizeof(Transcripts));
        kv_init(*r->tscripts);
    }

    return r;
}

static void free_records(Records *rs) {
    for (size_t i = 0; i < kv_size(*rs); i++) Record_destroy(kv_A(*rs, i));
    kv_destroy(*rs);
}

START_TEST(check_record_select_for_phase__keeps_the_records_on_the_haplotype_of_the_phase) {
    bcf_hdr_t *hdr = mock_hdr();
    Record *first = make_record(hdr, 100, 1, 0, true, true);   // First haplotype only
    Record *second = make_record(hdr, 110, 0, 1, true, true);  // Second haplotype only
    Record *both = make_record(hdr, 120, 1, 1, true, true);    // Both
    Records block; kv_init(block);
    kv_push(Record *, block, first); kv_push(Record *, block, second); kv_push(Record *, block, both);

    Records phase_0 = Record_select_for_phase(&block, hdr, 0, false);
    ck_assert_int_eq(kv_size(phase_0), 2);
    ck_assert_ptr_eq(kv_A(phase_0, 0), first);
    ck_assert_ptr_eq(kv_A(phase_0, 1), both);

    Records phase_1 = Record_select_for_phase(&block, hdr, 1, false);
    ck_assert_int_eq(kv_size(phase_1), 2);
    ck_assert_ptr_eq(kv_A(phase_1, 0), second);
    ck_assert_ptr_eq(kv_A(phase_1, 1), both);

    kv_destroy(phase_0);
    kv_destroy(phase_1);
    free_records(&block);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_record_select_for_phase__keeps_the_order_of_the_block) {
    bcf_hdr_t *hdr = mock_hdr();
    Records block; kv_init(block);
    for (int i = 0; i < 4; i++) kv_push(Record *, block, make_record(hdr, 100 + 10 * i, 1, 1, true, true));

    Records selected = Record_select_for_phase(&block, hdr, 0, false);

    ck_assert_int_eq(kv_size(selected), 4);
    for (int i = 0; i < 4; i++) ck_assert_ptr_eq(kv_A(selected, i), kv_A(block, i));

    kv_destroy(selected);
    free_records(&block);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_record_select_for_phase__leaves_out_records_without_coding_transcripts) {
    bcf_hdr_t *hdr = mock_hdr();
    Record *without = make_record(hdr, 100, 1, 1, true, false);
    Record *with = make_record(hdr, 110, 1, 1, true, true);
    Records block; kv_init(block);
    kv_push(Record *, block, without); kv_push(Record *, block, with);

    Records selected = Record_select_for_phase(&block, hdr, 0, false);

    ck_assert_int_eq(kv_size(selected), 1);
    ck_assert_ptr_eq(kv_A(selected, 0), with);

    kv_destroy(selected);
    free_records(&block);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_record_select_for_phase__unphased_records_need_include_unphased) {
    bcf_hdr_t *hdr = mock_hdr();
    Record *unphased = make_record(hdr, 100, 0, 1, false, true);
    Records block; kv_init(block);
    kv_push(Record *, block, unphased);

    Records excluded = Record_select_for_phase(&block, hdr, 0, false);
    ck_assert_int_eq(kv_size(excluded), 0);

    Records included = Record_select_for_phase(&block, hdr, 0, true);
    ck_assert_int_eq(kv_size(included), 1);

    kv_destroy(excluded);
    kv_destroy(included);
    free_records(&block);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_record_select_for_phase__empty_block) {
    bcf_hdr_t *hdr = mock_hdr();
    Records block; kv_init(block);

    Records selected = Record_select_for_phase(&block, hdr, 0, false);
    ck_assert_int_eq(kv_size(selected), 0);

    kv_destroy(selected);
    kv_destroy(block);
    bcf_hdr_destroy(hdr);
}
END_TEST

/*
* filter_overlapping_variants
*/

// The filter as it was before the anchor-free view (in rescueranger.c), without its log line, kept to compare against.
// Returns the records it kept, in order.
static Records old_filter_overlapping_variants(const Records *vs) {
    Range last_variant_range = {0};
    int last_variant_type = 0;

    Range current_variant_range = {0};
    int current_variant_type = 0;

    Records kept; kv_init(kept);
    for (int i = 0; i < vs->n; i++) {
        bcf1_t *v = kv_A(*vs, i)->v;

        if (current_variant_range.end > last_variant_range.end) {
            last_variant_range = current_variant_range;
            last_variant_type = current_variant_type;
        }

        current_variant_range = variant_to_range(v, 0);
        Range intersection;
        int overlaps = intersect_ranges(current_variant_range, last_variant_range, &intersection);

        bool remove = true;
        if (range_len(variant_to_range(v, 1)) > 1) { // Insertion
            current_variant_type = 1;
            current_variant_range.start++;

            if (last_variant_type == 0) {
                remove = false;
            } else if (last_variant_type == -1) {
                if (current_variant_range.start <= last_variant_range.start | current_variant_range.end >= last_variant_range.end) remove = false;
            } else {
                if (last_variant_range.start != current_variant_range.start | last_variant_range.end != current_variant_range.end) remove = false;
            }
        } else if (range_len(current_variant_range) > 1) { // Deletion
            current_variant_type = -1;
            current_variant_range.start++;
            if (current_variant_range.start >= last_variant_range.end) remove = false;
        } else { // SNV
            current_variant_type = 0;
            if (!overlaps || range_len(intersection) == 0) remove = false;
        }

        if (!remove) kv_push(Record *, kept, kv_A(*vs, i));
    }

    return kept;
}

// SNVs and normalised pure insertions and deletions, as REF and ALT
#define N_FILTER_VARIANTS 5
static char *filter_variants[N_FILTER_VARIANTS][2] = {
    { "N", "x" }, { "N", "Nx" }, { "N", "Nxy" }, { "NN", "N" }, { "NNNN", "N" }
};

// Filters the records with the new filter and the old one and checks that both keep the same records
static void assert_filter_keeps_the_old_records(const int64_t *positions, const int *kinds, const int n, bcf_hdr_t *hdr) {
    Records vs; kv_init(vs);
    for (int i = 0; i < n; i++) kv_push(Record *, vs, Record_init(mock_bcf1_t(positions[i], filter_variants[kinds[i]][0], filter_variants[kinds[i]][1]), NULL));
    Records all; kv_init(all);
    kv_copy(Record *, all, vs);

    Records old = old_filter_overlapping_variants(&vs);
    filter_overlapping_variants(&vs, hdr);

    char what[128] = {0};
    for (int i = 0; i < n; i++) snprintf(what + strlen(what), sizeof(what) - strlen(what), "%s>%s at %ld; ", filter_variants[kinds[i]][0], filter_variants[kinds[i]][1], (long) positions[i]);

    ck_assert_msg(kv_size(vs) == kv_size(old), "%s: kept %zu, old %zu", what, kv_size(vs), kv_size(old));
    for (size_t i = 0; i < kv_size(vs); i++) ck_assert_msg(kv_A(vs, i) == kv_A(old, i), "%s: record %zu differs", what, i);

    for (size_t i = 0; i < kv_size(all); i++) Record_destroy(kv_A(all, i));
    kv_destroy(all); kv_destroy(vs); kv_destroy(old);
}

START_TEST(check_filter_overlapping_variants__pairs_give_the_old_result) {
    bcf_hdr_t *hdr = mock_hdr();

    for (int first = 0; first < N_FILTER_VARIANTS; first++) {
        for (int second = 0; second < N_FILTER_VARIANTS; second++) {
            for (int64_t pos = 110; pos <= 116; pos++) {
                const int64_t positions[2] = { 110, pos };
                const int kinds[2] = { first, second };
                assert_filter_keeps_the_old_records(positions, kinds, 2, hdr);
            }
        }
    }

    bcf_hdr_destroy(hdr);
} END_TEST

// Three records, so a removed record that reaches further than the one before it still counts for the next one
START_TEST(check_filter_overlapping_variants__triples_give_the_old_result) {
    bcf_hdr_t *hdr = mock_hdr();

    for (int k0 = 0; k0 < N_FILTER_VARIANTS; k0++) {
        for (int k1 = 0; k1 < N_FILTER_VARIANTS; k1++) {
            for (int k2 = 0; k2 < N_FILTER_VARIANTS; k2++) {
                for (int64_t p1 = 110; p1 <= 113; p1++) {
                    for (int64_t p2 = p1; p2 <= 116; p2++) {
                        const int64_t positions[3] = { 110, p1, p2 };
                        const int kinds[3] = { k0, k1, k2 };
                        assert_filter_keeps_the_old_records(positions, kinds, 3, hdr);
                    }
                }
            }
        }
    }

    bcf_hdr_destroy(hdr);
} END_TEST

// Keeps the records of the given REF/ALT pairs that the filter keeps, as indices into them
static int filter_kept(char *alleles[][2], const int64_t *positions, const int n, int *kept_indices) {
    bcf_hdr_t *hdr = mock_hdr();

    Records vs; kv_init(vs);
    for (int i = 0; i < n; i++) kv_push(Record *, vs, Record_init(mock_bcf1_t(positions[i], alleles[i][0], alleles[i][1]), NULL));
    Records all; kv_init(all);
    kv_copy(Record *, all, vs);

    filter_overlapping_variants(&vs, hdr);

    int n_kept = 0;
    for (size_t i = 0; i < kv_size(vs); i++) {
        for (int j = 0; j < n; j++) if (kv_A(vs, i) == kv_A(all, j)) kept_indices[n_kept++] = j;
    }

    for (size_t i = 0; i < kv_size(all); i++) Record_destroy(kv_A(all, i));
    kv_destroy(all); kv_destroy(vs);
    bcf_hdr_destroy(hdr);

    return n_kept;
}

// An MNV is a substitution of all its REF bases. The old filter took it for an insertion, so it kept an MNV on an SNV's base.
START_TEST(check_filter_overlapping_variants__mnv_on_an_snv_base_is_removed) {
    char *alleles[2][2] = { { "N", "x" }, { "NN", "xy" } };
    const int64_t positions[2] = { 110, 110 };
    int kept[2];

    ck_assert_int_eq(filter_kept(alleles, positions, 2, kept), 1);
    ck_assert_int_eq(kept[0], 0);
} END_TEST

// A complex variant that replaces a deleted base is removed. NN>Nx at 110 replaces 111, which the deletion at 109 removed;
// the old filter took it for an insertion and kept it.
START_TEST(check_filter_overlapping_variants__complex_variant_on_a_deleted_base_is_removed) {
    char *alleles[2][2] = { { "NNN", "N" }, { "NN", "Nx" } };
    const int64_t positions[2] = { 109, 110 };
    int kept[2];

    ck_assert_int_eq(filter_kept(alleles, positions, 2, kept), 1);
    ck_assert_int_eq(kept[0], 0);
} END_TEST

// Separate MNVs and complex variants are kept
START_TEST(check_filter_overlapping_variants__separate_mnv_and_complex_variant_are_kept) {
    char *alleles[3][2] = { { "NN", "xy" }, { "NNN", "Nx" }, { "N", "x" } };
    const int64_t positions[3] = { 110, 113, 117 };
    int kept[3];

    ck_assert_int_eq(filter_kept(alleles, positions, 3, kept), 3);
} END_TEST

/*
* flag_phasing_in_haplotype: the focal record gets SKIPPED_VARIANT_IN_PREDICTION / UNPHASED_VARIANT_IN_PREDICTION, on every
* transcript for the phase, whenever the haplotype holds a skipped / unphased variant in the prediction window.
*/

// The prediction radius of the window {100, 201} used in these tests
#define TEST_RADIUS 50

// A record with the given alleles, n_transcripts transcript entries and a flag slot for every transcript and phase
static Record *record_with_alleles_and_flag_slots(const int64_t pos, const char *ref, const char *alt, const int n_transcripts, const bool unphased) {
    Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);

    r->transcript_genes = malloc(sizeof(Genes));
    kv_init(*r->transcript_genes);
    for (int i = 0; i < n_transcripts; i++) kv_push(Gene *, *r->transcript_genes, NULL);

    r->flags = calloc(n_transcripts * NUM_PHASES, sizeof(uint32_t));
    r->unphased = unphased;
    return r;
}

// A substitution of N by x. A skipped record is one that is too long for the radius (51 ALT bases), without r->skipped set: 8b
// only sets that when it reaches the record, so the flag must not depend on it.
static Record *record_with_flag_slots(const int64_t pos, const int n_transcripts, const bool skipped, const bool unphased) {
    char too_long[52];
    memset(too_long, 'x', 51);
    too_long[51] = '\0';

    return record_with_alleles_and_flag_slots(pos, "N", skipped ? too_long : "x", n_transcripts, unphased);
}

static uint32_t flags_of(const Record *r, const int transcript, const int phase) {
    return r->flags[transcript * NUM_PHASES + phase];
}

// Whether an unphased variant with the alleles gets the focal record flagged, for the window {100, 201}
static bool flags_the_focal_record(const int64_t pos, const char *ref, const char *alt) {
    Record *focal = record_with_flag_slots(150, 1, false, false);
    Record *other = record_with_alleles_and_flag_slots(pos, ref, alt, 1, true);
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal); kv_push(Record *, hap, other);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);
    const bool flagged = flags_of(focal, 0, 0) == UNPHASED_VARIANT_IN_PREDICTION;

    Record_destroy(focal); Record_destroy(other);
    kv_destroy(hap);
    return flagged;
}

// A variant is in the window when a base it changes is in it, not when only the base it is attached to is. The window is
// the bases 100 to 200.
START_TEST(check_flag_phasing_in_haplotype__a_deletion_counts_by_the_bases_it_deletes) {
    ck_assert(flags_the_focal_record(99, "NN", "N"));    // Attached before the window, deletes its first base
    ck_assert(!flags_the_focal_record(98, "NN", "N"));   // Deletes 99
    ck_assert(flags_the_focal_record(199, "NNN", "N"));  // Deletes 200 and 201
    ck_assert(flags_the_focal_record(199, "NN", "N"));   // Deletes 200, the last base of the window
    ck_assert(!flags_the_focal_record(200, "NN", "N"));  // Attached on the last base, deletes 201
    ck_assert(!flags_the_focal_record(200, "NNN", "N"));
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__a_complex_indel_counts_by_the_bases_it_replaces) {
    ck_assert(flags_the_focal_record(99, "NNN", "Nx"));   // Replaces 100 and 101
    ck_assert(flags_the_focal_record(98, "NNN", "Nx"));   // Replaces 99 and 100
    ck_assert(!flags_the_focal_record(97, "NNN", "Nx"));  // Replaces 98 and 99
    ck_assert(flags_the_focal_record(199, "NNN", "Nx"));  // Replaces 200 and 201
    ck_assert(!flags_the_focal_record(200, "NNN", "Nx")); // Replaces 201 and 202
} END_TEST

// An insertion changes no base. It lies between two bases and is in the window when both are
START_TEST(check_flag_phasing_in_haplotype__an_insertion_is_in_when_both_its_neighbours_are) {
    ck_assert(flags_the_focal_record(100, "N", "Nx"));  // Between 100 and 101
    ck_assert(flags_the_focal_record(199, "N", "Nx"));  // Between 199 and 200
    ck_assert(!flags_the_focal_record(200, "N", "Nx")); // Between 200 and 201
    ck_assert(!flags_the_focal_record(99, "N", "Nx"));  // Between 99 and 100
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__a_substitution_counts_at_its_bases) {
    ck_assert(flags_the_focal_record(100, "NN", "xx"));
    ck_assert(flags_the_focal_record(99, "NN", "xx"));   // Changes 99 and 100
    ck_assert(!flags_the_focal_record(98, "NN", "xx"));
    ck_assert(flags_the_focal_record(200, "NN", "xx"));  // Changes 200 and 201
    ck_assert(!flags_the_focal_record(201, "NN", "xx"));
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__a_skipped_variant_flags_the_focal_record) {
    Record *focal = record_with_flag_slots(150, 2, false, false);
    Record *skipped = record_with_flag_slots(100, 1, true, false); // The first base of the window
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal); kv_push(Record *, hap, skipped);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 1);

    for (int t = 0; t < 2; t++) {
        ck_assert_uint_eq(flags_of(focal, t, 1), SKIPPED_VARIANT_IN_PREDICTION);
        ck_assert_uint_eq(flags_of(focal, t, 0), 0); // The other phase is left alone
    }
    ck_assert_uint_eq(flags_of(skipped, 0, 1), 0); // Only the focal record is flagged

    Record_destroy(focal); Record_destroy(skipped);
    kv_destroy(hap);
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__skipped_and_unphased_give_both_flags) {
    Record *focal = record_with_flag_slots(150, 1, false, false);
    Record *skipped = record_with_flag_slots(160, 1, true, false);
    Record *unphased = record_with_flag_slots(170, 1, false, true);
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal); kv_push(Record *, hap, skipped); kv_push(Record *, hap, unphased);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);

    ck_assert_uint_eq(flags_of(focal, 0, 0), SKIPPED_VARIANT_IN_PREDICTION | UNPHASED_VARIANT_IN_PREDICTION);

    Record_destroy(focal); Record_destroy(skipped); Record_destroy(unphased);
    kv_destroy(hap);
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__an_unphased_variant_flags_the_focal_record) {
    Record *focal = record_with_flag_slots(150, 2, false, false);
    Record *unphased = record_with_flag_slots(200, 1, false, true); // The last base of the window
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal); kv_push(Record *, hap, unphased);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 1);

    for (int t = 0; t < 2; t++) {
        ck_assert_uint_eq(flags_of(focal, t, 1), UNPHASED_VARIANT_IN_PREDICTION);
        ck_assert_uint_eq(flags_of(focal, t, 0), 0); // The other phase is left alone
    }
    ck_assert_uint_eq(flags_of(unphased, 0, 1), 0); // Only the focal record is flagged

    Record_destroy(focal); Record_destroy(unphased);
    kv_destroy(hap);
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__the_focal_record_being_unphased_or_skipped_counts) {
    Record *focal = record_with_flag_slots(150, 1, true, true);
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);

    ck_assert_uint_eq(flags_of(focal, 0, 0), SKIPPED_VARIANT_IN_PREDICTION | UNPHASED_VARIANT_IN_PREDICTION);

    Record_destroy(focal);
    kv_destroy(hap);
} END_TEST

// A skipped or unphased variant outside the prediction window is not relevant to the variant
START_TEST(check_flag_phasing_in_haplotype__a_variant_outside_the_window_is_ignored) {
    Record *focal = record_with_flag_slots(150, 1, false, false);
    Record *before = record_with_flag_slots(99, 1, true, false);
    Record *after = record_with_flag_slots(201, 1, false, true);
    Records hap; kv_init(hap);
    kv_push(Record *, hap, before); kv_push(Record *, hap, focal); kv_push(Record *, hap, after);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);

    ck_assert_uint_eq(flags_of(focal, 0, 0), 0);

    Record_destroy(focal); Record_destroy(before); Record_destroy(after);
    kv_destroy(hap);
} END_TEST

START_TEST(check_flag_phasing_in_haplotype__a_clean_haplotype_gets_nothing) {
    Record *focal = record_with_flag_slots(150, 1, false, false);
    Record *other = record_with_flag_slots(160, 1, false, false);
    Records hap; kv_init(hap);
    kv_push(Record *, hap, focal); kv_push(Record *, hap, other);

    flag_phasing_in_haplotype(focal, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);

    ck_assert_uint_eq(flags_of(focal, 0, 0), 0);

    Record_destroy(focal); Record_destroy(other);
    kv_destroy(hap);
} END_TEST

// A focal record without flag slots (no coding state) is left alone, not dereferenced
START_TEST(check_flag_phasing_in_haplotype__a_focal_record_without_flag_slots_is_left_alone) {
    Record *bare = Record_init(mock_bcf1_t(150, "N", "x"), NULL);
    bare->unphased = 1;
    Records hap; kv_init(hap);
    kv_push(Record *, hap, bare);

    flag_phasing_in_haplotype(bare, &hap, (Range) { 100, 201 }, TEST_RADIUS, 0);

    Record_destroy(bare);
    kv_destroy(hap);
} END_TEST

Suite *phase_suite() {
    Suite *s = suite_create("Phase");

    TCase *tc_hmp = mocked_tcase("has_matching_phase");
    suite_add_tcase(s, tc_hmp);

    tcase_add_test(tc_hmp, check_has_matching_phase__homozygous_is_on_both_haplotypes);
    tcase_add_test(tc_hmp, check_has_matching_phase__homozygous_reference_also_matches);
    tcase_add_test(tc_hmp, check_has_matching_phase__phased_heterozygous_is_on_the_haplotype_with_the_alt);
    tcase_add_test(tc_hmp, check_has_matching_phase__unphased_heterozygous_follows_include_unphased);

    TCase *tc_select = mocked_tcase("record_select_for_phase");
    suite_add_tcase(s, tc_select);

    tcase_add_test(tc_select, check_record_select_for_phase__keeps_the_records_on_the_haplotype_of_the_phase);
    tcase_add_test(tc_select, check_record_select_for_phase__keeps_the_order_of_the_block);
    tcase_add_test(tc_select, check_record_select_for_phase__leaves_out_records_without_coding_transcripts);
    tcase_add_test(tc_select, check_record_select_for_phase__unphased_records_need_include_unphased);
    tcase_add_test(tc_select, check_record_select_for_phase__empty_block);

    TCase *tc_filter = mocked_tcase("filter_overlapping_variants");
    suite_add_tcase(s, tc_filter);

    tcase_add_test(tc_filter, check_filter_overlapping_variants__pairs_give_the_old_result);
    tcase_add_test(tc_filter, check_filter_overlapping_variants__triples_give_the_old_result);
    tcase_add_test(tc_filter, check_filter_overlapping_variants__mnv_on_an_snv_base_is_removed);
    tcase_add_test(tc_filter, check_filter_overlapping_variants__complex_variant_on_a_deleted_base_is_removed);
    tcase_add_test(tc_filter, check_filter_overlapping_variants__separate_mnv_and_complex_variant_are_kept);

    TCase *tc_flags = mocked_tcase("flag_phasing_in_haplotype");
    suite_add_tcase(s, tc_flags);

    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_skipped_variant_flags_the_focal_record);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__skipped_and_unphased_give_both_flags);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__an_unphased_variant_flags_the_focal_record);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__the_focal_record_being_unphased_or_skipped_counts);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_variant_outside_the_window_is_ignored);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_deletion_counts_by_the_bases_it_deletes);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_complex_indel_counts_by_the_bases_it_replaces);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__an_insertion_is_in_when_both_its_neighbours_are);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_substitution_counts_at_its_bases);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_clean_haplotype_gets_nothing);
    tcase_add_test(tc_flags, check_flag_phasing_in_haplotype__a_focal_record_without_flag_slots_is_left_alone);

    return s;
}
