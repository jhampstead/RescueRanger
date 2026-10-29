#include "check_haplotype.h"

#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    for (size_t i = 0; i < kv_size(rs); i++) Record_destroy(kv_A(rs, i));
    kv_destroy(rs);

    free(seq.s);
    if (hap.s != NULL) free(hap.s);
}


/*
 * SNVs
 */
START_TEST(check_get_haplotype_sequence__snvs_within_bounds) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(6, "A", "T"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "G"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(19, "A", "C"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAATAAAGAAAAAAAAC");
} END_TEST

START_TEST(check_get_haplotype_sequence__snv_outside_bounds) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(20, "A", "T"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, seq.s);
} END_TEST

/*
 * Insertions
 */
START_TEST(check_get_haplotype_sequence__insertion_within_bounds) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "AT"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAATAAAAAAAAA");
} END_TEST

// The insertion sits between base 0 (outside) and base 1 (the gene start) and is attached to base 0, so it is not in the
// gene, as allele_offset and haplotype_offset already counted it. The version before the anchor-free view included it
// ("TAAAAAAAAAAAAAAAAAAA"); no caller passes such a record, they only pass records that overlap the gene.
START_TEST(check_get_haplotype_sequence__insertion_overlaps_bounds_start) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "A", "AT"), NULL));

    seq.s[19] = '\0';
    get_haplotype_sequence(seq, &rs, (Range) {1, 20}, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAA");
} END_TEST

// The insertion sits between base 19 (the last gene base) and base 20 (outside). It lies at the edge of the gene, so like the one
// at the gene start it is not in the gene; the version before the symmetric rule included it ("AAAAAAAAAAAAAAAAAAAATT").
START_TEST(check_get_haplotype_sequence__insertion_overlaps_bound_end) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(19, "A", "ATT"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, seq.s);
} END_TEST

// One base further in, the insertion is between two gene bases and is included
START_TEST(check_get_haplotype_sequence__insertion_next_to_bound_end) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(18, "A", "ATT"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAA" "AAAAAAAAA" "TT" "A"); // 19 bases, the insertion, the last base
} END_TEST

/*
* Deletions
*/
START_TEST(check_get_haplotype_sequence__deletion_within_bounds) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "AA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__deletion_overlaps_bound_start) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(0, "AAA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, (Range) {1, 20}, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAA");
} END_TEST

START_TEST(check_get_haplotype_sequence__deletion_overlaps_bound_end) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(18, "AAA", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "AAAAAAAAAAAAAAAAAAA");
} END_TEST

/*
* The anchor-free version against the version with the anchor
*/

// get_haplotype_sequence as it was before the anchor-free view, kept to compare against
static void old_get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap) {
    kstring_t haplotype_sequence = {0};

    int ref_idx = 0;
    int hap_idx = 0;

    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        Range ref = variant_to_range(r->v, 0);
        Range intersected_ref = {0};
        if (intersect_ranges(ref, gene, &intersected_ref) == 0) continue;
        shift_range(&ref, -gene.start);
        shift_range(&intersected_ref, -gene.start);
        int ref_segment = intersected_ref.start - ref_idx;
        if (ref_segment < 0) ref_segment = 0;

        Range alt = variant_to_range(r->v, 1);
        Range intersected_alt = {0};
        int alt_segment = 0;
        if (intersect_ranges(alt, gene, &intersected_alt) && range_len(intersected_alt) > 0) {
            shift_range(&alt, -gene.start);
            shift_range(&intersected_alt, -gene.start);

            alt_segment = range_len(alt);
        }

        char *ref_seq_pos = reference_sequence.s + ref_idx;
        int skipped_alt_bases = (intersected_alt.start - alt.start);
        char *alt_seq_pos = r->v->d.allele[1] + skipped_alt_bases;
        if (alt_segment > 0) alt_segment -= skipped_alt_bases;

        int skipped_ref_base = (range_len(ref) != range_len(alt)) & (ref.start == intersected_ref.start);

        if (ref_segment + skipped_ref_base > 0) {
            kputsn(ref_seq_pos, ref_segment + skipped_ref_base, &haplotype_sequence);
            ref_idx += ref_segment + skipped_ref_base;
            hap_idx += ref_segment + skipped_ref_base;
        }

        if (alt_segment - skipped_ref_base > 0) kputsn(alt_seq_pos + skipped_ref_base, alt_segment - skipped_ref_base, &haplotype_sequence);

        ref_idx += range_len(intersected_ref) - skipped_ref_base;
        hap_idx += alt_segment - skipped_ref_base;
    }

    int remainder = reference_sequence.l - ref_idx;
    if (remainder > 0) {
        haplotype_sequence.l = hap_idx;
        kputsn(reference_sequence.s + ref_idx, remainder, &haplotype_sequence);
    }

    *out_hap = haplotype_sequence;
}

// A gene of 36 unique bases on [100, 136), so a duplicated, missing or misplaced base shows
#define HAP_GENE_SEQ "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghij"
static const Range hap_gene = { 100, 136 };

static kstring_t hap_gene_sequence(void) {
    kstring_t s = {0};
    kputs(HAP_GENE_SEQ, &s);
    return s;
}

static void assert_same_sequence(const kstring_t *new, const kstring_t *old, const char *what) {
    ck_assert_msg(new->l == old->l && (new->l == 0 || memcmp(new->s, old->s, new->l) == 0), "%s: %.*s, old %.*s", what, (int) new->l, new->s, (int) old->l, old->s);
}

// Normalised variants of every kind (REF bases are placeholders, only their number matters)
#define N_HAP_VARIANTS 8
static const char *hap_variants[N_HAP_VARIANTS][2] = {
    { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "N", "Nxyz" }, { "NN", "N" }, { "NNNN", "N" }, { "NNN", "Nx" }, { "NN", "Nxyz" }
};

// Every single variant that overlaps the gene, at every position around both gene edges
START_TEST(check_get_haplotype_sequence__single_variants_give_the_old_sequence) {
    kstring_t ref = hap_gene_sequence();

    for (int i = 0; i < N_HAP_VARIANTS; i++) {
        for (int pos = 95; pos <= 140; pos++) {
            Records one; kv_init(one);
            kv_push(Record *, one, Record_init(mock_bcf1_t(pos, (char *) hap_variants[i][0], (char *) hap_variants[i][1]), NULL));
            if (!variant_in_region(&kv_A(one, 0)->alleles, hap_gene)) {
                Record_destroy(kv_A(one, 0));
                kv_destroy(one);
                continue;
            }

            kstring_t new = {0}, old = {0};
            get_haplotype_sequence(ref, &one, hap_gene, &new);
            old_get_haplotype_sequence(ref, &one, hap_gene, &old);

            char what[64];
            snprintf(what, sizeof(what), "%s>%s at %d", hap_variants[i][0], hap_variants[i][1], pos);
            assert_same_sequence(&new, &old, what);

            free(new.s); free(old.s);
            Record_destroy(kv_A(one, 0));
            kv_destroy(one);
        }
    }

    free(ref.s);
} END_TEST

// A variant anchored on the last gene base that only changes bases after the gene is not in the gene, so it leaves the gene
// sequence alone. The version with the anchor counted it: NN>Nxyz at 135 (the last base) added xyz after the gene.
START_TEST(check_get_haplotype_sequence__a_variant_changing_only_bases_after_the_gene_is_ignored) {
    kstring_t ref = hap_gene_sequence();
    const char *variants[][2] = { { "NN", "N" }, { "NNNN", "N" }, { "NN", "Nxyz" }, { "NNN", "Nx" } };

    for (int i = 0; i < 4; i++) {
        Records one; kv_init(one);
        kv_push(Record *, one, Record_init(mock_bcf1_t(135, (char *) variants[i][0], (char *) variants[i][1]), NULL));

        kstring_t new = {0};
        get_haplotype_sequence(ref, &one, hap_gene, &new);
        ck_assert_msg(new.l == ref.l && memcmp(new.s, ref.s, ref.l) == 0, "%s>%s at 135: %.*s", variants[i][0], variants[i][1], (int) new.l, new.s);

        free(new.s);
        Record_destroy(kv_A(one, 0));
        kv_destroy(one);
    }

    free(ref.s);
} END_TEST

// Two variants that do not share a base give the old sequence too
START_TEST(check_get_haplotype_sequence__separate_variants_give_the_old_sequence) {
    kstring_t ref = hap_gene_sequence();

    for (int i = 0; i < N_HAP_VARIANTS; i++) {
        for (int j = 0; j < N_HAP_VARIANTS; j++) {
            Records two; kv_init(two);
            kv_push(Record *, two, Record_init(mock_bcf1_t(105, (char *) hap_variants[i][0], (char *) hap_variants[i][1]), NULL));
            kv_push(Record *, two, Record_init(mock_bcf1_t(120, (char *) hap_variants[j][0], (char *) hap_variants[j][1]), NULL));

            kstring_t new = {0}, old = {0};
            get_haplotype_sequence(ref, &two, hap_gene, &new);
            old_get_haplotype_sequence(ref, &two, hap_gene, &old);

            char what[64];
            snprintf(what, sizeof(what), "%s>%s then %s>%s", hap_variants[i][0], hap_variants[i][1], hap_variants[j][0], hap_variants[j][1]);
            assert_same_sequence(&new, &old, what);

            free(new.s); free(old.s);
            Record_destroy(kv_A(two, 0));
            Record_destroy(kv_A(two, 1));
            kv_destroy(two);
        }
    }

    free(ref.s);
} END_TEST

// Where the old version went wrong: a deletion anchored on a base the deletion before it removed. AT>A at 110 removes 111
// (L), then TC>T at 111 removes 112 (M). The old version copied base 112 and removed 113 (N) instead.
START_TEST(check_get_haplotype_sequence__deletion_anchored_on_a_deleted_base) {
    kstring_t ref = hap_gene_sequence();

    Records two; kv_init(two);
    kv_push(Record *, two, Record_init(mock_bcf1_t(110, "NN", "N"), NULL));
    kv_push(Record *, two, Record_init(mock_bcf1_t(111, "NN", "N"), NULL));

    kstring_t new = {0}, old = {0};
    get_haplotype_sequence(ref, &two, hap_gene, &new);
    old_get_haplotype_sequence(ref, &two, hap_gene, &old);

    ck_assert_str_eq(new.s, "ABCDEFGHIJKNOPQRSTUVWXYZabcdefghij");
    ck_assert_str_eq(old.s, "ABCDEFGHIJKMOPQRSTUVWXYZabcdefghij");

    free(new.s); free(old.s);
    Record_destroy(kv_A(two, 0));
    Record_destroy(kv_A(two, 1));
    kv_destroy(two);
    free(ref.s);
} END_TEST

// The same for an insertion anchored on a deleted base: AT>A at 110 removes 111 (L), then T>Txy at 111 inserts xy after it.
// The inserted bases belong where L was; the old version put M before them.
START_TEST(check_get_haplotype_sequence__insertion_anchored_on_a_deleted_base) {
    kstring_t ref = hap_gene_sequence();

    Records two; kv_init(two);
    kv_push(Record *, two, Record_init(mock_bcf1_t(110, "NN", "N"), NULL));
    kv_push(Record *, two, Record_init(mock_bcf1_t(111, "N", "Nxy"), NULL));

    kstring_t new = {0}, old = {0};
    get_haplotype_sequence(ref, &two, hap_gene, &new);
    old_get_haplotype_sequence(ref, &two, hap_gene, &old);

    ck_assert_str_eq(new.s, "ABCDEFGHIJKxyMNOPQRSTUVWXYZabcdefghij");
    ck_assert_str_eq(old.s, "ABCDEFGHIJKMxyNOPQRSTUVWXYZabcdefghij");

    free(new.s); free(old.s);
    Record_destroy(kv_A(two, 0));
    Record_destroy(kv_A(two, 1));
    kv_destroy(two);
    free(ref.s);
} END_TEST

/*
* Complex cases
*/
START_TEST(check_get_haplotype_sequence__consecutive_indels) {
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "A"), NULL));
    kv_push(Record *, rs, Record_init(mock_bcf1_t(10, "A", "A"), NULL));

    get_haplotype_sequence(seq, &rs, range, &hap);

    ck_assert_str_eq(hap.s, "");
} END_TEST

/*
* Splice haplotype sequence
*
* The gene starts at 100 and has the bases A-Z then a-j on positions 100..135, so every base is unique and a
* duplicated or missing base is visible. The coding exons are [102, 110) and [120, 130), with an intron [110, 120)
* between them, and the UTRs are [100, 102) and [130, 136).
*/
#define SPLICE_GENE_START 100
#define SPLICE_GENE_SEQ "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghij"

static NonCodingRegion region(const int is_used, const Range coords, const enum UTRType utr) {
    NonCodingRegion ncr = {0};
    ncr.is_used = is_used;
    ncr.coords = coords;
    ncr.canonical_coords = coords;
    ncr.is_utr = utr;
    return ncr;
}

START_TEST(check_splice_haplotype_sequence__forward_removes_used_introns_and_utrs) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 120 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &ncrs, SPLICE_GENE_START, STRAND_FWD, &out);

    ck_assert_str_eq(out.s, "CDEFGHIJUVWXYZabcd");

    free(out.s);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_splice_haplotype_sequence__retained_intron_stays_in_the_sequence_once) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(0, (Range) { 110, 120 }, 0)); // Retained
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &ncrs, SPLICE_GENE_START, STRAND_FWD, &out);

    // Exon, intron and exon, each base exactly once
    ck_assert_str_eq(out.s, "CDEFGHIJKLMNOPQRSTUVWXYZabcd");

    free(out.s);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_splice_haplotype_sequence__reverse_gives_the_same_ascending_sequence) {
    // Regions are in transcript order, so highest first on the reverse strand
    NonCodingRegions used; kv_init(used);
    kv_push(NonCodingRegion, used, region(1, (Range) { 130, 136 }, UTR5));
    kv_push(NonCodingRegion, used, region(1, (Range) { 110, 120 }, 0));
    kv_push(NonCodingRegion, used, region(1, (Range) { 100, 102 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &used, SPLICE_GENE_START, STRAND_REV, &out);
    ck_assert_str_eq(out.s, "CDEFGHIJUVWXYZabcd");
    free(out.s);
    kv_destroy(used);

    NonCodingRegions retained; kv_init(retained);
    kv_push(NonCodingRegion, retained, region(1, (Range) { 130, 136 }, UTR5));
    kv_push(NonCodingRegion, retained, region(0, (Range) { 110, 120 }, 0)); // Retained
    kv_push(NonCodingRegion, retained, region(1, (Range) { 100, 102 }, UTR3));

    kstring_t out_retained = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &retained, SPLICE_GENE_START, STRAND_REV, &out_retained);
    ck_assert_str_eq(out_retained.s, "CDEFGHIJKLMNOPQRSTUVWXYZabcd");
    free(out_retained.s);
    kv_destroy(retained);
} END_TEST

START_TEST(check_splice_haplotype_sequence__one_retained_and_one_spliced_intron) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 115 }, 0)); // Spliced out: K-O
    kv_push(NonCodingRegion, ncrs, region(0, (Range) { 120, 125 }, 0)); // Retained: U-Y
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &ncrs, SPLICE_GENE_START, STRAND_FWD, &out);

    ck_assert_str_eq(out.s, "CDEFGHIJPQRSTUVWXYZabcd");

    free(out.s);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_splice_haplotype_sequence__empty_utr_removes_nothing) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 100 }, UTR5)); // The CDS starts at the gene start
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 120 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &ncrs, SPLICE_GENE_START, STRAND_FWD, &out);

    ck_assert_str_eq(out.s, "ABCDEFGHIJUVWXYZabcd");

    free(out.s);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_splice_haplotype_sequence__unknown_strand_appends_nothing) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));

    kstring_t out = {0};
    splice_haplotype_sequence(SPLICE_GENE_SEQ, &ncrs, SPLICE_GENE_START, STRAND_UNK, &out);

    ck_assert_int_eq(out.l, 0);

    free(out.s);
    kv_destroy(ncrs);
} END_TEST

// The builder that was in main() before the ranges were made half-open, kept here to show what it produced. Its regions
// are inclusive: an intron is {donor base, acceptor base}, a UTR3 starts at the last coding base, and a retained intron
// is the zero-length range {donor base, donor base}.
static void old_forward_splice_sequence(const char *hap_seq, const Range *regions, const int n_regions, const int gene_start, kstring_t *out) {
    int seq_idx = 1;
    for (int j = 0; j < n_regions; j++) {
        int intron_start = regions[j].start - gene_start + 1;
        int intron_end = regions[j].end - gene_start;

        if (intron_start < seq_idx) continue;

        kputsn(hap_seq + seq_idx, intron_start - seq_idx, out);
        seq_idx = intron_end;
    }
}

START_TEST(check_old_forward_splice_sequence__duplicated_the_donor_base_of_a_retained_intron) {
    const Range spliced[] = { {100, 102}, {109, 120}, {129, 135} };
    kstring_t out_spliced = {0};
    old_forward_splice_sequence(SPLICE_GENE_SEQ, spliced, 3, SPLICE_GENE_START, &out_spliced);
    ck_assert_str_eq(out_spliced.s, "CDEFGHIJUVWXYZabcd"); // Same as the current builder
    free(out_spliced.s);

    const Range retained[] = { {100, 102}, {109, 109}, {129, 135} };
    kstring_t out_retained = {0};
    old_forward_splice_sequence(SPLICE_GENE_SEQ, retained, 3, SPLICE_GENE_START, &out_retained);

    // The last base of the upstream exon, J, appears twice. The current builder gives CDEFGHIJKLMNOPQRSTUVWXYZabcd
    ck_assert_str_eq(out_retained.s, "CDEFGHIJJKLMNOPQRSTUVWXYZabcd");
    free(out_retained.s);
} END_TEST

/*
* Enumerate region combinations
*
* Every region is identified by its start; regions are one base long.
*/
static NonCodingRegions group_of(const int64_t *starts, const size_t n) {
    NonCodingRegions group; kv_init(group);
    for (size_t i = 0; i < n; i++) kv_push(NonCodingRegion, group, region(1, (Range) { starts[i], starts[i] + 1 }, 0));
    return group;
}

static void assert_set_starts(const NonCodingRegions *set, const int64_t *starts, const size_t n) {
    ck_assert_uint_eq(kv_size(*set), n);
    for (size_t i = 0; i < n; i++) ck_assert_int_eq(kv_A(*set, i).coords.start, starts[i]);
}

static void destroy_groups(NonCodingRegions *groups, const size_t n) {
    for (size_t i = 0; i < n; i++) kv_destroy(groups[i]);
}

START_TEST(check_enumerate_region_combinations__no_groups_gives_no_combinations) {
    NonCodingRegionSets sets;
    enumerate_region_combinations(NULL, 0, &sets);

    ck_assert_uint_eq(kv_size(sets), 0);

    NonCodingRegionSets_destroy(&sets);
} END_TEST

START_TEST(check_enumerate_region_combinations__only_empty_groups_gives_no_combinations) {
    NonCodingRegions groups[2] = { group_of(NULL, 0), group_of(NULL, 0) };

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups, 2, &sets);

    ck_assert_uint_eq(kv_size(sets), 0);

    NonCodingRegionSets_destroy(&sets);
    destroy_groups(groups, 2);
} END_TEST

START_TEST(check_enumerate_region_combinations__one_group_gives_one_combination_per_region) {
    NonCodingRegions groups[1] = { group_of((int64_t[]) { 10, 20, 30 }, 3) };

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups, 1, &sets);

    ck_assert_uint_eq(kv_size(sets), 3);
    assert_set_starts(&kv_A(sets, 0), (int64_t[]) { 10 }, 1);
    assert_set_starts(&kv_A(sets, 1), (int64_t[]) { 20 }, 1);
    assert_set_starts(&kv_A(sets, 2), (int64_t[]) { 30 }, 1);

    NonCodingRegionSets_destroy(&sets);
    destroy_groups(groups, 1);
} END_TEST

// The order of the code this replaced: the choice in the first group varies fastest
START_TEST(check_enumerate_region_combinations__first_group_varies_fastest) {
    NonCodingRegions groups[3] = {
        group_of((int64_t[]) { 1 }, 1),
        group_of((int64_t[]) { 10, 11 }, 2),
        group_of((int64_t[]) { 20, 21, 22 }, 3)
    };

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups, 3, &sets);

    ck_assert_uint_eq(kv_size(sets), 6);
    assert_set_starts(&kv_A(sets, 0), (int64_t[]) { 1, 10, 20 }, 3);
    assert_set_starts(&kv_A(sets, 1), (int64_t[]) { 1, 11, 20 }, 3);
    assert_set_starts(&kv_A(sets, 2), (int64_t[]) { 1, 10, 21 }, 3);
    assert_set_starts(&kv_A(sets, 3), (int64_t[]) { 1, 11, 21 }, 3);
    assert_set_starts(&kv_A(sets, 4), (int64_t[]) { 1, 10, 22 }, 3);
    assert_set_starts(&kv_A(sets, 5), (int64_t[]) { 1, 11, 22 }, 3);

    NonCodingRegionSets_destroy(&sets);
    destroy_groups(groups, 3);
} END_TEST

START_TEST(check_enumerate_region_combinations__empty_groups_are_skipped) {
    NonCodingRegions groups[4] = {
        group_of(NULL, 0),
        group_of((int64_t[]) { 1 }, 1),
        group_of(NULL, 0),
        group_of((int64_t[]) { 10, 11 }, 2)
    };

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups, 4, &sets);

    ck_assert_uint_eq(kv_size(sets), 2);
    assert_set_starts(&kv_A(sets, 0), (int64_t[]) { 1, 10 }, 2);
    assert_set_starts(&kv_A(sets, 1), (int64_t[]) { 1, 11 }, 2);

    NonCodingRegionSets_destroy(&sets);
    destroy_groups(groups, 4);
} END_TEST

START_TEST(check_enumerate_region_combinations__regions_are_copied_by_value) {
    NonCodingRegions groups[2] = {
        group_of((int64_t[]) { 1 }, 1),
        group_of((int64_t[]) { 10, 11 }, 2)
    };

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups, 2, &sets);

    kv_A(kv_A(sets, 0), 0).is_used = 0;
    ck_assert_int_eq(kv_A(groups[0], 0).is_used, 1);
    ck_assert_int_eq(kv_A(kv_A(sets, 1), 0).is_used, 1);

    NonCodingRegionSets_destroy(&sets);
    destroy_groups(groups, 2);
} END_TEST

/*
* Isoform helpers
*/
START_TEST(check_codon_start__rounds_down_to_a_multiple_of_three) {
    ck_assert_int_eq(codon_start(0), 0);
    ck_assert_int_eq(codon_start(2), 0);
    ck_assert_int_eq(codon_start(3), 3);
    ck_assert_int_eq(codon_start(7), 6);
} END_TEST

START_TEST(check_codon_window__runs_from_the_codon_start_to_offset_plus_length) {
    const Range window = codon_window(7, 2);
    ck_assert_int_eq(window.start, 6);
    ck_assert_int_eq(window.end, 9);

    const Range one_base = codon_window(4, 1);
    ck_assert_int_eq(one_base.start, 3);
    ck_assert_int_eq(one_base.end, 5);
} END_TEST

// A deletion of the first or last coding bases puts its position just outside the coding sequence, at -1 or at its length
START_TEST(check_codon_window_within__keeps_the_window_inside_the_sequence) {
    // Inside the sequence it is codon_window
    Range window = codon_window_within(7, 2, 50);
    ck_assert_int_eq(window.start, 6);
    ck_assert_int_eq(window.end, 9);

    // Just before the sequence: from its first codon
    window = codon_window_within(-1, 3, 50);
    ck_assert_int_eq(window.start, 0);
    ck_assert_int_eq(window.end, 3);

    // Far before it, as the old positions of a deletion in the UTR were
    window = codon_window_within(-11, 1, 50);
    ck_assert_int_eq(window.start, 0);
    ck_assert_int_eq(window.end, 1);

    // Just after the sequence: the last base's codon, and no further than the end
    window = codon_window_within(50, 3, 50);
    ck_assert_int_eq(window.start, 48);
    ck_assert_int_eq(window.end, 50);

    // A window that would run past the end stops there
    window = codon_window_within(47, 6, 50);
    ck_assert_int_eq(window.start, 45);
    ck_assert_int_eq(window.end, 50);

    // No sequence
    window = codon_window_within(0, 3, 0);
    ck_assert_int_eq(window.start, 0);
    ck_assert_int_eq(window.end, 0);
} END_TEST

// The sequence is the whole string
#define HAS_STOP(seq, start, length) has_stop_codon((seq), strlen(seq), (start), (length))

START_TEST(check_has_stop_codon__finds_each_stop_codon) {
    ck_assert(HAS_STOP("TAG", 0, 3));
    ck_assert(HAS_STOP("TGA", 0, 3));
    ck_assert(HAS_STOP("TAA", 0, 3));
    ck_assert(!HAS_STOP("ATG", 0, 3));
} END_TEST

START_TEST(check_has_stop_codon__reads_codons_from_start_in_steps_of_three) {
    ck_assert(HAS_STOP("ATGTAGCCC", 0, 6));
    ck_assert(HAS_STOP("ATGTAGCCC", 3, 3));
    ck_assert(!HAS_STOP("ATGTAGCCC", 0, 3));
    ck_assert(!HAS_STOP("ATGTAGCCC", 6, 3));
} END_TEST

START_TEST(check_has_stop_codon__ignores_stops_out_of_frame) {
    ck_assert(!HAS_STOP("ATAGCC", 0, 6)); // TAG at 1
} END_TEST

START_TEST(check_has_stop_codon__reads_the_whole_last_codon) {
    ck_assert(HAS_STOP("ATGTAG", 3, 1));
    ck_assert(HAS_STOP("ATGTAG", 0, 4));
} END_TEST

START_TEST(check_has_stop_codon__zero_length_reads_nothing) {
    ck_assert(!HAS_STOP("TAG", 0, 0));
} END_TEST

// A scan that runs past the end of the sequence stops there: the bases after it (here a TAG behind the terminator) are not read
START_TEST(check_has_stop_codon__does_not_read_past_the_sequence) {
    const char seq[] = "ATGCC\0TAGTAG";

    ck_assert(!has_stop_codon(seq, 5, 3, 9));
    ck_assert(has_stop_codon(seq, 12, 3, 9)); // The same bytes, as part of the sequence
} END_TEST

// A codon that is not entirely inside the sequence is not read
START_TEST(check_has_stop_codon__a_codon_over_the_end_is_not_read) {
    ck_assert(!has_stop_codon("ATGTAG", 5, 3, 3)); // TAG needs the base at 5
    ck_assert(has_stop_codon("ATGTAG", 6, 3, 3));
} END_TEST

START_TEST(check_has_stop_codon__a_start_outside_the_sequence_reads_nothing) {
    ck_assert(!has_stop_codon("TAGTAG", 6, 6, 3));
    ck_assert(!has_stop_codon("TAGTAG", 9, 9, 3));
    ck_assert(!has_stop_codon("TAG", 3, -3, 3)); // And not before it
    ck_assert(!has_stop_codon("", 0, 0, 3));
} END_TEST

/*
* first_stop_codon: the same scan as has_stop_codon, with the position of the stop. The sequence is the whole string.
*/
#define FIRST_STOP(seq, start, length, ignore_final) first_stop_codon((seq), strlen(seq), (start), (length), (ignore_final))

START_TEST(check_first_stop_codon__is_the_index_of_the_first_base_of_the_stop) {
    ck_assert_int_eq(FIRST_STOP("TAA", 0, 3, false), 0);
    ck_assert_int_eq(FIRST_STOP("ATGTAGCCC", 0, 9, false), 3);
    ck_assert_int_eq(FIRST_STOP("ATGCCCTGA", 0, 9, false), 6);
    ck_assert_int_eq(FIRST_STOP("ATGCCCTAA", 0, 9, false), 6);
} END_TEST

START_TEST(check_first_stop_codon__is_minus_one_without_a_stop) {
    ck_assert_int_eq(FIRST_STOP("ATGCCCGGG", 0, 9, false), -1);
    ck_assert_int_eq(FIRST_STOP("", 0, 3, false), -1);
} END_TEST

START_TEST(check_first_stop_codon__gives_the_first_of_several) {
    ck_assert_int_eq(FIRST_STOP("ATGTAGTGATAA", 0, 12, false), 3);
    ck_assert_int_eq(FIRST_STOP("ATGTAGTGATAA", 6, 6, false), 6);
} END_TEST

START_TEST(check_first_stop_codon__reads_from_start_in_steps_of_three_up_to_the_length) {
    ck_assert_int_eq(FIRST_STOP("ATGTAGCCC", 3, 3, false), 3);
    ck_assert_int_eq(FIRST_STOP("ATGTAGCCC", 6, 3, false), -1);   // The stop is before start
    ck_assert_int_eq(FIRST_STOP("ATGCCCTAA", 0, 6, false), -1);   // And after the window
    ck_assert_int_eq(FIRST_STOP("ATAGCC", 0, 6, false), -1);      // TAG at 1 is out of frame
} END_TEST

START_TEST(check_first_stop_codon__does_not_read_a_codon_that_is_not_entirely_in_the_sequence) {
    ck_assert_int_eq(FIRST_STOP("ATGTAG", 3, 3, false), 3);
    ck_assert_int_eq(first_stop_codon("ATGTAG", 5, 3, 3, false), -1);    // TAG needs the base at 5
    ck_assert_int_eq(first_stop_codon("TAGTAG", 6, 6, 3, false), -1);    // Starts outside
    ck_assert_int_eq(first_stop_codon("TAG", 3, -3, 3, false), -1);      // Before the sequence
} END_TEST

// The original stop codon is the last three bases of the isoform sequence. It does not count as a new stop
START_TEST(check_first_stop_codon__can_leave_out_the_stop_codon_that_ends_the_sequence) {
    ck_assert_int_eq(FIRST_STOP("ATGCCCTAA", 0, 9, true), -1);
    ck_assert_int_eq(FIRST_STOP("ATGCCCTAA", 0, 9, false), 6);

    // A stop before it still counts, and so does one that only touches the end when the sequence is longer than the scan
    ck_assert_int_eq(FIRST_STOP("ATGTAGCCCTAA", 0, 12, true), 3);
    ck_assert_int_eq(FIRST_STOP("ATGCCCTAAGGG", 0, 9, true), 6);
    ck_assert_int_eq(first_stop_codon("ATGCCCTAAGGG", 9, 0, 9, true), -1);   // The sequence ends at 9: TAA is its last codon
} END_TEST

START_TEST(check_first_stop_codon__a_frameshifted_sequence_that_reads_on_to_the_original_stop_has_no_new_stop) {
    // After a frameshift the codons from the variant read CCC GGG TAA: the stop is the original one at the end
    ck_assert_int_eq(FIRST_STOP("ATGCCCGGGTAA", 3, 9, true), -1);
    ck_assert_int_eq(FIRST_STOP("ATGCCCGGGTAA", 3, 9, false), 9);
} END_TEST

START_TEST(check_has_stop_codon__agrees_with_first_stop_codon) {
    const char *seqs[] = { "", "TAG", "ATGTAGCCC", "ATAGCC", "ATGCCCTAA", "ATGTAGTGATAA" };

    for (size_t si = 0; si < sizeof(seqs) / sizeof(seqs[0]); si++) {
        const int64_t n = (int64_t) strlen(seqs[si]);
        for (int64_t start = -4; start <= n + 4; start++) {
            for (int64_t length = 0; length <= n + 4; length++) {
                ck_assert_msg(has_stop_codon(seqs[si], n, start, length) == (first_stop_codon(seqs[si], n, start, length, false) >= 0),
                              "%s, start %ld, length %ld", seqs[si], (long) start, (long) length);
            }
        }
    }
} END_TEST

/*
* The coordinates of the isoform sequence. The layout of the splice tests: the gene starts at 100 and has 36 bases, the UTRs are
* [100, 102) and [130, 136) and the intron is [110, 120). The regions are given in transcript order.
*/
static NonCodingRegions splice_layout(const int intron_used, const bool reverse) {
    NonCodingRegions ncrs; kv_init(ncrs);

    if (!reverse) {
        kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
        kv_push(NonCodingRegion, ncrs, region(intron_used, (Range) { 110, 120 }, 0));
        kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));
    } else {
        kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR5));
        kv_push(NonCodingRegion, ncrs, region(intron_used, (Range) { 110, 120 }, 0));
        kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR3));
    }

    return ncrs;
}

static void check_segment(const Ranges *segments, const size_t i, const int64_t start, const int64_t end) {
    ck_assert_uint_gt(kv_size(*segments), i);
    ck_assert_int_eq(kv_A(*segments, i).start, start);
    ck_assert_int_eq(kv_A(*segments, i).end, end);
}

START_TEST(check_isoform_segments__spliced_introns_and_utrs_leave_the_exons) {
    NonCodingRegions ncrs = splice_layout(1, false);
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_FWD);

    ck_assert_int_eq(kv_size(segments), 2);
    check_segment(&segments, 0, 2, 10);    // C-J
    check_segment(&segments, 1, 20, 30);   // U-d
    ck_assert_int_eq(isoform_length(&segments), 18);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__a_retained_intron_is_part_of_one_segment) {
    NonCodingRegions ncrs = splice_layout(0, false);
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_FWD);

    ck_assert_int_eq(kv_size(segments), 1);
    check_segment(&segments, 0, 2, 30);
    ck_assert_int_eq(isoform_length(&segments), 28);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__one_retained_and_one_spliced_intron) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 115 }, 0));
    kv_push(NonCodingRegion, ncrs, region(0, (Range) { 120, 125 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_FWD);

    ck_assert_int_eq(kv_size(segments), 2);
    check_segment(&segments, 0, 2, 10);
    check_segment(&segments, 1, 15, 30);
    ck_assert_int_eq(isoform_length(&segments), 23);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__reverse_gives_the_same_ascending_segments) {
    NonCodingRegions ncrs = splice_layout(1, true);
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_REV);

    ck_assert_int_eq(kv_size(segments), 2);
    check_segment(&segments, 0, 2, 10);
    check_segment(&segments, 1, 20, 30);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__an_empty_utr_at_the_gene_start_leaves_nothing_out) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 100 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 120 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_FWD);

    ck_assert_int_eq(kv_size(segments), 2);
    check_segment(&segments, 0, 0, 10);
    check_segment(&segments, 1, 20, 30);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__a_region_that_starts_before_the_end_of_the_previous_one_is_skipped) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 120 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 115, 125 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_FWD);

    ck_assert_int_eq(kv_size(segments), 2);
    check_segment(&segments, 0, 2, 10);
    check_segment(&segments, 1, 20, 30);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_segments__an_unknown_strand_has_none) {
    NonCodingRegions ncrs = splice_layout(1, false);
    Ranges segments = isoform_segments(&ncrs, SPLICE_GENE_START, STRAND_UNK);

    ck_assert_int_eq(kv_size(segments), 0);
    ck_assert_int_eq(isoform_length(&segments), 0);

    kv_destroy(segments); kv_destroy(ncrs);
} END_TEST

// hap_range_to_isoform on the layout above: the isoform of the forward strand is C-J then U-d (18 bases)
static Range spliced_range(const NonCodingRegions *ncrs, const uint32_t strand, const int64_t start, const int64_t end) {
    Ranges segments = isoform_segments(ncrs, SPLICE_GENE_START, strand);
    const Range range = hap_range_to_isoform(&segments, SPLICE_GENE_START, strand, (Range) { start, end });
    kv_destroy(segments);
    return range;
}

// The macro parameters must not be called start or end: they would replace the members of got
#define CHECK_SPLICED(ncrs, strand, from, to, expected_start, expected_end) do { \
    const Range got = spliced_range((ncrs), (strand), (from), (to)); \
    ck_assert_msg(got.start == (expected_start) && got.end == (expected_end), "[%d, %d): got [%ld, %ld), expected [%d, %d)", \
                  (int) (from), (int) (to), (long) got.start, (long) got.end, (int) (expected_start), (int) (expected_end)); \
} while (0)

START_TEST(check_hap_range_to_isoform__forward) {
    NonCodingRegions ncrs = splice_layout(1, false);

    CHECK_SPLICED(&ncrs, STRAND_FWD, 103, 106, 1, 4);      // Inside the first exon
    CHECK_SPLICED(&ncrs, STRAND_FWD, 108, 122, 6, 10);     // Across the spliced intron: I, J and U, V
    CHECK_SPLICED(&ncrs, STRAND_FWD, 112, 125, 8, 13);     // Starting in the intron: from U
    CHECK_SPLICED(&ncrs, STRAND_FWD, 112, 118, 8, 8);      // Only intron: empty, where its bases would have been
    CHECK_SPLICED(&ncrs, STRAND_FWD, 100, 136, 0, 18);     // The whole gene
    CHECK_SPLICED(&ncrs, STRAND_FWD, 100, 102, 0, 0);      // Only the 5' UTR
    CHECK_SPLICED(&ncrs, STRAND_FWD, 130, 136, 18, 18);    // Only the 3' UTR
    CHECK_SPLICED(&ncrs, STRAND_FWD, 95, 140, 0, 18);      // Past both ends of the gene
    CHECK_SPLICED(&ncrs, STRAND_FWD, 105, 105, 3, 3);      // Empty
    CHECK_SPLICED(&ncrs, STRAND_FWD, 106, 103, 4, 4);      // Inverted

    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__reverse_counts_from_the_other_end) {
    NonCodingRegions ncrs = splice_layout(1, true);

    CHECK_SPLICED(&ncrs, STRAND_REV, 103, 106, 14, 17);    // Bases 1 to 3 from the left are 14 to 16 from the right
    CHECK_SPLICED(&ncrs, STRAND_REV, 108, 122, 8, 12);
    CHECK_SPLICED(&ncrs, STRAND_REV, 100, 136, 0, 18);
    CHECK_SPLICED(&ncrs, STRAND_REV, 100, 102, 18, 18);
    CHECK_SPLICED(&ncrs, STRAND_REV, 112, 118, 10, 10);

    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__a_retained_intron_is_counted) {
    NonCodingRegions ncrs = splice_layout(0, false);

    CHECK_SPLICED(&ncrs, STRAND_FWD, 108, 122, 6, 20);     // I to V, with the whole intron between them
    CHECK_SPLICED(&ncrs, STRAND_FWD, 110, 120, 8, 18);     // The intron itself
    CHECK_SPLICED(&ncrs, STRAND_FWD, 100, 136, 0, 28);

    kv_destroy(ncrs);
} END_TEST

// Every start and end around a layout against a probe: the base at offset i of the gene is the character 33 + i, so the isoform
// sequence says where each of its bases came from. The expected range is read from it, in the direction of the transcript.
static void check_conversion_against_the_probe(const NonCodingRegions *ncrs, const uint32_t strand) {
    char probe[37];
    for (int i = 0; i < 36; i++) probe[i] = (char) (33 + i);
    probe[36] = '\0';

    kstring_t spliced = {0};
    splice_haplotype_sequence(probe, ncrs, SPLICE_GENE_START, strand, &spliced);
    const int64_t length = (int64_t) spliced.l;

    // The bases in the direction of the transcript
    char *oriented = malloc(length + 1);
    for (int64_t p = 0; p < length; p++) oriented[p] = strand == STRAND_REV ? spliced.s[length - 1 - p] : spliced.s[p];

    Ranges segments = isoform_segments(ncrs, SPLICE_GENE_START, strand);
    ck_assert_int_eq(isoform_length(&segments), length);

    for (int64_t start = SPLICE_GENE_START - 4; start <= SPLICE_GENE_START + 40; start++) {
        for (int64_t end = start - 2; end <= SPLICE_GENE_START + 40; end++) {
            const Range got = hap_range_to_isoform(&segments, SPLICE_GENE_START, strand, (Range) { start, end });

            const int64_t lo = start - SPLICE_GENE_START, hi = end - SPLICE_GENE_START;
            int64_t first = -1, last = -1, before = 0;
            for (int64_t p = 0; p < length; p++) {
                const int64_t origin = oriented[p] - 33;
                if (origin < lo) before++;
                if (origin >= lo && origin < hi) {
                    if (first < 0) first = p;
                    last = p;
                }
            }

            Range expected;
            if (first >= 0) expected = (Range) { first, last + 1 };
            else expected = (Range) { strand == STRAND_REV ? length - before : before, strand == STRAND_REV ? length - before : before };

            ck_assert_msg(got.start == expected.start && got.end == expected.end,
                          "strand %u, [%ld, %ld): got [%ld, %ld), expected [%ld, %ld)", strand, (long) start, (long) end,
                          (long) got.start, (long) got.end, (long) expected.start, (long) expected.end);
        }
    }

    free(oriented); free(spliced.s); kv_destroy(segments);
}

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_spliced_forward_gene) {
    NonCodingRegions ncrs = splice_layout(1, false);
    check_conversion_against_the_probe(&ncrs, STRAND_FWD);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_forward_gene_with_a_retained_intron) {
    NonCodingRegions ncrs = splice_layout(0, false);
    check_conversion_against_the_probe(&ncrs, STRAND_FWD);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_spliced_reverse_gene) {
    NonCodingRegions ncrs = splice_layout(1, true);
    check_conversion_against_the_probe(&ncrs, STRAND_REV);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_reverse_gene_with_a_retained_intron) {
    NonCodingRegions ncrs = splice_layout(0, true);
    check_conversion_against_the_probe(&ncrs, STRAND_REV);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_with_one_retained_and_one_spliced_intron) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 110, 115 }, 0));
    kv_push(NonCodingRegion, ncrs, region(0, (Range) { 120, 125 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 130, 136 }, UTR3));
    check_conversion_against_the_probe(&ncrs, STRAND_FWD);
    kv_destroy(ncrs);
} END_TEST

START_TEST(check_hap_range_to_isoform__every_range_agrees_with_the_probe_with_two_spliced_introns_and_an_empty_utr) {
    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 100 }, UTR5));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 106, 112 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 118, 126 }, 0));
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 134, 136 }, UTR3));
    check_conversion_against_the_probe(&ncrs, STRAND_FWD);
    kv_destroy(ncrs);

    NonCodingRegions reverse; kv_init(reverse);
    kv_push(NonCodingRegion, reverse, region(1, (Range) { 134, 136 }, UTR5));
    kv_push(NonCodingRegion, reverse, region(1, (Range) { 118, 126 }, 0));
    kv_push(NonCodingRegion, reverse, region(1, (Range) { 106, 112 }, 0));
    kv_push(NonCodingRegion, reverse, region(1, (Range) { 100, 100 }, UTR3));
    check_conversion_against_the_probe(&reverse, STRAND_REV);
    kv_destroy(reverse);
} END_TEST

START_TEST(check_position_in_ranges__is_half_open) {
    Ranges ranges; kv_init(ranges);
    kv_push(Range, ranges, ((Range) { 10, 20 }));
    kv_push(Range, ranges, ((Range) { 30, 40 }));

    ck_assert(!position_in_ranges(&ranges, 9));
    ck_assert(position_in_ranges(&ranges, 10));
    ck_assert(position_in_ranges(&ranges, 19));
    ck_assert(!position_in_ranges(&ranges, 20));
    ck_assert(!position_in_ranges(&ranges, 25));
    ck_assert(position_in_ranges(&ranges, 30));
    ck_assert(position_in_ranges(&ranges, 39));
    ck_assert(!position_in_ranges(&ranges, 40));
    ck_assert(!position_in_ranges(&ranges, -1));

    kv_destroy(ranges);

    Ranges none; kv_init(none);
    ck_assert(!position_in_ranges(&none, 0));
    kv_destroy(none);
} END_TEST

/*
* The NMD escape regions in the coordinates of an isoform. Forward transcript: exons [1000, 1200), [2000, 2600), [3000, 3300) with
* the coding regions [1100, 1200), [2000, 2600), [3000, 3150) and the gene [1000, 3300). Reverse transcript: exons [1000, 1300),
* [2000, 2600), [3000, 3200) with [1150, 1300), [2000, 2600), [3000, 3100) and the gene [1000, 3200); the transcript starts at the
* high end. The coding sequence of both has 850 bases: 100, 600 and 150 in transcript order.
*/
typedef struct NmdLayout {
    Transcript *transcript;
    Range gene;
    NonCodingRegions ncrs;   // In transcript order, in reference coordinates
} NmdLayout;

// retained: 0 splices out both introns, 1 retains the first intron in transcript order, 2 retains the last one
static NmdLayout nmd_layout(const bool reverse, const int retained) {
    NmdLayout layout = {0};
    kv_init(layout.ncrs);

    if (!reverse) {
        Range cds[] = { {1100, 1199}, {2000, 2599}, {3000, 3149} };
        int constitutive[] = { 1, 1, 1 };
        Range exons[] = { {1000, 1199}, {2000, 2599}, {3000, 3299} };
        TranscriptData data = { .id = 1, .start = 1000, .end = 3299, .strand = STRAND_FWD,
                                .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };
        layout.transcript = Transcript_build(&data);
        layout.gene = (Range) { 1000, 3300 };

        kv_push(NonCodingRegion, layout.ncrs, region(1, (Range) { 1000, 1100 }, UTR5));
        kv_push(NonCodingRegion, layout.ncrs, region(retained != 1, (Range) { 1200, 2000 }, 0));
        kv_push(NonCodingRegion, layout.ncrs, region(retained != 2, (Range) { 2600, 3000 }, 0));
        kv_push(NonCodingRegion, layout.ncrs, region(1, (Range) { 3150, 3300 }, UTR3));
    } else {
        Range cds[] = { {1150, 1299}, {2000, 2599}, {3000, 3099} };
        int constitutive[] = { 1, 1, 1 };
        Range exons[] = { {1000, 1299}, {2000, 2599}, {3000, 3199} };
        TranscriptData data = { .id = 2, .start = 1000, .end = 3199, .strand = STRAND_REV,
                                .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };
        layout.transcript = Transcript_build(&data);
        layout.gene = (Range) { 1000, 3200 };

        kv_push(NonCodingRegion, layout.ncrs, region(1, (Range) { 3100, 3200 }, UTR5));
        kv_push(NonCodingRegion, layout.ncrs, region(retained != 1, (Range) { 2600, 3000 }, 0));
        kv_push(NonCodingRegion, layout.ncrs, region(retained != 2, (Range) { 1300, 2000 }, 0));
        kv_push(NonCodingRegion, layout.ncrs, region(1, (Range) { 1000, 1150 }, UTR3));
    }

    return layout;
}

static void nmd_layout_destroy(NmdLayout *layout) {
    Transcript_destroy(layout->transcript);
    kv_destroy(layout->ncrs);
}

// What evaluate_isoform does: the non-coding regions go to the haplotype, they make the segments, and the escape regions follow
static NmdRegions nmd_regions_in_isoform(const NmdLayout *layout, const Records *hap_in) {
    NonCodingRegions hap_ncrs; kv_init(hap_ncrs);
    for (size_t k = 0; k < kv_size(layout->ncrs); k++) {
        NonCodingRegion ncr = kv_A(layout->ncrs, k);
        ncr.coords = region_to_haplotype(hap_in, layout->gene, ncr.coords, ncr.is_utr != 0);
        kv_push(NonCodingRegion, hap_ncrs, ncr);
    }

    Ranges segments = isoform_segments(&hap_ncrs, layout->gene.start, layout->transcript->strand);
    NmdRegions regions = Transcript_nmd_regions(layout->transcript);
    NmdRegions escape;
    nmd_escape_ranges_in_isoform(&regions, hap_in, layout->gene, &segments, layout->transcript->strand, &escape);

    NmdRegions_destroy(&regions);
    kv_destroy(segments);
    kv_destroy(hap_ncrs);

    return escape;
}

// The same, with the ranges of the three kinds one after the other: start-proximal, long exon, distal
static void nmd_ranges(const NmdLayout *layout, const Records *hap_in, Ranges *out) {
    NmdRegions escape = nmd_regions_in_isoform(layout, hap_in);

    kv_init(*out);
    const Ranges *lists[3] = { &escape.start_proximal, &escape.long_exon, &escape.distal };
    for (int l = 0; l < 3; l++) {
        for (size_t i = 0; i < kv_size(*lists[l]); i++) kv_push(Range, *out, kv_A(*lists[l], i));
    }

    NmdRegions_destroy(&escape);
}

static void check_isoform_range(const Ranges *ranges, const size_t i, const int64_t start, const int64_t end) {
    ck_assert_uint_gt(kv_size(*ranges), i);
    ck_assert_msg(kv_A(*ranges, i).start == start && kv_A(*ranges, i).end == end, "range %zu: got [%ld, %ld), expected [%ld, %ld)",
                  i, (long) kv_A(*ranges, i).start, (long) kv_A(*ranges, i).end, (long) start, (long) end);
}

START_TEST(check_nmd_escape_ranges_in_isoform__forward_without_variants) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // The isoform is the first exon at 0 to 99, the second at 100 to 699 and the third at 700 to 849
    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 0, 100);     // Start-proximal: the first exon...
    check_isoform_range(&ranges, 1, 100, 150);   // ...and 50 bases of the second
    check_isoform_range(&ranges, 2, 100, 299);   // Long exon: the first 199 bases of the second exon
    check_isoform_range(&ranges, 3, 650, 850);   // Distal: the last 50 bases of the second exon and all of the last

    kv_destroy(ranges); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__every_kind_keeps_its_own_list) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    NmdRegions escape = nmd_regions_in_isoform(&layout, &hap_in);

    ck_assert_int_eq(kv_size(escape.start_proximal), 2);
    check_isoform_range(&escape.start_proximal, 0, 0, 100);
    check_isoform_range(&escape.start_proximal, 1, 100, 150);
    ck_assert_int_eq(kv_size(escape.long_exon), 1);
    check_isoform_range(&escape.long_exon, 0, 100, 299);
    ck_assert_int_eq(kv_size(escape.distal), 1);
    check_isoform_range(&escape.distal, 0, 650, 850);

    NmdRegions_destroy(&escape); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__reverse_without_variants_is_the_same_in_transcript_order) {
    NmdLayout layout = nmd_layout(true, 0);
    Records hap_in; kv_init(hap_in);
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 100, 150);   // Start-proximal: the first 50 coding bases of the second exon, which is the top of it
    check_isoform_range(&ranges, 1, 0, 100);     // The first exon in transcript order
    check_isoform_range(&ranges, 2, 100, 299);
    check_isoform_range(&ranges, 3, 650, 850);

    kv_destroy(ranges); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_deletion_inside_a_region_shortens_it) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, Record_init(mock_bcf1_t(2100, "NNNN", "N"), NULL));   // Deletes 2101 to 2103, in the long exon region
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 0, 100);
    check_isoform_range(&ranges, 1, 100, 150);   // Before the deletion
    check_isoform_range(&ranges, 2, 100, 296);   // 3 bases shorter
    check_isoform_range(&ranges, 3, 647, 847);   // Starts 3 earlier, and the isoform is 3 shorter

    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in); kv_destroy(ranges);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_deletion_before_every_region_moves_them_but_keeps_the_exon_ends) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, Record_init(mock_bcf1_t(1150, "NNNNN", "N"), NULL));   // Deletes 1151 to 1154, in the first exon
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 0, 96);      // The first exon has 4 bases less
    check_isoform_range(&ranges, 1, 96, 146);    // Still 50 bases of the second exon, which now start at 96
    check_isoform_range(&ranges, 2, 96, 295);
    check_isoform_range(&ranges, 3, 646, 846);

    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in); kv_destroy(ranges);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__an_insertion_inside_a_region_lengthens_it) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, Record_init(mock_bcf1_t(2020, "N", "Nxy"), NULL));   // Inserts 2 bases between 2020 and 2021
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 0, 100);
    check_isoform_range(&ranges, 1, 100, 152);   // The inserted bases lie in the start-proximal stretch
    check_isoform_range(&ranges, 2, 100, 301);
    check_isoform_range(&ranges, 3, 652, 852);

    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in); kv_destroy(ranges);
    nmd_layout_destroy(&layout);
} END_TEST

// An insertion at a border belongs to the side that holds the base it is attached to: the one before it in the direction of the
// transcript, so the last base of a region takes the insertion after it in that direction, and the first base does not take the
// one before it. The reverse layout is the same transcript on the other strand, so it gives the same ranges with the records
// attached to the other neighbour (the lower base of the two it lies between).
START_TEST(check_nmd_escape_ranges_in_isoform__reverse_an_insertion_at_the_end_of_a_region_is_inside_and_at_the_start_outside) {
    NmdLayout layout = nmd_layout(true, 0);

    // The second start-proximal stretch is the coding bases 100 to 149 of the transcript, the genomic [2550, 2600). Coding base
    // 149 is the genomic base 2550, and the insertion between the coding bases 149 and 150 is between the genomic 2549 and 2550:
    // the VCF record is attached to 2549, and in the direction of the transcript to 2550, which is inside the stretch
    Records at_end; kv_init(at_end);
    kv_push(Record *, at_end, Record_init(mock_bcf1_t(2549, "N", "Nxy"), NULL));
    Ranges ranges;
    nmd_ranges(&layout, &at_end, &ranges);
    check_isoform_range(&ranges, 0, 100, 152);
    Record_destroy(kv_A(at_end, 0)); kv_destroy(at_end); kv_destroy(ranges);

    // The coding base 100 is the genomic base 2599, and the insertion between the coding bases 99 and 100 is between the
    // genomic 2599 and 2600: attached to 2599 in the VCF, and in the direction of the transcript to 2600, the last base of the
    // intron. Outside the stretch, and in the isoform at the start of the exon (an insertion at the end of an intron goes to the exon)
    Records at_start; kv_init(at_start);
    kv_push(Record *, at_start, Record_init(mock_bcf1_t(2599, "N", "Nxy"), NULL));
    nmd_ranges(&layout, &at_start, &ranges);
    check_isoform_range(&ranges, 0, 102, 152);
    Record_destroy(kv_A(at_start, 0)); kv_destroy(at_start); kv_destroy(ranges);

    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__an_insertion_at_the_end_of_a_region_is_inside_and_at_the_start_outside) {
    NmdLayout layout = nmd_layout(false, 0);

    // Attached to 2049: the last base of the start-proximal stretch of the second exon, [2000, 2050)
    Records at_end; kv_init(at_end);
    kv_push(Record *, at_end, Record_init(mock_bcf1_t(2049, "N", "Nxy"), NULL));
    Ranges ranges;
    nmd_ranges(&layout, &at_end, &ranges);
    check_isoform_range(&ranges, 1, 100, 152);
    Record_destroy(kv_A(at_end, 0)); kv_destroy(at_end); kv_destroy(ranges);

    // Attached to 1999, the last base of the intron before the stretch: the inserted bases are not in the stretch. They are in
    // the isoform, at the start of the exon (an insertion at the end of an intron goes to the exon), so the stretch starts after them
    Records at_start; kv_init(at_start);
    kv_push(Record *, at_start, Record_init(mock_bcf1_t(1999, "N", "Nxy"), NULL));
    nmd_ranges(&layout, &at_start, &ranges);
    check_isoform_range(&ranges, 1, 102, 152);
    Record_destroy(kv_A(at_start, 0)); kv_destroy(at_start); kv_destroy(ranges);

    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_retained_last_intron_is_distal) {
    NmdLayout layout = nmd_layout(false, 2);
    Records hap_in; kv_init(hap_in);
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // The isoform is 100, then 600, then the 400 bases of the intron at 700 to 1099, then 150 at 1100 to 1249
    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 3, 650, 1250);   // The end of the second exon, the intron and the last exon

    ck_assert(position_in_ranges(&ranges, 700));    // In the retained intron
    ck_assert(position_in_ranges(&ranges, 1099));
    ck_assert(position_in_ranges(&ranges, 1100));
    ck_assert(!position_in_ranges(&ranges, 649));   // Before the last 50 bases of the second exon

    kv_destroy(ranges); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_retained_intron_elsewhere_is_not_an_escape_region) {
    NmdLayout layout = nmd_layout(false, 1);
    Records hap_in; kv_init(hap_in);
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // The isoform is 100 bases, then the 800 bases of the first intron at 100 to 899, then 600 at 900 to 1499, then 150
    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 0, 0, 100);
    check_isoform_range(&ranges, 1, 900, 950);    // Start-proximal is the coding bases, not the intron before them
    check_isoform_range(&ranges, 2, 900, 1099);
    check_isoform_range(&ranges, 3, 1450, 1650);

    ck_assert(!position_in_ranges(&ranges, 100));   // The retained intron: not an escape region
    ck_assert(!position_in_ranges(&ranges, 500));
    ck_assert(!position_in_ranges(&ranges, 899));
    ck_assert(position_in_ranges(&ranges, 900));

    kv_destroy(ranges); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__reverse_with_a_retained_last_intron_is_distal) {
    NmdLayout layout = nmd_layout(true, 2);
    Records hap_in; kv_init(hap_in);
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // In transcript order: 100 bases, 600, then the 700 bases of the retained intron at 700 to 1399, then 150 at 1400 to 1549
    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 3, 650, 1550);
    ck_assert(position_in_ranges(&ranges, 700));
    ck_assert(position_in_ranges(&ranges, 1399));
    ck_assert(position_in_ranges(&ranges, 1400));
    ck_assert(!position_in_ranges(&ranges, 649));

    kv_destroy(ranges); kv_destroy(hap_in);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_region_that_keeps_one_base_is_one_base_long) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    char ref[52]; memset(ref, 'N', 51); ref[51] = '\0';
    kv_push(Record *, hap_in, Record_init(mock_bcf1_t(2000, ref, "N"), NULL));   // Deletes 2001 to 2050
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // 2001 to 2049 are gone from the start-proximal stretch [2000, 2050): the base 2000 is left
    ck_assert_int_eq(kv_size(ranges), 4);
    check_isoform_range(&ranges, 1, 100, 101);

    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in); kv_destroy(ranges);
    nmd_layout_destroy(&layout);
} END_TEST

START_TEST(check_nmd_escape_ranges_in_isoform__a_region_without_bases_in_the_isoform_is_left_out) {
    NmdLayout layout = nmd_layout(false, 0);
    Records hap_in; kv_init(hap_in);
    char ref[52]; memset(ref, 'N', 51); ref[51] = '\0';
    kv_push(Record *, hap_in, Record_init(mock_bcf1_t(1999, ref, "N"), NULL));   // Deletes 2000 to 2049, all of the stretch
    Ranges ranges;
    nmd_ranges(&layout, &hap_in, &ranges);

    // The second start-proximal stretch is gone. The long exon region lost the same 50 bases, and the distal one starts 50 earlier
    ck_assert_int_eq(kv_size(ranges), 3);
    check_isoform_range(&ranges, 0, 0, 100);
    check_isoform_range(&ranges, 1, 100, 249);
    check_isoform_range(&ranges, 2, 600, 800);

    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in); kv_destroy(ranges);
    nmd_layout_destroy(&layout);
} END_TEST

/*
* PTV stop location. The reference coding sequences here have no stop before their own, so after a change the first stop in the
* reading frame of the new sequence (not counting the last codon, the original stop) is the new stop. The oracle is a plain loop
* over the codons of the changed sequence.
*/
static bool test_is_stop(const char *codon) {
    return codon[0] == 'T' && ((codon[1] == 'A' && (codon[2] == 'A' || codon[2] == 'G')) || (codon[1] == 'G' && codon[2] == 'A'));
}

// The first stop in frame from the start of the sequence that is not the last codon
static int64_t oracle_first_stop(const char *seq, const int64_t length) {
    for (int64_t i = 0; i + 3 <= length; i += 3) {
        if (i + 3 == length) break;
        if (test_is_stop(seq + i)) return i;
    }

    return -1;
}

// The first such stop among the codons that overlap the changed bases [position, position + changed)
static int64_t oracle_window_stop(const char *seq, const int64_t length, const int64_t position, int64_t changed) {
    if (changed < 1) changed = 1;

    for (int64_t i = 0; i + 3 <= length; i += 3) {
        if (i + 3 == length) break;
        if (i + 3 <= position) continue;
        if (i >= position + changed) break;
        if (test_is_stop(seq + i)) return i;
    }

    return -1;
}

static uint32_t lcg_state;
static uint32_t lcg_next(void) {
    lcg_state = lcg_state * 1664525u + 1013904223u;
    return lcg_state >> 8;
}

// n_codons codons: ATG, sense codons and TAA at the end
static char *random_cds(const int n_codons, const uint32_t seed) {
    static const char bases[] = "ACGT";
    lcg_state = seed;

    char *cds = malloc(3 * n_codons + 1);
    for (int c = 0; c < n_codons; c++) {
        char codon[4] = "ATG";
        if (c == n_codons - 1) memcpy(codon, "TAA", 4);
        else if (c > 0) {
            do {
                const uint32_t draw = lcg_next();
                codon[0] = bases[draw & 3]; codon[1] = bases[(draw >> 2) & 3]; codon[2] = bases[(draw >> 4) & 3];
            } while (test_is_stop(codon));
        }

        memcpy(cds + 3 * c, codon, 3);
    }
    cds[3 * n_codons] = '\0';

    return cds;
}

enum TestVariantKind { TV_SNV, TV_MNV, TV_DEL, TV_INS };

// A change to a coding sequence. SNV and MNV: bases replace the bases from index. DEL: length bases from index are deleted.
// INS: bases are inserted after the base at index, which is the base the insertion is attached to.
typedef struct TestVariant {
    enum TestVariantKind kind;
    int64_t index;
    int64_t length;
    const char *bases;
} TestVariant;

// The changed sequence, and where the change begins (position) and how many bases from there are changed (changed_length): the
// bases of an SNV or MNV, none for a deletion (the junction), and for an insertion the attached base and the inserted ones
static char *apply_variant(const char *ref, const TestVariant v, int64_t *position, int64_t *changed_length) {
    const int64_t n = (int64_t) strlen(ref);
    const int64_t added = v.kind == TV_INS ? (int64_t) strlen(v.bases) : 0;
    char *m = malloc(n + added + 1);

    switch (v.kind) {
        case TV_SNV:
        case TV_MNV: {
            const int64_t k = (int64_t) strlen(v.bases);
            memcpy(m, ref, n + 1);
            memcpy(m + v.index, v.bases, k);
            *position = v.index;
            *changed_length = k;
            break;
        }
        case TV_DEL:
            memcpy(m, ref, v.index);
            memcpy(m + v.index, ref + v.index + v.length, n - v.index - v.length + 1);
            *position = v.index;
            *changed_length = 0;
            break;
        case TV_INS:
            memcpy(m, ref, v.index + 1);
            memcpy(m + v.index + 1, v.bases, added);
            memcpy(m + v.index + 1 + added, ref + v.index + 1, n - v.index);
            *position = v.index;
            *changed_length = 1 + added;
            break;
    }

    return m;
}

static void expect_stops(const char *ref, const TestVariant v, const int64_t expected_in_codons, const int64_t expected_downstream) {
    int64_t position, changed;
    char *m = apply_variant(ref, v, &position, &changed);
    const int64_t length = (int64_t) strlen(m);

    const int64_t in_codons = ptv_stop_in_codons(m, length, position, changed);
    const int64_t downstream = ptv_stop_downstream(m, length, position);

    ck_assert_msg(in_codons == expected_in_codons && downstream == expected_downstream,
                  "%s: in codons %ld (expected %ld), downstream %ld (expected %ld)", m,
                  (long) in_codons, (long) expected_in_codons, (long) downstream, (long) expected_downstream);
    free(m);
}

// ATG GCC AAA CTG GGC TTT CGA TAA
#define PTV_REF "ATGGCCAAACTGGGCTTTCGATAA"

START_TEST(check_ptv_stop__a_substitution_that_makes_a_stop_gives_the_index_of_its_codon) {
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 6, 1, "T" }, 6, 6);   // AAA to TAA

    // Each position of a codon: TGG to TGA, TAC to TAA and TAG, TTA to TAA
    expect_stops("ATGTGGAAACTGGGCTTTCGATAA", (TestVariant) { TV_SNV, 5, 1, "A" }, 3, 3);
    expect_stops("ATGTACAAACTGGGCTTTCGATAA", (TestVariant) { TV_SNV, 5, 1, "A" }, 3, 3);
    expect_stops("ATGTACAAACTGGGCTTTCGATAA", (TestVariant) { TV_SNV, 5, 1, "G" }, 3, 3);
    expect_stops("ATGTTAAAACTGGGCTTTCGATAA", (TestVariant) { TV_SNV, 4, 1, "A" }, 3, 3);
} END_TEST

START_TEST(check_ptv_stop__a_substitution_without_a_stop_has_none) {
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 6, 1, "G" }, -1, -1);
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 0, 1, "C" }, -1, -1);   // The start codon is not a stop either
} END_TEST

START_TEST(check_ptv_stop__the_original_stop_codon_never_counts) {
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 21, 1, "C" }, -1, -1);   // TAA to CAA: the stop is lost
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 23, 1, "G" }, -1, -1);   // TAA to TAG: still a stop, but the original one
    expect_stops(PTV_REF, (TestVariant) { TV_SNV, 22, 1, "G" }, -1, -1);   // TAA to TGA
} END_TEST

START_TEST(check_ptv_stop__a_substitution_that_runs_into_the_next_codon_is_not_cut_short) {
    // ATG GCC CAA AAA CTG TAA, with the bases 5 and 6 replaced: GCG TAA
    const char *ref = "ATGGCCCAAAAACTGTAA";
    const TestVariant v = { TV_MNV, 5, 2, "GT" };
    expect_stops(ref, v, 6, 6);

    // A window of the length of the change from the start of the codon of the first base reads only the codon GCG
    int64_t position, changed;
    char *m = apply_variant(ref, v, &position, &changed);
    ck_assert_int_eq(first_stop_codon(m, 18, codon_start(position), changed, true), -1);
    free(m);
} END_TEST

START_TEST(check_ptv_stop__an_insertion_in_frame_that_is_a_stop_gives_its_index) {
    expect_stops(PTV_REF, (TestVariant) { TV_INS, 5, 3, "TGA" }, 6, 6);   // ATG GCC TGA AAA ...
    expect_stops(PTV_REF, (TestVariant) { TV_INS, 5, 3, "TAG" }, 6, 6);
    expect_stops(PTV_REF, (TestVariant) { TV_INS, 5, 3, "GGG" }, -1, -1);
} END_TEST

START_TEST(check_ptv_stop__a_frameshifting_insertion_makes_a_stop_with_the_bases_after_it) {
    expect_stops(PTV_REF, (TestVariant) { TV_INS, 5, 1, "T" }, 6, 6);     // ATG GCC TAA ACT ...
    expect_stops(PTV_REF, (TestVariant) { TV_INS, 4, 2, "TG" }, -1, -1);  // ATG GCT GCA AAC TGG GCT TTC GAT AA: no stop in the new frame
} END_TEST

START_TEST(check_ptv_stop__a_deletion_makes_a_stop_when_the_bases_after_it_do) {
    expect_stops("ATGCTAAAAGGGTAA", (TestVariant) { TV_DEL, 3, 1, NULL }, 3, 3);   // ATG TAA AAG GGT AA
} END_TEST

// The stop can be some codons after the change: only the scan to the end of the sequence finds it
START_TEST(check_ptv_stop__a_frameshift_finds_a_stop_further_on_that_the_codons_of_the_change_do_not_have) {
    // ATG GCC CTA AGG ... without the G at 3: ATG CCC TAA GGA AAA AAC TGT AA
    expect_stops("ATGGCCCTAAGGAAAAAACTGTAA", (TestVariant) { TV_DEL, 3, 1, NULL }, -1, 6);
} END_TEST

START_TEST(check_ptv_stop__a_frameshift_that_reads_on_without_a_stop_has_none) {
    expect_stops(PTV_REF, (TestVariant) { TV_DEL, 4, 1, NULL }, -1, -1);   // ATG GCA AAC TGG GCT TTC GAT AA
} END_TEST

START_TEST(check_ptv_stop__a_deletion_in_frame_has_no_new_stop) {
    expect_stops(PTV_REF, (TestVariant) { TV_DEL, 6, 3, NULL }, -1, -1);    // ATG GCC CTG GGC TTT CGA TAA: only the original one
    expect_stops(PTV_REF, (TestVariant) { TV_DEL, 3, 6, NULL }, -1, -1);
} END_TEST

START_TEST(check_ptv_stop__a_position_outside_the_sequence_is_moved_into_it) {
    const char *seq = "ATGTAAGGGTAA";

    ck_assert_int_eq(ptv_stop_downstream(seq, 12, 100), -1);   // The last codon: the original stop
    ck_assert_int_eq(ptv_stop_downstream(seq, 12, 12), -1);    // A deletion of the last bases puts the position at the length
    ck_assert_int_eq(ptv_stop_downstream(seq, 12, -5), 3);     // Moved to the start: the scan finds the first stop
    ck_assert_int_eq(ptv_stop_downstream(seq, 12, 0), 3);
    ck_assert_int_eq(ptv_stop_downstream(seq, 12, 4), 3);      // Inside the codon of the stop
    ck_assert_int_eq(ptv_stop_downstream(seq, 12, 6), -1);     // Past it
    ck_assert_int_eq(ptv_stop_in_codons(seq, 12, 100, 1), -1);
    ck_assert_int_eq(ptv_stop_in_codons(seq, 12, -5, 1), -1);  // Moved to base 0: only the codon ATG
    ck_assert_int_eq(ptv_stop_in_codons(seq, 12, 3, 1), 3);
    ck_assert_int_eq(ptv_stop_in_codons(seq, 12, 3, 0), 3);    // At least one base is looked at
} END_TEST

START_TEST(check_ptv_stop__an_empty_sequence_has_no_stop) {
    ck_assert_int_eq(ptv_stop_downstream("", 0, 0), -1);
    ck_assert_int_eq(ptv_stop_in_codons("", 0, 0, 3), -1);
    ck_assert_int_eq(ptv_stop_downstream(NULL, 0, 5), -1);
    ck_assert_int_eq(ptv_stop_in_codons(NULL, 0, 5, 3), -1);
} END_TEST

START_TEST(check_nmd_flag_for_stop__names_the_region_of_the_stop_and_nothing_outside_them_or_without_a_stop) {
    NmdRegions regions;
    kv_init(regions.start_proximal); kv_init(regions.long_exon); kv_init(regions.distal);
    kv_push(Range, regions.start_proximal, ((Range) { 10, 20 }));
    kv_push(Range, regions.long_exon, ((Range) { 50, 60 }));
    kv_push(Range, regions.distal, ((Range) { 100, 110 }));

    ck_assert_int_eq(nmd_flag_for_stop(-1, &regions), 0);
    ck_assert_int_eq(nmd_flag_for_stop(9, &regions), 0);
    ck_assert_int_eq(nmd_flag_for_stop(10, &regions), NMD_START_PROXIMAL);
    ck_assert_int_eq(nmd_flag_for_stop(19, &regions), NMD_START_PROXIMAL);
    ck_assert_int_eq(nmd_flag_for_stop(20, &regions), 0);
    ck_assert_int_eq(nmd_flag_for_stop(50, &regions), NMD_LONG_EXON);
    ck_assert_int_eq(nmd_flag_for_stop(60, &regions), 0);
    ck_assert_int_eq(nmd_flag_for_stop(100, &regions), NMD_DISTAL);
    ck_assert_int_eq(nmd_flag_for_stop(110, &regions), 0);
    NmdRegions_destroy(&regions);

    NmdRegions none;
    kv_init(none.start_proximal); kv_init(none.long_exon); kv_init(none.distal);
    ck_assert_int_eq(nmd_flag_for_stop(5, &none), 0);
    ck_assert_int_eq(nmd_flag_for_stop(-1, &none), 0);
    NmdRegions_destroy(&none);
} END_TEST

// The start-proximal stretch can overlap a long exon or the distal region: it comes first, and then the long exon
START_TEST(check_nmd_flag_for_stop__overlapping_regions_give_one_flag_in_a_fixed_order) {
    NmdRegions regions;
    kv_init(regions.start_proximal); kv_init(regions.long_exon); kv_init(regions.distal);
    kv_push(Range, regions.start_proximal, ((Range) { 0, 150 }));
    kv_push(Range, regions.long_exon, ((Range) { 100, 300 }));
    kv_push(Range, regions.distal, ((Range) { 120, 400 }));

    ck_assert_int_eq(nmd_flag_for_stop(110, &regions), NMD_START_PROXIMAL);
    ck_assert_int_eq(nmd_flag_for_stop(130, &regions), NMD_START_PROXIMAL);
    ck_assert_int_eq(nmd_flag_for_stop(200, &regions), NMD_LONG_EXON);
    ck_assert_int_eq(nmd_flag_for_stop(350, &regions), NMD_DISTAL);
    NmdRegions_destroy(&regions);
} END_TEST

START_TEST(check_nmd_flag_for_stop__the_flags_are_distinct_bits) {
    ck_assert_int_ne(NMD_START_PROXIMAL, NMD_LONG_EXON);
    ck_assert_int_ne(NMD_START_PROXIMAL, NMD_DISTAL);
    ck_assert_int_ne(NMD_LONG_EXON, NMD_DISTAL);
    ck_assert_int_eq((NMD_START_PROXIMAL | NMD_LONG_EXON | NMD_DISTAL) & (UNPHASED_VARIANT_IN_PREDICTION | INFRAME | NON_CONSTITUTIVE | UPSTREAM_PTV), 0);
} END_TEST

static void check_variant_against_translation(const char *ref, const TestVariant v) {
    int64_t position, changed;
    char *m = apply_variant(ref, v, &position, &changed);
    const int64_t length = (int64_t) strlen(m);

    const int64_t downstream = ptv_stop_downstream(m, length, position);
    const int64_t in_codons = ptv_stop_in_codons(m, length, position, changed);
    const int64_t expected_downstream = oracle_first_stop(m, length);
    const int64_t expected_in_codons = oracle_window_stop(m, length, position, changed);

    ck_assert_msg(downstream == expected_downstream && in_codons == expected_in_codons && (in_codons < 0 || in_codons == downstream),
                  "%s -> %s (kind %d at %ld, length %ld): downstream %ld (expected %ld), in codons %ld (expected %ld)",
                  ref, m, (int) v.kind, (long) v.index, (long) v.length, (long) downstream, (long) expected_downstream,
                  (long) in_codons, (long) expected_in_codons);
    free(m);
}

// Every substitution, MNV, deletion of one to four bases and insertion on 24 random coding sequences
START_TEST(check_ptv_stop__agrees_with_a_translation_of_every_small_variant) {
    static const char *snv_bases[] = { "A", "C", "G", "T" };
    static const char *mnv_bases[] = { "TA", "TG", "AA", "CT", "TAA", "TGA", "TAG", "GGG" };
    static const char *inserted[] = { "T", "A", "TA", "TG", "TAA", "TGA", "TAG", "GCA", "TTTT", "TGAT" };

    for (uint32_t seed = 1; seed <= 24; seed++) {
        char *ref = random_cds(3 + (int) (seed % 9), seed * 7919u);
        const int64_t n = (int64_t) strlen(ref);

        for (int64_t index = 0; index < n; index++) {
            for (int b = 0; b < 4; b++) check_variant_against_translation(ref, (TestVariant) { TV_SNV, index, 1, snv_bases[b] });
        }

        for (int b = 0; b < 8; b++) {
            const int64_t k = (int64_t) strlen(mnv_bases[b]);
            for (int64_t index = 0; index + k <= n; index++) check_variant_against_translation(ref, (TestVariant) { TV_MNV, index, k, mnv_bases[b] });
        }

        for (int64_t length = 1; length <= 4; length++) {
            for (int64_t index = 0; index + length <= n; index++) check_variant_against_translation(ref, (TestVariant) { TV_DEL, index, length, NULL });
        }

        // Attached to the last base of the stop codon an insertion would lie outside the coding sequence
        for (int64_t index = 0; index <= n - 2; index++) {
            for (int b = 0; b < 10; b++) check_variant_against_translation(ref, (TestVariant) { TV_INS, index, (int64_t) strlen(inserted[b]), inserted[b] });
        }

        free(ref);
    }
} END_TEST

/*
* The whole way from a change to its flag, on a gene of three exons whose coding sequence has 651 bases (217 codons). The change is
* given on the transcript, applied to the coding sequence by the oracle, and also written as a record on the gene (the reverse
* strand is the same transcript on the other strand, so the record has the complemented bases and the order turned around). The
* gene sequence goes through get_haplotype_sequence and splice_haplotype_sequence, the stop is found in the result, and the
* escape regions of the transcript go through the conversions. The oracle finds the stop in the changed coding sequence and tests
* the reference base it came from against the escape regions of the reference: the coding bases 0 to 149 (start-proximal) and 450
* to 650 (the last 50 bases of the second exon, 450 to 499, and all of the last exon).
*
* Forward: exons at gene offsets [0, 300), [400, 600), [700, 851) (so coding bases 0-299, 300-499 and 500-650).
* Reverse: exons at [551, 851), [251, 451), [0, 151); coding base t is at offset 850 - t, 450 - (t - 300) and 150 - (t - 500).
*/
#define NMD_CDS_CODONS 217
#define NMD_CDS_LENGTH 651
#define NMD_GENE_LENGTH 851
#define NMD_GENE_START 1000

typedef struct NmdGene {
    bool reverse;
    char *ref_cds;               // The coding sequence in the direction of the transcript
    kstring_t seq;               // The gene, ascending
    Range gene;
    Transcript *transcript;
    NonCodingRegions ncrs;       // In transcript order, with an empty UTR at each end
} NmdGene;

static char test_complement(const char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return base;
    }
}

// The offset from the gene start of the base at index t of the coding sequence
static int64_t nmd_offset(const bool reverse, const int64_t t) {
    if (!reverse) return t < 300 ? t : t < 500 ? 400 + (t - 300) : 700 + (t - 500);
    return t < 300 ? 850 - t : t < 500 ? 450 - (t - 300) : 150 - (t - 500);
}

static NmdGene nmd_gene(const bool reverse) {
    // Codons that a substitution turns into a stop (CAA to TAA) at both sides of the borders of the regions, and a run that has
    // no stop in the frame after the base at 639 is deleted
    static const int planted_stops[] = { 20, 49, 50, 70, 120, 149, 150, 160, 200 };

    NmdGene g = {0};
    g.reverse = reverse;
    g.ref_cds = random_cds(NMD_CDS_CODONS, 4242);
    for (size_t i = 0; i < sizeof(planted_stops) / sizeof(planted_stops[0]); i++) memcpy(g.ref_cds + 3 * planted_stops[i], "CAA", 3);
    memcpy(g.ref_cds + 3 * 213, "GCCGCCGCC", 9);

    // The introns have bases of their own that nothing reads
    lcg_state = 99;
    char *seq = malloc(NMD_GENE_LENGTH + 1);
    for (int i = 0; i < NMD_GENE_LENGTH; i++) seq[i] = "ACGT"[lcg_next() & 3];
    seq[NMD_GENE_LENGTH] = '\0';
    for (int64_t t = 0; t < NMD_CDS_LENGTH; t++) seq[nmd_offset(reverse, t)] = reverse ? test_complement(g.ref_cds[t]) : g.ref_cds[t];
    g.seq = (kstring_t) { NMD_GENE_LENGTH, NMD_GENE_LENGTH + 1, seq };
    g.gene = (Range) { NMD_GENE_START, NMD_GENE_START + NMD_GENE_LENGTH };

    kv_init(g.ncrs);
    if (!reverse) {
        Range cds[] = { {1000, 1299}, {1400, 1599}, {1700, 1850} };
        int constitutive[] = { 1, 1, 1 };
        TranscriptData data = { .id = 1, .start = 1000, .end = 1850, .strand = STRAND_FWD,
                                .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = cds, .nexons = 3 };
        g.transcript = Transcript_build(&data);

        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1000, 1000 }, UTR5));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1300, 1400 }, 0));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1600, 1700 }, 0));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1851, 1851 }, UTR3));
    } else {
        Range cds[] = { {1000, 1150}, {1251, 1450}, {1551, 1850} };
        int constitutive[] = { 1, 1, 1 };
        TranscriptData data = { .id = 2, .start = 1000, .end = 1850, .strand = STRAND_REV,
                                .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = cds, .nexons = 3 };
        g.transcript = Transcript_build(&data);

        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1851, 1851 }, UTR5));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1451, 1551 }, 0));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1151, 1251 }, 0));
        kv_push(NonCodingRegion, g.ncrs, region(1, (Range) { 1000, 1000 }, UTR3));
    }

    return g;
}

static void nmd_gene_destroy(NmdGene *g) {
    free(g->ref_cds);
    free(g->seq.s);
    Transcript_destroy(g->transcript);
    kv_destroy(g->ncrs);
}

// The change as a record on the gene. ref_buf and alt_buf are the alleles, which the record keeps pointing at: the caller keeps
// them alive for as long as the record is used (both 16 bases or more)
static Record *nmd_record(const NmdGene *g, const TestVariant v, char *ref_buf, char *alt_buf) {
    const char *seq = g->seq.s;
    int64_t anchor = 0;         // Offset of the first base of REF
    int64_t ref_length = 1;

    if (v.kind == TV_SNV) {
        anchor = nmd_offset(g->reverse, v.index);
        alt_buf[0] = g->reverse ? test_complement(v.bases[0]) : v.bases[0];
        alt_buf[1] = '\0';
    } else if (v.kind == TV_DEL) {
        // The anchor is the base before the deleted ones, in ascending order
        anchor = (g->reverse ? nmd_offset(true, v.index + v.length - 1) : nmd_offset(false, v.index)) - 1;
        ref_length = v.length + 1;
        alt_buf[0] = seq[anchor];
        alt_buf[1] = '\0';
    } else {
        // The anchor is the lower of the two bases the insertion lies between: the base before it on the forward strand, and
        // the base after it in the direction of the transcript on the reverse strand
        anchor = g->reverse ? nmd_offset(true, v.index + 1) : nmd_offset(false, v.index);
        const size_t added = strlen(v.bases);
        alt_buf[0] = seq[anchor];
        for (size_t i = 0; i < added; i++) {
            alt_buf[1 + i] = g->reverse ? test_complement(v.bases[added - 1 - i]) : v.bases[i];
        }
        alt_buf[1 + added] = '\0';
    }

    memcpy(ref_buf, seq + anchor, ref_length);
    ref_buf[ref_length] = '\0';

    return Record_init(mock_bcf1_t(NMD_GENE_START + anchor, ref_buf, alt_buf), NULL);
}

// What evaluate_isoform does with a gene and its haplotype: the sequence of the isoform in the direction of the transcript, the
// stop found from the position, and the flag of the stop against the escape regions
static void nmd_pipeline(const NmdGene *g, const Records *hap_in, const int64_t position, kstring_t *spliced, int64_t *stop, uint32_t *flag) {
    kstring_t hap_seq = {0};
    get_haplotype_sequence(g->seq, hap_in, g->gene, &hap_seq);

    NonCodingRegions hap_ncrs; kv_init(hap_ncrs);
    for (size_t k = 0; k < kv_size(g->ncrs); k++) {
        NonCodingRegion ncr = kv_A(g->ncrs, k);
        ncr.coords = region_to_haplotype(hap_in, g->gene, ncr.coords, ncr.is_utr != 0);
        kv_push(NonCodingRegion, hap_ncrs, ncr);
    }

    *spliced = (kstring_t) {0};
    splice_haplotype_sequence(hap_seq.s, &hap_ncrs, g->gene.start, g->transcript->strand, spliced);
    if (g->reverse) reverse_complement(spliced->s, spliced->l);

    *stop = ptv_stop_downstream(spliced->s, (int64_t) spliced->l, position);

    Ranges segments = isoform_segments(&hap_ncrs, g->gene.start, g->transcript->strand);
    NmdRegions regions = Transcript_nmd_regions(g->transcript);
    NmdRegions escape;
    nmd_escape_ranges_in_isoform(&regions, hap_in, g->gene, &segments, g->transcript->strand, &escape);
    *flag = nmd_flag_for_stop(*stop, &escape);

    NmdRegions_destroy(&escape);
    NmdRegions_destroy(&regions);
    kv_destroy(segments);
    kv_destroy(hap_ncrs);
    free(hap_seq.s);
}

// The base of the reference coding sequence that the base at index of the changed sequence comes from. Bases an insertion adds
// belong to the base they are attached to, which in the direction of the transcript is the one before them on both strands
// (on the reverse strand the VCF record is attached to the lower base, which is the one after them in that direction)
static int64_t oracle_origin(const TestVariant v, const int64_t index) {
    switch (v.kind) {
        case TV_INS: {
            const int64_t added = (int64_t) strlen(v.bases);
            if (index <= v.index) return index;
            if (index <= v.index + added) return v.index;
            return index - added;
        }
        case TV_DEL: return index < v.index ? index : index + v.length;
        default: return index;
    }
}

static void run_nmd_case(const NmdGene *g, const TestVariant v, const int case_number) {
    Records hap_in; kv_init(hap_in);
    char ref_buf[16], alt_buf[16];
    kv_push(Record *, hap_in, nmd_record(g, v, ref_buf, alt_buf));

    int64_t position, changed;
    char *expected_seq = apply_variant(g->ref_cds, v, &position, &changed);
    const int64_t expected_length = (int64_t) strlen(expected_seq);

    kstring_t spliced; int64_t stop; uint32_t flag;
    nmd_pipeline(g, &hap_in, position, &spliced, &stop, &flag);

    ck_assert_msg((int64_t) spliced.l == expected_length && memcmp(spliced.s, expected_seq, expected_length) == 0,
                  "case %d (kind %d at %ld, reverse %d): the isoform sequence is not the changed coding sequence", case_number, (int) v.kind, (long) v.index, (int) g->reverse);

    const int64_t expected_stop = oracle_first_stop(expected_seq, expected_length);
    ck_assert_msg(stop == expected_stop, "case %d (kind %d at %ld, reverse %d): stop %ld, expected %ld", case_number, (int) v.kind, (long) v.index,
                  (int) g->reverse, (long) stop, (long) expected_stop);

    // No exon is long enough for a long exon region, so the stop is start-proximal, distal or in neither
    uint32_t expected_flag = 0;
    if (expected_stop >= 0) {
        const int64_t origin = oracle_origin(v, expected_stop);
        expected_flag = origin < 150 ? NMD_START_PROXIMAL : origin >= 450 ? NMD_DISTAL : 0;
    }
    ck_assert_msg(flag == expected_flag, "case %d (kind %d at %ld, reverse %d): stop %ld gives flag %u, expected %u", case_number, (int) v.kind,
                  (long) v.index, (int) g->reverse, (long) stop, flag, expected_flag);

    free(spliced.s); free(expected_seq);
    Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in);
}

static const TestVariant NMD_CASES[] = {
    // A substitution turns a planted CAA into a stop: 60 and 147 are start-proximal, 150 the first base after it; 210 and 360 are
    // in the NMD target region, 447 is its last codon and 450 the first base of the last 50 of the second exon; 480 is there
    // too and 600 is in the last exon
    { TV_SNV, 60, 1, "T" }, { TV_SNV, 147, 1, "T" }, { TV_SNV, 150, 1, "T" }, { TV_SNV, 210, 1, "T" }, { TV_SNV, 360, 1, "T" },
    { TV_SNV, 447, 1, "T" }, { TV_SNV, 450, 1, "T" }, { TV_SNV, 480, 1, "T" }, { TV_SNV, 600, 1, "T" },
    // One base deleted: the frame shifts and the stop is wherever the sequence puts it
    { TV_DEL, 10, 1, NULL }, { TV_DEL, 40, 1, NULL }, { TV_DEL, 100, 1, NULL }, { TV_DEL, 160, 1, NULL }, { TV_DEL, 250, 1, NULL },
    { TV_DEL, 310, 1, NULL }, { TV_DEL, 420, 1, NULL }, { TV_DEL, 470, 1, NULL }, { TV_DEL, 520, 1, NULL }, { TV_DEL, 600, 1, NULL },
    // Two bases, and three: in frame
    { TV_DEL, 120, 2, NULL }, { TV_DEL, 330, 2, NULL }, { TV_DEL, 540, 2, NULL }, { TV_DEL, 200, 3, NULL }, { TV_DEL, 400, 3, NULL },
    // Insertions that shift the frame
    { TV_INS, 5, 1, "T" }, { TV_INS, 70, 1, "T" }, { TV_INS, 190, 2, "TG" }, { TV_INS, 305, 1, "A" }, { TV_INS, 440, 2, "TA" },
    { TV_INS, 470, 1, "T" }, { TV_INS, 560, 2, "GC" }, { TV_INS, 640, 1, "T" },
    // Insertions of a stop codon, in frame
    { TV_INS, 70, 3, "TGA" }, { TV_INS, 400, 3, "TAG" }, { TV_INS, 520, 3, "TAA" },
    // The base 639 starts GCC GCC GCC TAA: without it the frame has no stop before the end, and the original one is not in frame
    { TV_DEL, 639, 1, NULL },
    // A stop codon inserted at the border of a region: it belongs to the base it is attached to (before it in the transcript)
    { TV_INS, 149, 3, "TAA" }, { TV_INS, 449, 3, "TAA" },
};

static void run_all_nmd_cases(const bool reverse) {
    NmdGene g = nmd_gene(reverse);

    // The coding sequence is what the gene says, before any change
    Records none; kv_init(none);
    kstring_t spliced; int64_t stop; uint32_t flag;
    nmd_pipeline(&g, &none, 0, &spliced, &stop, &flag);
    ck_assert_int_eq(spliced.l, NMD_CDS_LENGTH);
    ck_assert_int_eq(memcmp(spliced.s, g.ref_cds, NMD_CDS_LENGTH), 0);
    ck_assert_int_eq(stop, -1);   // Only the original stop
    free(spliced.s); kv_destroy(none);

    for (size_t i = 0; i < sizeof(NMD_CASES) / sizeof(NMD_CASES[0]); i++) run_nmd_case(&g, NMD_CASES[i], (int) i);

    nmd_gene_destroy(&g);
}

START_TEST(check_nmd_location__forward_gene_from_the_change_to_the_flag) {
    run_all_nmd_cases(false);
} END_TEST

START_TEST(check_nmd_location__reverse_gene_from_the_change_to_the_flag) {
    run_all_nmd_cases(true);
} END_TEST

// The border of a region, with the first base of the stop codon as its position. The expected flags are written out, not derived
START_TEST(check_nmd_location__a_stop_is_placed_by_its_first_base) {
    const bool strands[2] = { false, true };

    for (int si = 0; si < 2; si++) {
        NmdGene g = nmd_gene(strands[si]);
        const struct { int64_t index; uint32_t flag; } expected[] = {
            { 60, NMD_START_PROXIMAL },
            { 147, NMD_START_PROXIMAL },   // The codon 147 to 149 begins before base 150
            { 150, 0 },                    // The codon 150 to 152 begins at 150: the first NMD target base
            { 210, 0 }, { 360, 0 },
            { 447, 0 },                    // The codon 447 to 449 ends where the last 50 bases of the second exon begin
            { 450, NMD_DISTAL },           // The first base of the last 50
            { 480, NMD_DISTAL }, { 600, NMD_DISTAL },
        };

        for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); i++) {
            Records hap_in; kv_init(hap_in);
            char ref_buf[16], alt_buf[16];
            const TestVariant v = { TV_SNV, expected[i].index, 1, "T" };
            kv_push(Record *, hap_in, nmd_record(&g, v, ref_buf, alt_buf));

            kstring_t spliced; int64_t stop; uint32_t flag;
            nmd_pipeline(&g, &hap_in, expected[i].index, &spliced, &stop, &flag);

            ck_assert_msg(stop == expected[i].index && flag == expected[i].flag, "reverse %d, base %ld: stop %ld, flag %u, expected flag %u",
                          (int) g.reverse, (long) expected[i].index, (long) stop, flag, expected[i].flag);

            free(spliced.s); Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in);
        }

        nmd_gene_destroy(&g);
    }
} END_TEST

// A stop codon inserted in frame (after the last base of a codon) has its first base right after the attached one, and belongs to
// the attached base, the one before it in the direction of the transcript: attached to 149 it is start-proximal, attached to 449,
// the last base of the NMD target region, it is not. On both strands. The attached bases are the last of an exon-internal codon
START_TEST(check_nmd_location__an_inserted_stop_belongs_to_the_base_it_is_attached_to) {
    const bool strands[2] = { false, true };

    for (int si = 0; si < 2; si++) {
        NmdGene g = nmd_gene(strands[si]);
        const struct { int64_t attached_to; uint32_t flag; } expected[] = {
            { 59, NMD_START_PROXIMAL },
            { 149, NMD_START_PROXIMAL },   // The last base of the start-proximal stretch
            { 152, 0 },                    // The first codon after it
            { 209, 0 },
            { 449, 0 },                    // The last base of the NMD target region
            { 452, NMD_DISTAL },           // The first codon of the last 50 bases of the second exon
            { 479, NMD_DISTAL },
            { 599, NMD_DISTAL },           // The last exon
        };

        for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); i++) {
            ck_assert_int_eq(expected[i].attached_to % 3, 2);   // In frame

            Records hap_in; kv_init(hap_in);
            char ref_buf[16], alt_buf[16];
            const TestVariant v = { TV_INS, expected[i].attached_to, 3, "TAA" };
            kv_push(Record *, hap_in, nmd_record(&g, v, ref_buf, alt_buf));

            kstring_t spliced; int64_t stop; uint32_t flag;
            nmd_pipeline(&g, &hap_in, expected[i].attached_to, &spliced, &stop, &flag);

            ck_assert_msg(stop == expected[i].attached_to + 1 && flag == expected[i].flag,
                          "reverse %d, attached to %ld: stop %ld, flag %u, expected the stop at %ld and flag %u", (int) g.reverse,
                          (long) expected[i].attached_to, (long) stop, flag, (long) expected[i].attached_to + 1, expected[i].flag);

            free(spliced.s); Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in);
        }

        nmd_gene_destroy(&g);
    }
} END_TEST

// No stop before the end of the sequence: no flag, and the original stop codon is not found instead
START_TEST(check_nmd_location__no_new_stop_has_no_flag) {
    const bool strands[2] = { false, true };

    for (int si = 0; si < 2; si++) {
        NmdGene g = nmd_gene(strands[si]);
        Records hap_in; kv_init(hap_in);
        char ref_buf[16], alt_buf[16];
        const TestVariant v = { TV_DEL, 639, 1, NULL };
        kv_push(Record *, hap_in, nmd_record(&g, v, ref_buf, alt_buf));

        kstring_t spliced; int64_t stop; uint32_t flag;
        nmd_pipeline(&g, &hap_in, 639, &spliced, &stop, &flag);

        ck_assert_int_eq(stop, -1);
        ck_assert_int_eq(flag, 0);
        ck_assert_int_eq(spliced.l, NMD_CDS_LENGTH - 1);

        free(spliced.s); Record_destroy(kv_A(hap_in, 0)); kv_destroy(hap_in);
        nmd_gene_destroy(&g);
    }
} END_TEST

START_TEST(check_reverse_complement__complements_and_reverses) {
    char seq[] = "AACGTN";
    reverse_complement(seq, strlen(seq));

    ck_assert_str_eq(seq, "NACGTT");
} END_TEST

START_TEST(check_reverse_complement__leaves_other_characters_unchanged) {
    char seq[] = "acgT";
    reverse_complement(seq, strlen(seq));

    ck_assert_str_eq(seq, "Agca");
} END_TEST

START_TEST(check_reverse_complement__empty_sequence) {
    char seq[] = "";
    reverse_complement(seq, 0);

    ck_assert_str_eq(seq, "");
} END_TEST

static SpliceSite site_with(const float prior, const float posterior) {
    SpliceSite site = {0};
    site.prior_probability = prior;
    site.posterior_probability = posterior;
    return site;
}

START_TEST(check_isoform_prior__multiplies_the_priors_of_used_regions) {
    SpliceSite donor = site_with(0.8, 0.5);
    SpliceSite acceptor = site_with(0.6, 0.0); // Lost, counts as 1 - 0.6
    SpliceSite retained_donor = site_with(0.1, 0.0);

    NonCodingRegions ncrs; kv_init(ncrs);
    kv_push(NonCodingRegion, ncrs, region(1, (Range) { 100, 102 }, UTR5)); // No sites
    NonCodingRegion spliced = region(1, (Range) { 110, 120 }, 0);
    spliced.donor = &donor;
    spliced.acceptor = &acceptor;
    kv_push(NonCodingRegion, ncrs, spliced);
    NonCodingRegion retained = region(0, (Range) { 125, 128 }, 0); // Not used, so not counted
    retained.donor = &retained_donor;
    kv_push(NonCodingRegion, ncrs, retained);

    ck_assert_float_eq_tol(isoform_prior(&ncrs), 0.8f * 0.4f, 1e-6);

    kv_destroy(ncrs);
} END_TEST

START_TEST(check_isoform_prior__no_regions_gives_one) {
    NonCodingRegions ncrs; kv_init(ncrs);

    ck_assert_float_eq(isoform_prior(&ncrs), 1.0f);

    kv_destroy(ncrs);
} END_TEST

/*
* record_overlaps_gene and find_haplotype_records: which records make up the haplotype of a gene.
*/

// A record with the given genes (none for a record that was not annotated). Record_destroy frees the genes vector, not the genes.
static Record *record_in_genes(const int64_t pos, Gene **genes, const int n_genes, const bool unphased) {
    Record *r = Record_init(mock_bcf1_t(pos, "N", "x"), NULL);
    r->unphased = unphased;

    if (n_genes >= 0) {
        r->genes = malloc(sizeof(Genes));
        kv_init(*r->genes);
        for (int i = 0; i < n_genes; i++) kv_push(Gene *, *r->genes, genes[i]);
    }

    return r;
}

START_TEST(check_record_overlaps_gene__follows_the_changed_bases) {
    Gene gene = { .gff_id = 1, .coords = { 100, 200 } };

    Record *inside = Record_init(mock_bcf1_t(150, "N", "x"), NULL);
    Record *first = Record_init(mock_bcf1_t(100, "N", "x"), NULL);
    Record *after = Record_init(mock_bcf1_t(200, "N", "x"), NULL);
    Record *insertion_at_the_end = Record_init(mock_bcf1_t(199, "N", "Nx"), NULL); // Between 199 and 200
    Record *insertion_inside = Record_init(mock_bcf1_t(198, "N", "Nx"), NULL);
    Record *deletion_after_the_end = Record_init(mock_bcf1_t(199, "NN", "N"), NULL); // Deletes 200

    ck_assert(record_overlaps_gene(inside, &gene));
    ck_assert(record_overlaps_gene(first, &gene));
    ck_assert(!record_overlaps_gene(after, &gene));
    ck_assert(!record_overlaps_gene(insertion_at_the_end, &gene));
    ck_assert(record_overlaps_gene(insertion_inside, &gene));
    ck_assert(!record_overlaps_gene(deletion_after_the_end, &gene));

    Record_destroy(inside); Record_destroy(first); Record_destroy(after);
    Record_destroy(insertion_at_the_end); Record_destroy(insertion_inside); Record_destroy(deletion_after_the_end);
} END_TEST

START_TEST(check_find_haplotype_records__the_records_in_the_gene_in_order) {
    Gene gene = { .gff_id = 1, .coords = { 100, 200 } };
    Gene *genes[] = { &gene };
    Records phase; kv_init(phase);
    Record *a = record_in_genes(110, genes, 1, false);
    Record *b = record_in_genes(120, genes, 1, false);
    Record *c = record_in_genes(130, genes, 1, false);
    kv_push(Record *, phase, a); kv_push(Record *, phase, b); kv_push(Record *, phase, c);

    Records hap_in, hap_ex;
    find_haplotype_records(&phase, b, &gene, true, &hap_in, &hap_ex);

    ck_assert_int_eq(kv_size(hap_in), 3);
    ck_assert_ptr_eq(kv_A(hap_in, 0), a); ck_assert_ptr_eq(kv_A(hap_in, 1), b); ck_assert_ptr_eq(kv_A(hap_in, 2), c);
    ck_assert_int_eq(kv_size(hap_ex), 2);
    ck_assert_ptr_eq(kv_A(hap_ex, 0), a); ck_assert_ptr_eq(kv_A(hap_ex, 1), c);

    kv_destroy(hap_in); kv_destroy(hap_ex); kv_destroy(phase);
    Record_destroy(a); Record_destroy(b); Record_destroy(c);
} END_TEST

START_TEST(check_find_haplotype_records__other_genes_other_places_and_unannotated_records_are_left_out) {
    Gene gene = { .gff_id = 1, .coords = { 100, 200 } };
    Gene other = { .gff_id = 2, .coords = { 100, 200 } };
    Gene *genes[] = { &gene };
    Gene *other_genes[] = { &other };
    Gene *both[] = { &other, &gene };
    Records phase; kv_init(phase);
    Record *focal = record_in_genes(150, genes, 1, false);
    Record *outside = record_in_genes(250, genes, 1, false);          // Not in the gene
    Record *elsewhere = record_in_genes(160, other_genes, 1, false);  // In another gene
    Record *unannotated = record_in_genes(170, NULL, -1, false);      // No genes at all
    Record *listed_with_another = record_in_genes(180, both, 2, false);
    kv_push(Record *, phase, focal); kv_push(Record *, phase, outside); kv_push(Record *, phase, elsewhere);
    kv_push(Record *, phase, unannotated); kv_push(Record *, phase, listed_with_another);

    Records hap_in;
    find_haplotype_records(&phase, focal, &gene, true, &hap_in, NULL);

    ck_assert_int_eq(kv_size(hap_in), 2);
    ck_assert_ptr_eq(kv_A(hap_in, 0), focal);
    ck_assert_ptr_eq(kv_A(hap_in, 1), listed_with_another);

    kv_destroy(hap_in); kv_destroy(phase);
    Record_destroy(focal); Record_destroy(outside); Record_destroy(elsewhere);
    Record_destroy(unannotated); Record_destroy(listed_with_another);
} END_TEST

// An unphased record is on the haplotype of itself, and of another record only when asked for
START_TEST(check_find_haplotype_records__other_unphased_records_follow_the_flag) {
    Gene gene = { .gff_id = 1, .coords = { 100, 200 } };
    Gene *genes[] = { &gene };
    Records phase; kv_init(phase);
    Record *phased = record_in_genes(110, genes, 1, false);
    Record *unphased = record_in_genes(120, genes, 1, true);
    kv_push(Record *, phase, phased); kv_push(Record *, phase, unphased);

    Records hap_in, hap_ex;

    find_haplotype_records(&phase, phased, &gene, true, &hap_in, &hap_ex);
    ck_assert_int_eq(kv_size(hap_in), 1);
    ck_assert_int_eq(kv_size(hap_ex), 0);
    kv_destroy(hap_in); kv_destroy(hap_ex);

    find_haplotype_records(&phase, unphased, &gene, true, &hap_in, &hap_ex); // For the unphased record the phased one and itself
    ck_assert_int_eq(kv_size(hap_in), 2);
    ck_assert_int_eq(kv_size(hap_ex), 1);
    kv_destroy(hap_in); kv_destroy(hap_ex);

    find_haplotype_records(&phase, phased, &gene, false, &hap_in, &hap_ex);
    ck_assert_int_eq(kv_size(hap_in), 2);
    ck_assert_ptr_eq(kv_A(hap_in, 1), unphased);
    ck_assert_int_eq(kv_size(hap_ex), 1);
    kv_destroy(hap_in); kv_destroy(hap_ex);

    kv_destroy(phase);
    Record_destroy(phased); Record_destroy(unphased);
} END_TEST

Suite *haplotype_suite() {
    Suite *s = suite_create("Haplotype");

    TCase *tc_ghs = mocked_tcase("get_haplotype_sequence");

    suite_add_tcase(s, tc_ghs);
    tcase_add_checked_fixture(tc_ghs, checked_setup, checked_teardown);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__snvs_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__snv_outside_bounds);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_overlaps_bounds_start);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_overlaps_bound_end);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__insertion_next_to_bound_end);

    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_within_bounds);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_overlaps_bound_start);
    tcase_add_test(tc_ghs, check_get_haplotype_sequence__deletion_overlaps_bound_end);

    TCase *tc_view = mocked_tcase("get_haplotype_sequence_anchor_free");
    suite_add_tcase(s, tc_view);

    tcase_add_test(tc_view, check_get_haplotype_sequence__single_variants_give_the_old_sequence);
    tcase_add_test(tc_view, check_get_haplotype_sequence__a_variant_changing_only_bases_after_the_gene_is_ignored);
    tcase_add_test(tc_view, check_get_haplotype_sequence__separate_variants_give_the_old_sequence);
    tcase_add_test(tc_view, check_get_haplotype_sequence__deletion_anchored_on_a_deleted_base);
    tcase_add_test(tc_view, check_get_haplotype_sequence__insertion_anchored_on_a_deleted_base);

    TCase *tc_shs = mocked_tcase("splice_haplotype_sequence");
    suite_add_tcase(s, tc_shs);

    tcase_add_test(tc_shs, check_splice_haplotype_sequence__forward_removes_used_introns_and_utrs);
    tcase_add_test(tc_shs, check_splice_haplotype_sequence__retained_intron_stays_in_the_sequence_once);
    tcase_add_test(tc_shs, check_splice_haplotype_sequence__reverse_gives_the_same_ascending_sequence);
    tcase_add_test(tc_shs, check_splice_haplotype_sequence__one_retained_and_one_spliced_intron);
    tcase_add_test(tc_shs, check_splice_haplotype_sequence__empty_utr_removes_nothing);
    tcase_add_test(tc_shs, check_splice_haplotype_sequence__unknown_strand_appends_nothing);

    TCase *tc_coords = mocked_tcase("isoform_coordinates");
    suite_add_tcase(s, tc_coords);

    tcase_add_test(tc_coords, check_isoform_segments__spliced_introns_and_utrs_leave_the_exons);
    tcase_add_test(tc_coords, check_isoform_segments__a_retained_intron_is_part_of_one_segment);
    tcase_add_test(tc_coords, check_isoform_segments__one_retained_and_one_spliced_intron);
    tcase_add_test(tc_coords, check_isoform_segments__reverse_gives_the_same_ascending_segments);
    tcase_add_test(tc_coords, check_isoform_segments__an_empty_utr_at_the_gene_start_leaves_nothing_out);
    tcase_add_test(tc_coords, check_isoform_segments__a_region_that_starts_before_the_end_of_the_previous_one_is_skipped);
    tcase_add_test(tc_coords, check_isoform_segments__an_unknown_strand_has_none);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__forward);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__reverse_counts_from_the_other_end);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__a_retained_intron_is_counted);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_spliced_forward_gene);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_forward_gene_with_a_retained_intron);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_spliced_reverse_gene);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_on_a_reverse_gene_with_a_retained_intron);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_with_one_retained_and_one_spliced_intron);
    tcase_add_test(tc_coords, check_hap_range_to_isoform__every_range_agrees_with_the_probe_with_two_spliced_introns_and_an_empty_utr);
    tcase_add_test(tc_coords, check_position_in_ranges__is_half_open);

    TCase *tc_nmd_iso = mocked_tcase("nmd_escape_ranges_in_isoform");
    suite_add_tcase(s, tc_nmd_iso);

    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__forward_without_variants);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__every_kind_keeps_its_own_list);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__reverse_without_variants_is_the_same_in_transcript_order);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_deletion_inside_a_region_shortens_it);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_deletion_before_every_region_moves_them_but_keeps_the_exon_ends);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__an_insertion_inside_a_region_lengthens_it);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__an_insertion_at_the_end_of_a_region_is_inside_and_at_the_start_outside);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__reverse_an_insertion_at_the_end_of_a_region_is_inside_and_at_the_start_outside);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_retained_last_intron_is_distal);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_retained_intron_elsewhere_is_not_an_escape_region);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__reverse_with_a_retained_last_intron_is_distal);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_region_that_keeps_one_base_is_one_base_long);
    tcase_add_test(tc_nmd_iso, check_nmd_escape_ranges_in_isoform__a_region_without_bases_in_the_isoform_is_left_out);

    TCase *tc_ptv = mocked_tcase("ptv_stop_location");
    suite_add_tcase(s, tc_ptv);

    tcase_add_test(tc_ptv, check_ptv_stop__a_substitution_that_makes_a_stop_gives_the_index_of_its_codon);
    tcase_add_test(tc_ptv, check_ptv_stop__a_substitution_without_a_stop_has_none);
    tcase_add_test(tc_ptv, check_ptv_stop__the_original_stop_codon_never_counts);
    tcase_add_test(tc_ptv, check_ptv_stop__a_substitution_that_runs_into_the_next_codon_is_not_cut_short);
    tcase_add_test(tc_ptv, check_ptv_stop__an_insertion_in_frame_that_is_a_stop_gives_its_index);
    tcase_add_test(tc_ptv, check_ptv_stop__a_frameshifting_insertion_makes_a_stop_with_the_bases_after_it);
    tcase_add_test(tc_ptv, check_ptv_stop__a_deletion_makes_a_stop_when_the_bases_after_it_do);
    tcase_add_test(tc_ptv, check_ptv_stop__a_frameshift_finds_a_stop_further_on_that_the_codons_of_the_change_do_not_have);
    tcase_add_test(tc_ptv, check_ptv_stop__a_frameshift_that_reads_on_without_a_stop_has_none);
    tcase_add_test(tc_ptv, check_ptv_stop__a_deletion_in_frame_has_no_new_stop);
    tcase_add_test(tc_ptv, check_ptv_stop__a_position_outside_the_sequence_is_moved_into_it);
    tcase_add_test(tc_ptv, check_ptv_stop__an_empty_sequence_has_no_stop);
    tcase_add_test(tc_ptv, check_nmd_flag_for_stop__names_the_region_of_the_stop_and_nothing_outside_them_or_without_a_stop);
    tcase_add_test(tc_ptv, check_nmd_flag_for_stop__overlapping_regions_give_one_flag_in_a_fixed_order);
    tcase_add_test(tc_ptv, check_nmd_flag_for_stop__the_flags_are_distinct_bits);
    tcase_add_test(tc_ptv, check_ptv_stop__agrees_with_a_translation_of_every_small_variant);
    tcase_add_test(tc_ptv, check_nmd_location__forward_gene_from_the_change_to_the_flag);
    tcase_add_test(tc_ptv, check_nmd_location__reverse_gene_from_the_change_to_the_flag);
    tcase_add_test(tc_ptv, check_nmd_location__a_stop_is_placed_by_its_first_base);
    tcase_add_test(tc_ptv, check_nmd_location__an_inserted_stop_belongs_to_the_base_it_is_attached_to);
    tcase_add_test(tc_ptv, check_nmd_location__no_new_stop_has_no_flag);

    TCase *tc_old = mocked_tcase("old_splice_sequence");
    suite_add_tcase(s, tc_old);

    tcase_add_test(tc_old, check_old_forward_splice_sequence__duplicated_the_donor_base_of_a_retained_intron);

    TCase *tc_erc = mocked_tcase("enumerate_region_combinations");
    suite_add_tcase(s, tc_erc);

    tcase_add_test(tc_erc, check_enumerate_region_combinations__no_groups_gives_no_combinations);
    tcase_add_test(tc_erc, check_enumerate_region_combinations__only_empty_groups_gives_no_combinations);
    tcase_add_test(tc_erc, check_enumerate_region_combinations__one_group_gives_one_combination_per_region);
    tcase_add_test(tc_erc, check_enumerate_region_combinations__first_group_varies_fastest);
    tcase_add_test(tc_erc, check_enumerate_region_combinations__empty_groups_are_skipped);
    tcase_add_test(tc_erc, check_enumerate_region_combinations__regions_are_copied_by_value);

    TCase *tc_iso = mocked_tcase("isoform_helpers");
    suite_add_tcase(s, tc_iso);

    tcase_add_test(tc_iso, check_codon_start__rounds_down_to_a_multiple_of_three);
    tcase_add_test(tc_iso, check_codon_window__runs_from_the_codon_start_to_offset_plus_length);
    tcase_add_test(tc_iso, check_codon_window_within__keeps_the_window_inside_the_sequence);
    tcase_add_test(tc_iso, check_has_stop_codon__finds_each_stop_codon);
    tcase_add_test(tc_iso, check_has_stop_codon__reads_codons_from_start_in_steps_of_three);
    tcase_add_test(tc_iso, check_has_stop_codon__ignores_stops_out_of_frame);
    tcase_add_test(tc_iso, check_has_stop_codon__reads_the_whole_last_codon);
    tcase_add_test(tc_iso, check_has_stop_codon__zero_length_reads_nothing);
    tcase_add_test(tc_iso, check_has_stop_codon__does_not_read_past_the_sequence);
    tcase_add_test(tc_iso, check_has_stop_codon__a_codon_over_the_end_is_not_read);
    tcase_add_test(tc_iso, check_has_stop_codon__a_start_outside_the_sequence_reads_nothing);
    tcase_add_test(tc_iso, check_first_stop_codon__is_the_index_of_the_first_base_of_the_stop);
    tcase_add_test(tc_iso, check_first_stop_codon__is_minus_one_without_a_stop);
    tcase_add_test(tc_iso, check_first_stop_codon__gives_the_first_of_several);
    tcase_add_test(tc_iso, check_first_stop_codon__reads_from_start_in_steps_of_three_up_to_the_length);
    tcase_add_test(tc_iso, check_first_stop_codon__does_not_read_a_codon_that_is_not_entirely_in_the_sequence);
    tcase_add_test(tc_iso, check_first_stop_codon__can_leave_out_the_stop_codon_that_ends_the_sequence);
    tcase_add_test(tc_iso, check_first_stop_codon__a_frameshifted_sequence_that_reads_on_to_the_original_stop_has_no_new_stop);
    tcase_add_test(tc_iso, check_has_stop_codon__agrees_with_first_stop_codon);
    tcase_add_test(tc_iso, check_reverse_complement__complements_and_reverses);
    tcase_add_test(tc_iso, check_reverse_complement__leaves_other_characters_unchanged);
    tcase_add_test(tc_iso, check_reverse_complement__empty_sequence);
    tcase_add_test(tc_iso, check_isoform_prior__multiplies_the_priors_of_used_regions);
    tcase_add_test(tc_iso, check_isoform_prior__no_regions_gives_one);

    TCase *tc_hap = mocked_tcase("find_haplotype_records");
    suite_add_tcase(s, tc_hap);

    tcase_add_test(tc_hap, check_record_overlaps_gene__follows_the_changed_bases);
    tcase_add_test(tc_hap, check_find_haplotype_records__the_records_in_the_gene_in_order);
    tcase_add_test(tc_hap, check_find_haplotype_records__other_genes_other_places_and_unannotated_records_are_left_out);
    tcase_add_test(tc_hap, check_find_haplotype_records__other_unphased_records_follow_the_flag);

    return s;
}

