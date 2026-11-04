#include <check.h>
#include <cpliceai/range.h>
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>

#include "klib/kvec.h"
#include "mocks.h"
#include "record.h"
#include "triggers.h"
#include "utils.h"

#include "check_record.h"

// Builds a transcript whose exons have the given ref_bounds, copied as is.
static Transcript *make_transcript(const Range *ref_bounds, size_t n_exons) {
    Transcript *t = calloc(1, sizeof(Transcript));
    kv_init(t->exons);
    kv_init(t->introns);

    for (size_t i = 0; i < n_exons; i++) {
        Exon *e = calloc(1, sizeof(Exon));
        e->ref_bounds = malloc(sizeof(Range));
        *e->ref_bounds = ref_bounds[i];
        kv_push(Exon *, t->exons, e);
    }

    return t;
}

static void free_transcript(Transcript *t) {
    for (size_t i = 0; i < kv_size(t->exons); i++) {
        free(kv_A(t->exons, i)->ref_bounds);
        free(kv_A(t->exons, i));
    }

    kv_destroy(t->exons);
    kv_destroy(t->introns);
    free(t);
}

// Introns_from_exons: the exon ranges are half-open and ascending.

START_TEST(check_introns_from_exons__no_exons) {
    Introns introns = Introns_from_exons(NULL, 0);

    ck_assert_int_eq(kv_size(introns), 0);

    kv_destroy(introns);
}
END_TEST

START_TEST(check_introns_from_exons__single_exon_has_no_introns) {
    Range exons[] = { {100, 200} };
    Introns introns = Introns_from_exons(exons, 1);

    ck_assert_int_eq(kv_size(introns), 0);

    kv_destroy(introns);
}
END_TEST

START_TEST(check_introns_from_exons__two_exons) {
    // The intron covers bases 200..299, so 100 bases
    Range exons[] = { {100, 200}, {300, 400} };
    Introns introns = Introns_from_exons(exons, 2);

    ck_assert_int_eq(kv_size(introns), 1);
    ck_assert_int_eq(kv_A(introns, 0).start, 200);
    ck_assert_int_eq(kv_A(introns, 0).end, 300);
    ck_assert_int_eq(kv_A(introns, 0).len, 100);
    ck_assert_int_eq(kv_A(introns, 0).coding_offset, 100);

    kv_destroy(introns);
}
END_TEST

START_TEST(check_introns_from_exons__three_exons_accumulate_offset) {
    Range exons[] = { {100, 200}, {300, 350}, {500, 600} };
    Introns introns = Introns_from_exons(exons, 3);

    ck_assert_int_eq(kv_size(introns), 2);

    ck_assert_int_eq(kv_A(introns, 0).start, 200);
    ck_assert_int_eq(kv_A(introns, 0).end, 300);
    ck_assert_int_eq(kv_A(introns, 0).len, 100);
    ck_assert_int_eq(kv_A(introns, 0).coding_offset, 100);

    ck_assert_int_eq(kv_A(introns, 1).start, 350);
    ck_assert_int_eq(kv_A(introns, 1).end, 500);
    ck_assert_int_eq(kv_A(introns, 1).len, 150);
    ck_assert_int_eq(kv_A(introns, 1).coding_offset, 250);

    kv_destroy(introns);
}
END_TEST

START_TEST(check_introns_from_exons__one_intron_per_exon_gap) {
    Range exons[] = { {0, 10}, {20, 30}, {40, 50}, {60, 70}, {80, 90} };
    Introns introns = Introns_from_exons(exons, 5);

    ck_assert_int_eq(kv_size(introns), 4);
    for (int i = 0; i < 4; i++) {
        // Intron i lies between exon i and exon i+1
        ck_assert_int_eq(kv_A(introns, i).start, exons[i].end);
        ck_assert_int_eq(kv_A(introns, i).end, exons[i+1].start);
        ck_assert_int_eq(kv_A(introns, i).len, 10);
        ck_assert_int_eq(kv_A(introns, i).coding_offset, 10 * (i + 1));
    }

    kv_destroy(introns);
}
END_TEST

START_TEST(check_introns_from_exons__touching_exons_give_empty_intron) {
    // NOTE: not expected from a GFF, but the one-intron-per-gap rule keeps indices aligned with the exons
    Range exons[] = { {100, 200}, {200, 300}, {400, 500} };
    Introns introns = Introns_from_exons(exons, 3);

    ck_assert_int_eq(kv_size(introns), 2);
    ck_assert_int_eq(kv_A(introns, 0).len, 0);
    ck_assert_int_eq(kv_A(introns, 0).coding_offset, 0);
    ck_assert_int_eq(kv_A(introns, 1).len, 100);
    ck_assert_int_eq(kv_A(introns, 1).coding_offset, 100);

    kv_destroy(introns);
}
END_TEST

// Splice site bases and introns. Exons are half-open, so an exon [100, 200) has its first base at 100 and its last at 199.

START_TEST(check_exon_bases__forward_donor_is_last_base_and_acceptor_is_first) {
    Range bounds[] = { {100, 200} };
    Transcript *t = make_transcript(bounds, 1);

    ck_assert_int_eq(Exon_donor_base(kv_A(t->exons, 0), STRAND_FWD), 199);
    ck_assert_int_eq(Exon_acceptor_base(kv_A(t->exons, 0), STRAND_FWD), 100);

    free_transcript(t);
}
END_TEST

START_TEST(check_exon_bases__reverse_donor_is_first_base_and_acceptor_is_last) {
    Range bounds[] = { {100, 200} };
    Transcript *t = make_transcript(bounds, 1);

    ck_assert_int_eq(Exon_donor_base(kv_A(t->exons, 0), STRAND_REV), 100);
    ck_assert_int_eq(Exon_acceptor_base(kv_A(t->exons, 0), STRAND_REV), 199);

    free_transcript(t);
}
END_TEST

START_TEST(check_exon_bases__unknown_strand_behaves_as_forward) {
    Range bounds[] = { {100, 200} };
    Transcript *t = make_transcript(bounds, 1);

    ck_assert_int_eq(Exon_donor_base(kv_A(t->exons, 0), STRAND_UNK), 199);
    ck_assert_int_eq(Exon_acceptor_base(kv_A(t->exons, 0), STRAND_UNK), 100);

    free_transcript(t);
}
END_TEST

START_TEST(check_intron_between_sites__forward) {
    // Exons [100, 200) and [300, 400): the donor is base 199 and the acceptor is base 300
    Range intron = Intron_between_sites(STRAND_FWD, 199, 300);

    ck_assert_int_eq(intron.start, 200);
    ck_assert_int_eq(intron.end, 300);
}
END_TEST

START_TEST(check_intron_between_sites__reverse_gives_the_same_ascending_range) {
    // The same exons on the reverse strand: the donor is the first base of the upper exon, the acceptor the last base of the lower one
    Range intron = Intron_between_sites(STRAND_REV, 300, 199);

    ck_assert_int_eq(intron.start, 200);
    ck_assert_int_eq(intron.end, 300);
}
END_TEST

START_TEST(check_intron_bases__are_the_inverse_of_intron_between_sites) {
    uint32_t strands[] = { STRAND_FWD, STRAND_REV };
    for (int i = 0; i < 2; i++) {
        const int64_t donor = strands[i] == STRAND_FWD ? 199 : 300;
        const int64_t acceptor = strands[i] == STRAND_FWD ? 300 : 199;

        Range intron = Intron_between_sites(strands[i], donor, acceptor);

        ck_assert_int_eq(Intron_donor_base(intron, strands[i]), donor);
        ck_assert_int_eq(Intron_acceptor_base(intron, strands[i]), acceptor);
    }
}
END_TEST

START_TEST(check_utr_range__forward) {
    Range gene = { 100, 1000 };
    Range cds = { 300, 400 };

    Range utr5 = Utr_range(UTR5, STRAND_FWD, gene, cds);
    ck_assert_int_eq(utr5.start, 100);
    ck_assert_int_eq(utr5.end, 300);

    Range utr3 = Utr_range(UTR3, STRAND_FWD, gene, cds);
    ck_assert_int_eq(utr3.start, 400);
    ck_assert_int_eq(utr3.end, 1000);
}
END_TEST

START_TEST(check_utr_range__reverse_is_mirrored_and_still_ascending) {
    Range gene = { 100, 1000 };
    Range cds = { 300, 400 };

    Range utr5 = Utr_range(UTR5, STRAND_REV, gene, cds);
    ck_assert_int_eq(utr5.start, 400);
    ck_assert_int_eq(utr5.end, 1000);

    Range utr3 = Utr_range(UTR3, STRAND_REV, gene, cds);
    ck_assert_int_eq(utr3.start, 100);
    ck_assert_int_eq(utr3.end, 300);
}
END_TEST

START_TEST(check_utr_range__empty_when_cds_touches_the_gene_edge) {
    Range gene = { 100, 1000 };

    Range utr5 = Utr_range(UTR5, STRAND_FWD, gene, (Range) { 100, 400 });
    ck_assert_int_eq(range_len(utr5), 0);

    Range utr3 = Utr_range(UTR3, STRAND_FWD, gene, (Range) { 300, 1000 });
    ck_assert_int_eq(range_len(utr3), 0);
}
END_TEST

// Exon_splice_window: three exons [100, 200), [300, 400), [500, 600). The windows come from the reference bounds of the exons:
// the sites they are compared with carry reference coordinates, so the variants of a haplotype do not move them.

static Transcript *make_windowed_transcript(uint32_t strand) {
    Range bounds[] = { {100, 200}, {300, 400}, {500, 600} };
    Transcript *t = make_transcript(bounds, 3);
    t->strand = strand;

    return t;
}

START_TEST(check_exon_splice_window__forward_donor_looks_to_the_exon_above) {
    Transcript *t = make_windowed_transcript(STRAND_FWD);
    Range w;

    ck_assert(Exon_splice_window(t, 0, true, &w));
    ck_assert_int_eq(w.start, 101);
    ck_assert_int_eq(w.end, 300);

    ck_assert(Exon_splice_window(t, 1, true, &w));
    ck_assert_int_eq(w.start, 301);
    ck_assert_int_eq(w.end, 500);

    ck_assert(!Exon_splice_window(t, 2, true, &w)); // The last exon has no downstream acceptor to pair with

    free_transcript(t);
}
END_TEST

START_TEST(check_exon_splice_window__forward_acceptor_looks_to_the_exon_below) {
    Transcript *t = make_windowed_transcript(STRAND_FWD);
    Range w;

    ck_assert(!Exon_splice_window(t, 0, false, &w)); // The first exon has no acceptor

    ck_assert(Exon_splice_window(t, 1, false, &w));
    ck_assert_int_eq(w.start, 199);
    ck_assert_int_eq(w.end, 398);

    ck_assert(Exon_splice_window(t, 2, false, &w));
    ck_assert_int_eq(w.start, 399);
    ck_assert_int_eq(w.end, 598);

    free_transcript(t);
}
END_TEST

START_TEST(check_exon_splice_window__reverse_sites_use_the_window_of_the_same_border) {
    Transcript *t = make_windowed_transcript(STRAND_REV);
    Range w;

    // A reverse acceptor sits on the upper border, like a forward donor
    ck_assert(Exon_splice_window(t, 0, false, &w));
    ck_assert_int_eq(w.start, 101);
    ck_assert_int_eq(w.end, 300);
    ck_assert(!Exon_splice_window(t, 2, false, &w));

    // A reverse donor sits on the lower border, like a forward acceptor
    ck_assert(!Exon_splice_window(t, 0, true, &w));
    ck_assert(Exon_splice_window(t, 1, true, &w));
    ck_assert_int_eq(w.start, 199);
    ck_assert_int_eq(w.end, 398);

    free_transcript(t);
}
END_TEST

START_TEST(check_exon_splice_window__windows_are_ascending_on_both_strands) {
    uint32_t strands[] = { STRAND_FWD, STRAND_REV };
    for (int s = 0; s < 2; s++) {
        Transcript *t = make_windowed_transcript(strands[s]);

        for (size_t ei = 0; ei < 3; ei++) {
            Range w;
            if (Exon_splice_window(t, ei, true, &w)) ck_assert(w.start <= w.end);
            if (Exon_splice_window(t, ei, false, &w)) ck_assert(w.start <= w.end);
        }

        free_transcript(t);
    }
}
END_TEST

START_TEST(check_exon_splice_window__unknown_strand_has_no_windows) {
    Transcript *t = make_windowed_transcript(STRAND_UNK);
    Range w;

    ck_assert(!Exon_splice_window(t, 1, true, &w));
    ck_assert(!Exon_splice_window(t, 1, false, &w));

    free_transcript(t);
}
END_TEST

// Transcript_build: the input is inclusive GFF data, the output is half-open.

START_TEST(check_transcript_build__converts_coordinates_to_half_open) {
    Range cds[] = { {120, 199} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 199} };
    TranscriptData data = { .id = 7, .start = 100, .end = 199, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 1 };

    Transcript *t = Transcript_build(&data);

    ck_assert_int_eq(t->id, 7);
    ck_assert_int_eq(t->strand, STRAND_FWD);
    ck_assert_int_eq(t->start, 100);
    ck_assert_int_eq(t->end, 200);
    ck_assert_int_eq(t->ncds, 1);
    ck_assert_int_eq(t->cds[0].start, 120);
    ck_assert_int_eq(t->cds[0].end, 200);
    ck_assert_int_eq(kv_size(t->exons), 1);
    ck_assert_int_eq(kv_A(t->exons, 0)->ref_bounds->start, 100);
    ck_assert_int_eq(kv_A(t->exons, 0)->ref_bounds->end, 200);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__new_transcript_has_no_phase_state) {
    Range cds[] = { {100, 199} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 199} };
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = STRAND_REV,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 1 };

    Transcript *t = Transcript_build(&data);

    ck_assert_ptr_null(kv_A(t->exons, 0)->splice_sites);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__does_not_alias_input) {
    Range cds[] = { {100, 199} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 199} };
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 1 };

    Transcript *t = Transcript_build(&data);

    // The caller's data stays inclusive and untouched, and the transcript owns its own copies
    ck_assert_int_eq(exons[0].end, 199);
    ck_assert_int_eq(cds[0].end, 199);
    ck_assert_ptr_ne(kv_A(t->exons, 0)->ref_bounds, &exons[0]);
    ck_assert_ptr_ne(t->cds, cds);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__exon_without_cds_does_not_use_up_a_coding_region) {
    // The first exon is entirely UTR, so the two coding regions belong to the second and third exons
    Range cds[] = { {300, 349}, {500, 599} };
    int constitutive[] = { 1, 0 };
    Range exons[] = { {100, 199}, {300, 349}, {500, 599} };
    TranscriptData data = { .id = 1, .start = 100, .end = 599, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 3 };

    Transcript *t = Transcript_build(&data);

    ck_assert_int_eq(kv_size(t->exons), 3);

    Exon *utr = kv_A(t->exons, 0);
    ck_assert_ptr_null(utr->cds);
    ck_assert(!utr->constitutive);
    ck_assert(!utr->in_frame);

    ck_assert_ptr_eq(kv_A(t->exons, 1)->cds, &t->cds[0]);
    ck_assert(kv_A(t->exons, 1)->constitutive);

    ck_assert_ptr_eq(kv_A(t->exons, 2)->cds, &t->cds[1]);
    ck_assert(!kv_A(t->exons, 2)->constitutive);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__in_frame_depends_on_cds_length) {
    // Inclusive lengths: 102 bases is divisible by three, 100 is not
    Range cds[] = { {100, 201}, {300, 399} };
    int constitutive[] = { 1, 1 };
    Range exons[] = { {100, 201}, {300, 399} };
    TranscriptData data = { .id = 1, .start = 100, .end = 399, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 2 };

    Transcript *t = Transcript_build(&data);

    ck_assert(kv_A(t->exons, 0)->in_frame);
    ck_assert(!kv_A(t->exons, 1)->in_frame);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__derives_introns_from_the_half_open_exons) {
    Range cds[] = { {100, 199}, {300, 349}, {500, 599} };
    int constitutive[] = { 1, 1, 1 };
    Range exons[] = { {100, 199}, {300, 349}, {500, 599} };
    TranscriptData data = { .id = 1, .start = 100, .end = 599, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };

    Transcript *t = Transcript_build(&data);

    // Introns cover 200..299 and 350..499
    ck_assert_int_eq(kv_size(t->introns), 2);
    ck_assert_int_eq(kv_A(t->introns, 0).start, 200);
    ck_assert_int_eq(kv_A(t->introns, 0).end, 300);
    ck_assert_int_eq(kv_A(t->introns, 0).len, 100);
    ck_assert_int_eq(kv_A(t->introns, 0).coding_offset, 100);
    ck_assert_int_eq(kv_A(t->introns, 1).start, 350);
    ck_assert_int_eq(kv_A(t->introns, 1).end, 500);
    ck_assert_int_eq(kv_A(t->introns, 1).len, 150);
    ck_assert_int_eq(kv_A(t->introns, 1).coding_offset, 250);

    Transcript_destroy(t);
}
END_TEST

START_TEST(check_transcript_build__no_exons) {
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = STRAND_FWD,
                            .cds = NULL, .constitutive = NULL, .ncds = 0, .exons = NULL, .nexons = 0 };

    Transcript *t = Transcript_build(&data);

    ck_assert_int_eq(kv_size(t->exons), 0);
    ck_assert_int_eq(kv_size(t->introns), 0);

    Transcript_destroy(t);
}
END_TEST

/*
* Transcript_nmd_regions: the NMD escape regions of the reference transcript (Regional nonsense constraint, Nature Communications
* 2026): start-proximal (the first 150 coding bases), long exon (more than 400 bases before the last base of the exon, which
* leaves out the last exon) and distal (the last exon, the last 50 bases of the penultimate exon and the intron between them).
* Positions are those of the first base of a stop codon. Every range is half-open and ascending, whatever the strand.
*/
static void check_nmd_range(const Ranges *ranges, const size_t i, const int64_t start, const int64_t end) {
    ck_assert_uint_gt(kv_size(*ranges), i);
    ck_assert_int_eq(kv_A(*ranges, i).start, start);
    ck_assert_int_eq(kv_A(*ranges, i).end, end);
}

// Forward: exons [1000, 1200), [2000, 2600), [3000, 3300) with the coding regions [1100, 1200), [2000, 2600), [3000, 3150)
static Transcript *nmd_forward_transcript(void) {
    Range cds[] = { {1100, 1199}, {2000, 2599}, {3000, 3149} };
    int constitutive[] = { 1, 1, 1 };
    Range exons[] = { {1000, 1199}, {2000, 2599}, {3000, 3299} };
    TranscriptData data = { .id = 1, .start = 1000, .end = 3299, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };

    return Transcript_build(&data);
}

// Reverse: exons [1000, 1300), [2000, 2600), [3000, 3200) with the coding regions [1150, 1300), [2000, 2600), [3000, 3100).
// The transcript starts at the high end, so [1000, 1300) is the last exon.
static Transcript *nmd_reverse_transcript(void) {
    Range cds[] = { {1150, 1299}, {2000, 2599}, {3000, 3099} };
    int constitutive[] = { 1, 1, 1 };
    Range exons[] = { {1000, 1299}, {2000, 2599}, {3000, 3199} };
    TranscriptData data = { .id = 2, .start = 1000, .end = 3199, .strand = STRAND_REV,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };

    return Transcript_build(&data);
}

START_TEST(check_nmd_regions__forward_start_proximal_is_the_first_150_coding_bases) {
    Transcript *t = nmd_forward_transcript();
    NmdRegions regions = Transcript_nmd_regions(t);

    // 100 bases in the first coding region, the other 50 at the start of the second
    ck_assert_int_eq(kv_size(regions.start_proximal), 2);
    check_nmd_range(&regions.start_proximal, 0, 1100, 1200);
    check_nmd_range(&regions.start_proximal, 1, 2000, 2050);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__forward_long_exon_leaves_the_last_401_bases_of_an_exon_out) {
    Transcript *t = nmd_forward_transcript();
    NmdRegions regions = Transcript_nmd_regions(t);

    // The middle exon has 600 bases: the first 199 are more than 400 before its last base. The first exon is too short, and the
    // last exon has no donor
    ck_assert_int_eq(kv_size(regions.long_exon), 1);
    check_nmd_range(&regions.long_exon, 0, 2000, 2199);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__forward_distal_is_the_last_exon_the_intron_and_the_end_of_the_penultimate_exon) {
    Transcript *t = nmd_forward_transcript();
    NmdRegions regions = Transcript_nmd_regions(t);

    ck_assert_int_eq(kv_size(regions.distal), 1);
    check_nmd_range(&regions.distal, 0, 2550, 3300);   // The last 50 bases of [2000, 2600), the intron, and [3000, 3300)

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__reverse_regions_are_mirrored_and_still_ascending) {
    Transcript *t = nmd_reverse_transcript();
    NmdRegions regions = Transcript_nmd_regions(t);

    // The first 150 coding bases are 100 at [3000, 3100) and the 50 at the high end of [2000, 2600). Listed ascending
    ck_assert_int_eq(kv_size(regions.start_proximal), 2);
    check_nmd_range(&regions.start_proximal, 0, 2550, 2600);
    check_nmd_range(&regions.start_proximal, 1, 3000, 3100);

    // The middle exon is the penultimate one, not the last: the first 199 bases in transcript order are its top 199
    ck_assert_int_eq(kv_size(regions.long_exon), 1);
    check_nmd_range(&regions.long_exon, 0, 2401, 2600);

    // The last exon is [1000, 1300); the last 50 bases of the penultimate exon in transcript order are its lowest 50
    ck_assert_int_eq(kv_size(regions.distal), 1);
    check_nmd_range(&regions.distal, 0, 1000, 2050);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__an_exon_of_401_bases_has_no_long_exon_region_and_one_of_402_has_one_base) {
    Range cds[] = { {100, 1199} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 500}, {1000, 1401}, {2000, 2099} };   // 401, 402 and 100 bases
    TranscriptData data = { .id = 1, .start = 100, .end = 2099, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 3 };
    Transcript *t = Transcript_build(&data);
    NmdRegions regions = Transcript_nmd_regions(t);

    ck_assert_int_eq(kv_size(regions.long_exon), 1);
    check_nmd_range(&regions.long_exon, 0, 1000, 1001);   // Only the first base is 401 before the last base

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__a_single_exon_transcript_is_all_distal) {
    Range cds[] = { {550, 799} };
    int constitutive[] = { 1 };
    Range exons[] = { {500, 899} };
    TranscriptData data = { .id = 1, .start = 500, .end = 899, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 1 };
    Transcript *t = Transcript_build(&data);
    NmdRegions regions = Transcript_nmd_regions(t);

    ck_assert_int_eq(kv_size(regions.start_proximal), 1);
    check_nmd_range(&regions.start_proximal, 0, 550, 700);
    ck_assert_int_eq(kv_size(regions.long_exon), 0);   // The only exon is the last one
    ck_assert_int_eq(kv_size(regions.distal), 1);
    check_nmd_range(&regions.distal, 0, 500, 900);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__a_short_penultimate_exon_is_distal_as_a_whole_and_a_short_cds_is_start_proximal_as_a_whole) {
    Range cds[] = { {110, 129}, {200, 219} };   // 20 and 20 coding bases
    int constitutive[] = { 1, 1 };
    Range exons[] = { {100, 129}, {200, 299} };   // The penultimate exon has 30 bases
    TranscriptData data = { .id = 1, .start = 100, .end = 299, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 2 };
    Transcript *t = Transcript_build(&data);
    NmdRegions regions = Transcript_nmd_regions(t);

    ck_assert_int_eq(kv_size(regions.start_proximal), 2);
    check_nmd_range(&regions.start_proximal, 0, 110, 130);
    check_nmd_range(&regions.start_proximal, 1, 200, 220);
    ck_assert_int_eq(kv_size(regions.distal), 1);
    check_nmd_range(&regions.distal, 0, 100, 300);   // Not 50 bases before the intron: the exon has 30

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__a_transcript_without_exons_has_no_regions) {
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = STRAND_FWD,
                            .cds = NULL, .constitutive = NULL, .ncds = 0, .exons = NULL, .nexons = 0 };
    Transcript *t = Transcript_build(&data);
    NmdRegions regions = Transcript_nmd_regions(t);

    ck_assert_int_eq(kv_size(regions.start_proximal), 0);
    ck_assert_int_eq(kv_size(regions.long_exon), 0);
    ck_assert_int_eq(kv_size(regions.distal), 0);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

START_TEST(check_nmd_regions__an_unknown_strand_is_read_as_forward) {
    Range cds[] = { {1100, 1199}, {2000, 2599}, {3000, 3149} };
    int constitutive[] = { 1, 1, 1 };
    Range exons[] = { {1000, 1199}, {2000, 2599}, {3000, 3299} };
    TranscriptData data = { .id = 1, .start = 1000, .end = 3299, .strand = STRAND_UNK,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };
    Transcript *t = Transcript_build(&data);
    NmdRegions regions = Transcript_nmd_regions(t);

    check_nmd_range(&regions.distal, 0, 2550, 3300);

    NmdRegions_destroy(&regions);
    Transcript_destroy(t);
}
END_TEST

// Transcript_phase_copy: two coding exons [100, 200) and [300, 400), the first with a CDS.

static Transcript *make_shared_transcript(void) {
    Range cds[] = { {120, 199} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 199}, {300, 399} };
    TranscriptData data = { .id = 7, .start = 100, .end = 399, .strand = STRAND_REV,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 2 };

    return Transcript_build(&data);
}

// Builds splice sites the way main() does: arrays from init_splice_site and a Records vector owned by each site
static SpliceSites *make_splice_sites(void) {
    SpliceSites *sites = malloc(sizeof(SpliceSites));
    init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
    init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
    init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
    init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

    Records *donor_records = malloc(sizeof(Records));
    kv_init(*donor_records);
    sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { donor_records, 199, 199, 0, 1, 0, 1.0, 1.0, 1.0 };

    Records *cryptic_records = malloc(sizeof(Records));
    kv_init(*cryptic_records);
    sites->cryptic_donors[sites->n_cryptic_donors++] = (SpliceSite) { cryptic_records, 150, 199, 1, 0, 0, 0.5, 0.5, 0.5 };

    return sites;
}

START_TEST(check_transcript_phase_copy__has_the_values_of_the_shared_transcript) {
    Transcript *shared = make_shared_transcript();
    Transcript *copy = Transcript_phase_copy(shared);

    ck_assert_int_eq(copy->id, shared->id);
    ck_assert_int_eq(copy->start, shared->start);
    ck_assert_int_eq(copy->end, shared->end);
    ck_assert_int_eq(copy->strand, shared->strand);
    ck_assert_int_eq(copy->ncds, shared->ncds);
    ck_assert_int_eq(kv_size(copy->exons), 2);

    for (int i = 0; i < 2; i++) {
        Exon *c = kv_A(copy->exons, i), *s = kv_A(shared->exons, i);
        ck_assert_int_eq(c->ref_bounds->start, s->ref_bounds->start);
        ck_assert_int_eq(c->ref_bounds->end, s->ref_bounds->end);
        ck_assert(c->constitutive == s->constitutive);
        ck_assert(c->in_frame == s->in_frame);
    }

    Transcript_phase_destroy(copy);
    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_transcript_phase_copy__shares_the_read_only_data_but_not_the_exons) {
    Transcript *shared = make_shared_transcript();
    Transcript *copy = Transcript_phase_copy(shared);

    ck_assert_ptr_eq(copy->cds, shared->cds);
    ck_assert_ptr_eq(copy->introns.a, shared->introns.a);
    ck_assert_int_eq(kv_size(copy->introns), kv_size(shared->introns));

    for (int i = 0; i < 2; i++) {
        Exon *c = kv_A(copy->exons, i), *s = kv_A(shared->exons, i);
        ck_assert_ptr_ne(c, s);                         // The Exon structs are the copy's own
        ck_assert_ptr_eq(c->ref_bounds, s->ref_bounds); // The ranges are shared
        ck_assert_ptr_eq(c->cds, s->cds);
    }
    ck_assert_ptr_eq(kv_A(copy->exons, 0)->cds, &shared->cds[0]); // Exon.cds still points into Transcript.cds

    Transcript_phase_destroy(copy);
    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_transcript_phase_copy__starts_without_phase_state) {
    Transcript *shared = make_shared_transcript();

    // Even if the shared exon had phase state, the copy starts clean
    Exon *shared_exon = kv_A(shared->exons, 0);
    shared_exon->splice_sites = make_splice_sites();

    Transcript *copy = Transcript_phase_copy(shared);
    for (int i = 0; i < 2; i++) {
        ck_assert_ptr_null(kv_A(copy->exons, i)->splice_sites);
    }

    Transcript_phase_destroy(copy);
    SpliceSites_destroy(shared_exon->splice_sites);
    shared_exon->splice_sites = NULL;
    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_transcript_phase_copy__phase_state_is_independent) {
    Transcript *shared = make_shared_transcript();
    Transcript *first = Transcript_phase_copy(shared);
    Transcript *second = Transcript_phase_copy(shared);

    Exon *first_exon = kv_A(first->exons, 0);
    first_exon->splice_sites = make_splice_sites();

    ck_assert_ptr_null(kv_A(second->exons, 0)->splice_sites);
    ck_assert_ptr_null(kv_A(shared->exons, 0)->splice_sites);

    Transcript_phase_destroy(first);
    Transcript_phase_destroy(second);
    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_transcript_phase_copy__destroy_leaves_the_shared_transcript_intact) {
    Transcript *shared = make_shared_transcript();
    Transcript *copy = Transcript_phase_copy(shared);

    for (int i = 0; i < 2; i++) {
        Exon *e = kv_A(copy->exons, i);
        e->splice_sites = make_splice_sites();
    }

    Transcript_phase_destroy(copy);

    // Run this test under AddressSanitizer to also check that nothing leaked or was freed twice
    ck_assert_int_eq(kv_size(shared->exons), 2);
    ck_assert_int_eq(kv_A(shared->exons, 0)->ref_bounds->start, 100);
    ck_assert_int_eq(kv_A(shared->exons, 0)->ref_bounds->end, 200);
    ck_assert_int_eq(kv_A(shared->exons, 1)->ref_bounds->start, 300);
    ck_assert_int_eq(shared->cds[0].start, 120);
    ck_assert_int_eq(kv_size(shared->introns), 1);

    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_transcript_phase_copy__transcript_without_exons) {
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = STRAND_FWD,
                            .cds = NULL, .constitutive = NULL, .ncds = 0, .exons = NULL, .nexons = 0 };
    Transcript *shared = Transcript_build(&data);

    Transcript *copy = Transcript_phase_copy(shared);
    ck_assert_int_eq(kv_size(copy->exons), 0);

    Transcript_phase_destroy(copy);
    Transcript_destroy(shared);
}
END_TEST

START_TEST(check_splice_sites_destroy__empty_sites) {
    SpliceSites *sites = malloc(sizeof(SpliceSites));
    init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
    init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
    init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
    init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

    SpliceSites_destroy(sites);
}
END_TEST

// Record lifecycle. Frees cannot be asserted, so these tests are only meaningful under AddressSanitizer (or valgrind), where a
// leak, double free or free of an uninitialised pointer is reported.

// Isoforms as build_isoforms_for_phase leaves them: valid, n isoforms in the vector, and rescue probabilities as
// calculate_rescue_probabilities allocates them.
static Isoforms valid_isoforms(const uint32_t n) {
    Isoforms isoforms = {0};
    isoforms.valid = 1;
    isoforms.m_rescue_types = 1;

    kv_init(isoforms.isoforms);
    for (uint32_t i = 0; i < n; i++) {
        Isoform isoform = {0};
        kv_push(Isoform, isoforms.isoforms, isoform);
    }

    RescueProbabilities *probs = malloc(sizeof(RescueProbabilities));
    probs->n_isoforms = calloc(1, sizeof(uint32_t));
    probs->rescues = calloc(1, sizeof(uint8_t));
    probs->rescue_types = calloc(1, sizeof(enum RescueType));
    probs->rescue_probabilities = calloc(1, sizeof(float));
    probs->affected_bases = calloc(1, sizeof(uint32_t));
    isoforms.rescue_type_probabilities = probs;

    return isoforms;
}

// A record without VEP consequences; Record_init keeps its own copy of the variant
static Record *plain_record(void) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,T", 0, 1, true, -1);

    Record *r = Record_init(v, NULL);

    bcf_destroy(v);
    bcf_hdr_destroy(hdr);

    return r;
}

START_TEST(check_record_init__starts_without_flags_or_isoforms) {
    Record *r = plain_record();

    ck_assert_ptr_null(r->flags);
    ck_assert_ptr_null(r->isoforms);

    Record_destroy(r);
} END_TEST

START_TEST(check_isoforms_destroy__frees_isoforms_and_rescue_probabilities) {
    Isoforms isoforms = valid_isoforms(3);
    ck_assert_int_eq(kv_size(isoforms.isoforms), 3);

    Isoforms_destroy(&isoforms);
} END_TEST

START_TEST(check_isoforms_destroy__no_isoforms) {
    Isoforms isoforms = valid_isoforms(0);

    Isoforms_destroy(&isoforms);
} END_TEST

// One transcript, so NUM_PHASES isoform slots: the first is built, the others are not and must not be touched
START_TEST(check_record_destroy__frees_flags_and_built_isoform_slots) {
    Record *r = plain_record();

    r->tscripts = malloc(sizeof(Transcripts));
    kv_init(*r->tscripts);
    kv_push(Transcript *, *r->tscripts, NULL);

    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));
    r->isoforms = malloc(NUM_PHASES * sizeof(Isoforms));
    for (int k = 0; k < NUM_PHASES; k++) r->isoforms[k].valid = 0;
    r->isoforms[0] = valid_isoforms(2);

    Record_destroy(r);
} END_TEST

// One gene, so NUM_PHASES slots per score array: the first phase has scores, the others are NULL
START_TEST(check_record_destroy__frees_the_scores) {
    Record *r = plain_record();

    r->genes = malloc(sizeof(Genes));
    kv_init(*r->genes);
    kv_push(Gene *, *r->genes, NULL);

    r->native_scores = calloc(NUM_PHASES, sizeof(Score *));
    r->diff_scores = calloc(NUM_PHASES, sizeof(Score *));
    r->haplotype_scores = calloc(NUM_PHASES, sizeof(Score *));
    r->native_scores[0] = calloc(1, sizeof(Score));
    r->diff_scores[0] = calloc(1, sizeof(Score));
    r->haplotype_scores[0] = calloc(1, sizeof(Score));

    Record_destroy(r);
} END_TEST

/*
* parse_csq_hdr_fmt and parse_vep_transcript_consequences: the VEP annotation of a record, one entry per transcript feature.
*/

// The header of mock_hdr with an INFO/CSQ field of four VEP fields
static bcf_hdr_t *csq_hdr(void) {
    bcf_hdr_t *hdr = mock_hdr();

    ck_assert_int_eq(bcf_hdr_append(hdr, "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature\">"), 0);
    ck_assert_int_eq(bcf_hdr_sync(hdr), 0);

    return hdr;
}

// A record on chr1 with the given INFO/CSQ, or none when csq is NULL
static bcf1_t *csq_record(bcf_hdr_t *hdr, const char *csq) {
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, -1);
    if (csq != NULL) ck_assert_int_eq(bcf_update_info_string(hdr, v, "CSQ", csq), 0);

    return v;
}

START_TEST(check_parse_csq_hdr_fmt__the_position_of_every_field) {
    bcf_hdr_t *hdr = csq_hdr();

    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL);
    ck_assert_ptr_nonnull(fmt);
    ck_assert_str_eq(fmt->tag, "CSQ");
    ck_assert_int_eq(fmt->consequence_id, 1);
    ck_assert_int_eq(fmt->feature_type_id, 2);
    ck_assert_int_eq(fmt->feature_id, 3);

    int id;
    ck_assert_int_eq(khash_str2int_get(fmt->str2id, "Allele", &id), 0);
    ck_assert_int_eq(id, 0);
    ck_assert_int_eq(khash_str2int_get(fmt->str2id, "Consequence", &id), 0);
    ck_assert_int_eq(id, 1);
    ck_assert_int_eq(khash_str2int_get(fmt->str2id, "Feature_type", &id), 0);
    ck_assert_int_eq(id, 2);
    ck_assert_int_eq(khash_str2int_get(fmt->str2id, "Feature", &id), 0);
    ck_assert_int_eq(id, 3);
    ck_assert_str_eq(fmt->id2str[3], "Feature");

    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

START_TEST(check_parse_csq_hdr_fmt__a_header_without_csq_gives_null) {
    bcf_hdr_t *hdr = mock_hdr();

    ck_assert_ptr_null(parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL));

    bcf_hdr_destroy(hdr);
} END_TEST

// Only the transcript features are kept, with their consequences split on '&'
START_TEST(check_parse_vep_transcript_consequences__one_entry_per_transcript) {
    bcf_hdr_t *hdr = csq_hdr();
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL);
    bcf1_t *v = csq_record(hdr, "G|missense_variant&splice_region_variant|Transcript|ENST1,G|upstream_gene_variant|Transcript|ENST2,G|regulatory_region_variant|RegulatoryFeature|ENSR1");

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_ptr_nonnull(csqs);
    ck_assert_int_eq(kv_size(*csqs), 2); // The regulatory feature is not a transcript
    ck_assert_str_eq(kv_A(*csqs, 0)->tid, "ENST1");
    ck_assert_int_eq(kv_A(*csqs, 0)->num_consequences, 2);
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[0], "missense_variant");
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[1], "splice_region_variant");
    ck_assert_str_eq(kv_A(*csqs, 1)->tid, "ENST2");
    ck_assert_int_eq(kv_A(*csqs, 1)->num_consequences, 1);
    for (int i = 0; i < 3; i++) ck_assert_int_eq(kv_A(*csqs, 0)->triggers[i], 0); // No trigger yet
    for (int p = 0; p < NUM_PHASES; p++) ck_assert_int_eq(kv_A(*csqs, 0)->reading_frame[p], 0);

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

// A record that has no INFO/CSQ has no consequences: an empty list, as for a record that has no transcript feature. It used to be
// read through a NULL pointer.
START_TEST(check_parse_vep_transcript_consequences__a_record_without_csq_has_no_entries) {
    bcf_hdr_t *hdr = csq_hdr();
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL);
    bcf1_t *v = csq_record(hdr, NULL);

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_ptr_nonnull(csqs);
    ck_assert_int_eq(kv_size(*csqs), 0);

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

START_TEST(check_parse_vep_transcript_consequences__only_other_features_has_no_entries) {
    bcf_hdr_t *hdr = csq_hdr();
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL);
    bcf1_t *v = csq_record(hdr, "G|regulatory_region_variant|RegulatoryFeature|ENSR1");

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_ptr_nonnull(csqs);
    ck_assert_int_eq(kv_size(*csqs), 0);

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

// The empty list is a usable annotation: the record is made from it and is destroyed with it
START_TEST(check_record_init__an_empty_annotation_is_kept_and_freed) {
    bcf_hdr_t *hdr = csq_hdr();
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, NULL);
    bcf1_t *v = csq_record(hdr, NULL);

    Record *r = Record_init(v, parse_vep_transcript_consequences(v, fmt, hdr));

    ck_assert_ptr_nonnull(r->vep);
    ck_assert_int_eq(kv_size(*r->vep), 0);

    Record_destroy(r);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

/*
* The annotation of bcftools csq: another tag and columns, names of the consequences that differ from VEP's, entries that point to
* another record and entries without a transcript.
*/

// The header of mock_hdr with an INFO field `tag` in the format of bcftools csq
static bcf_hdr_t *bcsq_hdr(const char *tag) {
    bcf_hdr_t *hdr = mock_hdr();

    char line[512];
    snprintf(line, sizeof(line), "##INFO=<ID=%s,Number=.,Type=String,Description=\"Haplotype-aware consequence annotation from BCFtools/csq, see "
             "http://samtools.github.io/bcftools/howtos/csq-calling.html for details. "
             "Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change\">", tag);
    ck_assert_int_eq(bcf_hdr_append(hdr, line), 0);
    ck_assert_int_eq(bcf_hdr_sync(hdr), 0);

    return hdr;
}

// A record on chr1 with the value in INFO/<tag>
static bcf1_t *tagged_record(bcf_hdr_t *hdr, const char *tag, const char *value) {
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, -1);
    ck_assert_int_eq(bcf_update_info_string(hdr, v, tag, value), 0);

    return v;
}

START_TEST(check_parse_csq_hdr_fmt__bcftools_uses_bcsq_and_its_columns) {
    bcf_hdr_t *hdr = bcsq_hdr("BCSQ");

    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, NULL);

    ck_assert_ptr_nonnull(fmt);
    ck_assert_str_eq(fmt->tag, "BCSQ");
    ck_assert_int_eq(fmt->kind, ANNOTATION_BCFTOOLS);
    ck_assert_int_eq(fmt->consequence_id, 0);
    ck_assert_int_eq(fmt->feature_type_id, -1); // There is no such column
    ck_assert_int_eq(fmt->feature_id, 2); // The transcript

    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

START_TEST(check_parse_csq_hdr_fmt__a_custom_tag_is_used_only_when_given) {
    bcf_hdr_t *hdr = bcsq_hdr("MYCSQ");

    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, "MYCSQ");
    ck_assert_ptr_nonnull(fmt);
    ck_assert_str_eq(fmt->tag, "MYCSQ");
    destroy_hdr_fmt(fmt);

    ck_assert_ptr_null(parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, NULL)); // The default tag is not in this header
    ck_assert_ptr_null(parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, "")); // Empty is the default as well

    bcf_hdr_destroy(hdr);
} END_TEST

// A file of the other kind than the option says is refused, not read as one without transcripts
START_TEST(check_parse_csq_hdr_fmt__the_columns_of_the_other_kind_are_refused) {
    bcf_hdr_t *vep = csq_hdr();
    ck_assert_ptr_null(parse_csq_hdr_fmt(vep, ANNOTATION_BCFTOOLS, "CSQ")); // No transcript column
    bcf_hdr_destroy(vep);

    bcf_hdr_t *bcftools = bcsq_hdr("BCSQ");
    ck_assert_ptr_null(parse_csq_hdr_fmt(bcftools, ANNOTATION_VEP, "BCSQ")); // No Feature_type and Feature
    bcf_hdr_destroy(bcftools);
} END_TEST

START_TEST(check_parse_csq_hdr_fmt__an_info_field_without_a_format_is_refused) {
    bcf_hdr_t *hdr = mock_hdr();
    ck_assert_int_eq(bcf_hdr_append(hdr, "##INFO=<ID=BCSQ,Number=.,Type=String,Description=\"Consequences\">"), 0);
    ck_assert_int_eq(bcf_hdr_sync(hdr), 0);

    ck_assert_ptr_null(parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, NULL));

    bcf_hdr_destroy(hdr);
} END_TEST

// One record with every kind of entry: names of bcftools mapped to VEP's, a '*' dropped, a pointer ('@') and entries without a
// transcript left out, and the two entries of ENST1 joined in the order they come
START_TEST(check_parse_vep_transcript_consequences__bcftools_entries) {
    bcf_hdr_t *hdr = bcsq_hdr("BCSQ");
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, NULL);
    bcf1_t *v = tagged_record(hdr, "BCSQ",
        "frameshift&splice_region|G1|ENST1|protein_coding|+|10A>B|100A>T,"
        "*stop_gained|G1|ENST2|protein_coding|+|11A>*|101A>T,"
        "@874816,"
        "intergenic|||,"
        "intergenic,"
        "splice_donor|G1|ENST1|protein_coding|+|12A>B|102A>T,"
        "missense|G1|ENST3|protein_coding|+|13A>B|103A>T");

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_ptr_nonnull(csqs);
    ck_assert_int_eq(kv_size(*csqs), 3);

    ck_assert_str_eq(kv_A(*csqs, 0)->tid, "ENST1");
    ck_assert_int_eq(kv_A(*csqs, 0)->num_consequences, 3);
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[0], "frameshift_variant");
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[1], "splice_region_variant");
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[2], "splice_donor_variant");

    ck_assert_str_eq(kv_A(*csqs, 1)->tid, "ENST2");
    ck_assert_int_eq(kv_A(*csqs, 1)->num_consequences, 1);
    ck_assert_str_eq(kv_A(*csqs, 1)->consequences[0], "stop_gained");

    ck_assert_str_eq(kv_A(*csqs, 2)->tid, "ENST3");
    ck_assert_str_eq(kv_A(*csqs, 2)->consequences[0], "missense"); // A name that is not mapped is kept as it is

    for (int i = 0; i < 3; i++) ck_assert_int_eq(kv_A(*csqs, 0)->triggers[i], 0); // No trigger yet

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

START_TEST(check_parse_vep_transcript_consequences__bcftools_names_become_the_names_of_the_triggers) {
    bcf_hdr_t *hdr = bcsq_hdr("BCSQ");
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_BCFTOOLS, NULL);
    bcf1_t *v = tagged_record(hdr, "BCSQ", "splice_acceptor&frameshift&stop_gained|G1|ENST1|protein_coding|+|.|.");

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_int_eq(kv_size(*csqs), 1);
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[0], SPLICE_ACCEPTOR_VARIANT);
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[1], FRAMESHIFT_VARIANT);
    ck_assert_str_eq(kv_A(*csqs, 0)->consequences[2], STOP_GAINED_VARIANT);

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

// The tag of the format is read, not CSQ: a record with a value in both fields gives the one the format is for
START_TEST(check_parse_vep_transcript_consequences__reads_the_tag_of_the_format) {
    bcf_hdr_t *hdr = csq_hdr();
    ck_assert_int_eq(bcf_hdr_append(hdr, "##INFO=<ID=OTHER,Number=.,Type=String,Description=\"Format: Allele|Consequence|Feature_type|Feature\">"), 0);
    ck_assert_int_eq(bcf_hdr_sync(hdr), 0);
    HdrFormat *fmt = parse_csq_hdr_fmt(hdr, ANNOTATION_VEP, "OTHER");
    ck_assert_ptr_nonnull(fmt);

    bcf1_t *v = csq_record(hdr, "G|missense_variant|Transcript|FROM_CSQ");
    ck_assert_int_eq(bcf_update_info_string(hdr, v, "OTHER", "G|stop_gained|Transcript|FROM_OTHER"), 0);

    VepCsqs *csqs = parse_vep_transcript_consequences(v, fmt, hdr);

    ck_assert_int_eq(kv_size(*csqs), 1);
    ck_assert_str_eq(kv_A(*csqs, 0)->tid, "FROM_OTHER");

    destroy_vep_csqs(csqs);
    bcf_destroy(v);
    destroy_hdr_fmt(fmt);
    bcf_hdr_destroy(hdr);
} END_TEST

// Alleles_from_variant: the alleles without their shared prefix. Positions are 0-based.
static Alleles alleles_of(const char *alleles) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, alleles, 0, 1, true, -1);

    Alleles a = Alleles_from_variant(v);
    a.alt = NULL; // Points into v, which is freed below; tests that need it use their own record

    bcf_destroy(v);
    bcf_hdr_destroy(hdr);

    return a;
}

START_TEST(check_alleles_from_variant__snv_has_no_prefix) {
    const Alleles a = alleles_of("A,G");

    ck_assert_int_eq(a.pos, 100);
    ck_assert_int_eq(a.prefix_len, 0);
    ck_assert_int_eq(a.ref_len, 1);
    ck_assert_int_eq(a.alt_len, 1);
} END_TEST

START_TEST(check_alleles_from_variant__mnv_has_no_prefix) {
    const Alleles a = alleles_of("AT,GC");

    ck_assert_int_eq(a.pos, 100);
    ck_assert_int_eq(a.prefix_len, 0);
    ck_assert_int_eq(a.ref_len, 2);
    ck_assert_int_eq(a.alt_len, 2);
} END_TEST

START_TEST(check_alleles_from_variant__insertion_loses_its_anchor) {
    const Alleles a = alleles_of("A,ATT");

    ck_assert_int_eq(a.pos, 101);
    ck_assert_int_eq(a.prefix_len, 1);
    ck_assert_int_eq(a.ref_len, 0);
    ck_assert_int_eq(a.alt_len, 2);
} END_TEST

START_TEST(check_alleles_from_variant__deletion_loses_its_anchor) {
    const Alleles a = alleles_of("ATG,A");

    ck_assert_int_eq(a.pos, 101);
    ck_assert_int_eq(a.prefix_len, 1);
    ck_assert_int_eq(a.ref_len, 2);
    ck_assert_int_eq(a.alt_len, 0);
} END_TEST

START_TEST(check_alleles_from_variant__strips_a_longer_shared_prefix) {
    const Alleles deletion = alleles_of("ATGC,ATG");
    ck_assert_int_eq(deletion.pos, 103);
    ck_assert_int_eq(deletion.prefix_len, 3);
    ck_assert_int_eq(deletion.ref_len, 1);
    ck_assert_int_eq(deletion.alt_len, 0);

    const Alleles insertion = alleles_of("AT,ATT");
    ck_assert_int_eq(insertion.pos, 102);
    ck_assert_int_eq(insertion.prefix_len, 2);
    ck_assert_int_eq(insertion.ref_len, 0);
    ck_assert_int_eq(insertion.alt_len, 1);
} END_TEST

START_TEST(check_alleles_from_variant__complex_indel_keeps_the_rest) {
    const Alleles a = alleles_of("ATG,AC");

    ck_assert_int_eq(a.pos, 101);
    ck_assert_int_eq(a.prefix_len, 1);
    ck_assert_int_eq(a.ref_len, 2);
    ck_assert_int_eq(a.alt_len, 1);
} END_TEST

START_TEST(check_alleles_from_variant__reads_the_first_alt_only) {
    const Alleles a = alleles_of("A,AT,G");

    ck_assert_int_eq(a.pos, 101);
    ck_assert_int_eq(a.ref_len, 0);
    ck_assert_int_eq(a.alt_len, 1);
} END_TEST

START_TEST(check_alleles_from_variant__no_alt_reads_as_empty) {
    const Alleles a = alleles_of("A");

    ck_assert_int_eq(a.pos, 100);
    ck_assert_int_eq(a.prefix_len, 0);
    ck_assert_int_eq(a.ref_len, 1);
    ck_assert_int_eq(a.alt_len, 0);
} END_TEST

START_TEST(check_record_init__alleles_point_into_its_own_copy) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,ATT", 0, 1, true, -1);

    Record *r = Record_init(v, NULL);
    bcf_destroy(v);
    bcf_hdr_destroy(hdr);

    ck_assert_ptr_eq(r->alleles.alt, r->v->d.allele[1] + 1);
    ck_assert_str_eq(r->alleles.alt, "TT");
    ck_assert_int_eq(r->alleles.pos, 101);

    Record_destroy(r);
} END_TEST

// bcf1_variant_in_region: the load-time gene test. The gene is the bases 100 to 199 ([100, 200)). A variant belongs to it when
// a base it changes lies in it; an insertion changes none, so it belongs when both bases around it do.
static bool loaded_in_gene(const uint64_t pos, const char *ref, const char *alt) {
    bcf1_t *v = mock_bcf1_t(pos, ref, alt);
    return bcf1_variant_in_region(v, (Range) { 100, 200 });
}

START_TEST(check_bcf1_variant_in_region__substitutions) {
    ck_assert(loaded_in_gene(100, "N", "x"));
    ck_assert(loaded_in_gene(199, "N", "x"));
    ck_assert(!loaded_in_gene(99, "N", "x"));
    ck_assert(!loaded_in_gene(200, "N", "x"));
    ck_assert(loaded_in_gene(99, "NN", "xx"));  // Changes 99 and 100
    ck_assert(loaded_in_gene(199, "NN", "xx")); // Changes 199 and 200
    ck_assert(!loaded_in_gene(98, "NN", "xx"));
} END_TEST

// An insertion has no REF bases after the shared one: it is the point between the base it is attached to and the next
START_TEST(check_bcf1_variant_in_region__an_insertion_inside_the_gene_is_in) {
    ck_assert(loaded_in_gene(100, "N", "Nx"));    // Between 100 and 101
    ck_assert(loaded_in_gene(150, "N", "Nx"));
    ck_assert(loaded_in_gene(198, "N", "Nx"));    // Between 198 and 199
    ck_assert(loaded_in_gene(150, "N", "Nxyz"));  // Longer insertions too
    ck_assert(loaded_in_gene(150, "NN", "NNx"));  // With a longer shared prefix
} END_TEST

// The edges follow variant_in_region: an insertion needs both bases around it in the gene
START_TEST(check_bcf1_variant_in_region__an_insertion_at_the_gene_edge_is_out) {
    ck_assert(!loaded_in_gene(99, "N", "Nx"));    // Between 99 and 100
    ck_assert(!loaded_in_gene(199, "N", "Nx"));   // Between 199 and 200
    ck_assert(!loaded_in_gene(200, "N", "Nx"));
    ck_assert(!loaded_in_gene(98, "N", "Nx"));
} END_TEST

// Deletions and complex indels count by the bases they change, not by the base they are attached to
START_TEST(check_bcf1_variant_in_region__deletions_count_by_the_deleted_bases) {
    ck_assert(loaded_in_gene(99, "NN", "N"));      // Attached before the gene, deletes 100
    ck_assert(!loaded_in_gene(98, "NN", "N"));     // Deletes 99
    ck_assert(loaded_in_gene(198, "NN", "N"));     // Deletes 199
    ck_assert(!loaded_in_gene(199, "NN", "N"));    // Attached on the last base, deletes 200
    ck_assert(!loaded_in_gene(199, "NNN", "N"));   // Deletes 200 and 201
    ck_assert(loaded_in_gene(198, "NNN", "N"));    // Deletes 199 and 200
    ck_assert(loaded_in_gene(99, "NNN", "N"));     // Deletes 100 and 101
} END_TEST

START_TEST(check_bcf1_variant_in_region__complex_indels_count_by_the_changed_bases) {
    ck_assert(loaded_in_gene(198, "NNN", "Nx"));   // Replaces 199 and 200
    ck_assert(!loaded_in_gene(199, "NNN", "Nx"));  // Replaces 200 and 201
    ck_assert(loaded_in_gene(99, "NNN", "Nx"));    // Replaces 100 and 101
    ck_assert(!loaded_in_gene(97, "NNN", "Nx"));   // Replaces 98 and 99
    ck_assert(loaded_in_gene(98, "NNN", "Nx"));    // Replaces 99 and 100
} END_TEST

START_TEST(check_bcf1_variant_in_region__agrees_with_variant_in_region_on_the_view) {
    static const char *variants[][2] = { { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "N", "Nxyz" }, { "NN", "N" }, { "NNNN", "N" }, { "NNN", "Nx" } };

    for (int i = 0; i < 7; i++) {
        for (int pos = 90; pos <= 210; pos++) {
            bcf1_t *v = mock_bcf1_t(pos, variants[i][0], variants[i][1]);
            const Alleles a = Alleles_from_variant(v);

            ck_assert_msg(bcf1_variant_in_region(v, (Range) { 100, 200 }) == variant_in_region(&a, (Range) { 100, 200 }),
                "%s>%s at %d", variants[i][0], variants[i][1], pos);
        }
    }
} END_TEST

Suite *record_suite() {
    Suite *s = suite_create("Record");

    TCase *tc_from_exons = mocked_tcase("introns_from_exons");
    tcase_add_test(tc_from_exons, check_introns_from_exons__no_exons);
    tcase_add_test(tc_from_exons, check_introns_from_exons__single_exon_has_no_introns);
    tcase_add_test(tc_from_exons, check_introns_from_exons__two_exons);
    tcase_add_test(tc_from_exons, check_introns_from_exons__three_exons_accumulate_offset);
    tcase_add_test(tc_from_exons, check_introns_from_exons__one_intron_per_exon_gap);
    tcase_add_test(tc_from_exons, check_introns_from_exons__touching_exons_give_empty_intron);
    suite_add_tcase(s, tc_from_exons);

    TCase *tc_sites = mocked_tcase("splice_site_bases");
    tcase_add_test(tc_sites, check_exon_bases__forward_donor_is_last_base_and_acceptor_is_first);
    tcase_add_test(tc_sites, check_exon_bases__reverse_donor_is_first_base_and_acceptor_is_last);
    tcase_add_test(tc_sites, check_exon_bases__unknown_strand_behaves_as_forward);
    tcase_add_test(tc_sites, check_intron_between_sites__forward);
    tcase_add_test(tc_sites, check_intron_between_sites__reverse_gives_the_same_ascending_range);
    tcase_add_test(tc_sites, check_intron_bases__are_the_inverse_of_intron_between_sites);
    suite_add_tcase(s, tc_sites);

    TCase *tc_utr = mocked_tcase("utr_range");
    tcase_add_test(tc_utr, check_utr_range__forward);
    tcase_add_test(tc_utr, check_utr_range__reverse_is_mirrored_and_still_ascending);
    tcase_add_test(tc_utr, check_utr_range__empty_when_cds_touches_the_gene_edge);
    suite_add_tcase(s, tc_utr);

    TCase *tc_window = mocked_tcase("exon_splice_window");
    tcase_add_test(tc_window, check_exon_splice_window__forward_donor_looks_to_the_exon_above);
    tcase_add_test(tc_window, check_exon_splice_window__forward_acceptor_looks_to_the_exon_below);
    tcase_add_test(tc_window, check_exon_splice_window__reverse_sites_use_the_window_of_the_same_border);
    tcase_add_test(tc_window, check_exon_splice_window__windows_are_ascending_on_both_strands);
    tcase_add_test(tc_window, check_exon_splice_window__unknown_strand_has_no_windows);
    suite_add_tcase(s, tc_window);

    TCase *tc_phase = mocked_tcase("transcript_phase_copy");
    tcase_add_test(tc_phase, check_transcript_phase_copy__has_the_values_of_the_shared_transcript);
    tcase_add_test(tc_phase, check_transcript_phase_copy__shares_the_read_only_data_but_not_the_exons);
    tcase_add_test(tc_phase, check_transcript_phase_copy__starts_without_phase_state);
    tcase_add_test(tc_phase, check_transcript_phase_copy__phase_state_is_independent);
    tcase_add_test(tc_phase, check_transcript_phase_copy__destroy_leaves_the_shared_transcript_intact);
    tcase_add_test(tc_phase, check_transcript_phase_copy__transcript_without_exons);
    tcase_add_test(tc_phase, check_splice_sites_destroy__empty_sites);
    suite_add_tcase(s, tc_phase);

    TCase *tc_build = mocked_tcase("transcript_build");
    tcase_add_test(tc_build, check_transcript_build__converts_coordinates_to_half_open);
    tcase_add_test(tc_build, check_transcript_build__new_transcript_has_no_phase_state);
    tcase_add_test(tc_build, check_transcript_build__does_not_alias_input);
    tcase_add_test(tc_build, check_transcript_build__exon_without_cds_does_not_use_up_a_coding_region);
    tcase_add_test(tc_build, check_transcript_build__in_frame_depends_on_cds_length);
    tcase_add_test(tc_build, check_transcript_build__derives_introns_from_the_half_open_exons);
    tcase_add_test(tc_build, check_transcript_build__no_exons);
    suite_add_tcase(s, tc_build);

    TCase *tc_nmd = mocked_tcase("nmd_regions");
    tcase_add_test(tc_nmd, check_nmd_regions__forward_start_proximal_is_the_first_150_coding_bases);
    tcase_add_test(tc_nmd, check_nmd_regions__forward_long_exon_leaves_the_last_401_bases_of_an_exon_out);
    tcase_add_test(tc_nmd, check_nmd_regions__forward_distal_is_the_last_exon_the_intron_and_the_end_of_the_penultimate_exon);
    tcase_add_test(tc_nmd, check_nmd_regions__reverse_regions_are_mirrored_and_still_ascending);
    tcase_add_test(tc_nmd, check_nmd_regions__an_exon_of_401_bases_has_no_long_exon_region_and_one_of_402_has_one_base);
    tcase_add_test(tc_nmd, check_nmd_regions__a_single_exon_transcript_is_all_distal);
    tcase_add_test(tc_nmd, check_nmd_regions__a_short_penultimate_exon_is_distal_as_a_whole_and_a_short_cds_is_start_proximal_as_a_whole);
    tcase_add_test(tc_nmd, check_nmd_regions__a_transcript_without_exons_has_no_regions);
    tcase_add_test(tc_nmd, check_nmd_regions__an_unknown_strand_is_read_as_forward);
    suite_add_tcase(s, tc_nmd);

    TCase *tc_lifecycle = mocked_tcase("record_lifecycle");
    tcase_add_test(tc_lifecycle, check_record_init__starts_without_flags_or_isoforms);
    tcase_add_test(tc_lifecycle, check_isoforms_destroy__frees_isoforms_and_rescue_probabilities);
    tcase_add_test(tc_lifecycle, check_isoforms_destroy__no_isoforms);
    tcase_add_test(tc_lifecycle, check_record_destroy__frees_flags_and_built_isoform_slots);
    tcase_add_test(tc_lifecycle, check_record_destroy__frees_the_scores);
    suite_add_tcase(s, tc_lifecycle);

    TCase *tc_alleles = mocked_tcase("alleles");
    tcase_add_test(tc_alleles, check_alleles_from_variant__snv_has_no_prefix);
    tcase_add_test(tc_alleles, check_alleles_from_variant__mnv_has_no_prefix);
    tcase_add_test(tc_alleles, check_alleles_from_variant__insertion_loses_its_anchor);
    tcase_add_test(tc_alleles, check_alleles_from_variant__deletion_loses_its_anchor);
    tcase_add_test(tc_alleles, check_alleles_from_variant__strips_a_longer_shared_prefix);
    tcase_add_test(tc_alleles, check_alleles_from_variant__complex_indel_keeps_the_rest);
    tcase_add_test(tc_alleles, check_alleles_from_variant__reads_the_first_alt_only);
    tcase_add_test(tc_alleles, check_alleles_from_variant__no_alt_reads_as_empty);
    tcase_add_test(tc_alleles, check_record_init__alleles_point_into_its_own_copy);
    suite_add_tcase(s, tc_alleles);

    TCase *tc_vep = mocked_tcase("vep");
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__the_position_of_every_field);
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__a_header_without_csq_gives_null);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__one_entry_per_transcript);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__a_record_without_csq_has_no_entries);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__only_other_features_has_no_entries);
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__bcftools_uses_bcsq_and_its_columns);
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__a_custom_tag_is_used_only_when_given);
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__the_columns_of_the_other_kind_are_refused);
    tcase_add_test(tc_vep, check_parse_csq_hdr_fmt__an_info_field_without_a_format_is_refused);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__bcftools_entries);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__bcftools_names_become_the_names_of_the_triggers);
    tcase_add_test(tc_vep, check_parse_vep_transcript_consequences__reads_the_tag_of_the_format);
    tcase_add_test(tc_vep, check_record_init__an_empty_annotation_is_kept_and_freed);
    suite_add_tcase(s, tc_vep);

    TCase *tc_load = mocked_tcase("bcf1_variant_in_region");
    tcase_add_test(tc_load, check_bcf1_variant_in_region__substitutions);
    tcase_add_test(tc_load, check_bcf1_variant_in_region__an_insertion_inside_the_gene_is_in);
    tcase_add_test(tc_load, check_bcf1_variant_in_region__an_insertion_at_the_gene_edge_is_out);
    tcase_add_test(tc_load, check_bcf1_variant_in_region__deletions_count_by_the_deleted_bases);
    tcase_add_test(tc_load, check_bcf1_variant_in_region__complex_indels_count_by_the_changed_bases);
    tcase_add_test(tc_load, check_bcf1_variant_in_region__agrees_with_variant_in_region_on_the_view);
    suite_add_tcase(s, tc_load);

    return s;
}
