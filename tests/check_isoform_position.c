#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cpliceai/range.h>

#include "klib/kvec.h"
#include "haplotype.h"
#include "isoform_position.h"
#include "record.h"
#include "utils.h"

#include "check_isoform_position.h"
#include "mocks.h"

/*
* is_utr_variant as it was before it counted only changed bases, kept to check that the new version flags a subset
*/
static enum UTRType old_is_utr_variant(const Record *r, const NonCodingRegions *ncrs, const Records *rs, const Range bounds) {
    const int64_t position = r->v->pos;
    const int64_t isoform_position = position - haplotype_offset(rs, (Range) { bounds.start, position-1 });
    const Range variant = { isoform_position, isoform_position + r->v->rlen };

    for (int i = 0; i < ncrs->n; i++) {
        NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_utr) continue;

        if (ranges_overlap(ncr.coords, variant)) return ncr.is_utr;
    }

    return 0;
}

/*
* A gene on [100, 200) with a 5' UTR, an intron and a 3' UTR, in transcript order for the strand
*/
static const Range iso_gene = { 100, 200 };

static NonCodingRegion iso_region(const int is_used, const Range coords, const enum UTRType utr) {
    NonCodingRegion ncr = {0};
    ncr.is_used = is_used;
    ncr.coords = coords;
    ncr.canonical_coords = coords;
    ncr.is_utr = utr;
    return ncr;
}

// The non-coding regions of a gene [100, 200): a UTR at each end and up to two introns between them (ascending), given in
// transcript order for the strand
typedef struct RegionLayout {
    const char *name;
    Range low_utr, high_utr;
    Range introns[2];
    int n_introns;
} RegionLayout;

static NonCodingRegions layout_regions(const RegionLayout *layout, const int strand, const int introns_used) {
    NonCodingRegions ncrs; kv_init(ncrs);
    const NonCodingRegion low = iso_region(1, layout->low_utr, strand == STRAND_FWD ? UTR5 : UTR3);
    const NonCodingRegion high = iso_region(1, layout->high_utr, strand == STRAND_FWD ? UTR3 : UTR5);

    if (strand == STRAND_FWD) {
        kv_push(NonCodingRegion, ncrs, low);
        for (int i = 0; i < layout->n_introns; i++) kv_push(NonCodingRegion, ncrs, iso_region(introns_used, layout->introns[i], 0));
        kv_push(NonCodingRegion, ncrs, high);
    } else {
        kv_push(NonCodingRegion, ncrs, high);
        for (int i = layout->n_introns - 1; i >= 0; i--) kv_push(NonCodingRegion, ncrs, iso_region(introns_used, layout->introns[i], 0));
        kv_push(NonCodingRegion, ncrs, low);
    }
    return ncrs;
}

// The layout every other test uses: UTRs [100, 110) and [180, 200), one intron [130, 150)
static const RegionLayout default_layout = { "default", { 100, 110 }, { 180, 200 }, { { 130, 150 } }, 1 };

static NonCodingRegions iso_regions(const int strand, const int intron_used) {
    return layout_regions(&default_layout, strand, intron_used);
}

static NonCodingRegions to_haplotype(const NonCodingRegions *ncrs, const Records *rs) {
    NonCodingRegions hap; kv_init(hap);
    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        NonCodingRegion ncr = kv_A(*ncrs, i);
        ncr.coords = region_to_haplotype(rs, iso_gene, ncr.coords, ncr.is_utr != 0);
        ncr.canonical_coords = region_to_haplotype(rs, iso_gene, ncr.canonical_coords, ncr.is_utr != 0);
        kv_push(NonCodingRegion, hap, ncr);
    }
    return hap;
}

#define N_ISO_VARIANTS 8
static char *iso_variants[N_ISO_VARIANTS][2] = {
    { "N", "x" }, { "NN", "xy" }, { "N", "Nx" }, { "N", "Nxyz" }, { "NN", "N" }, { "NNNN", "N" }, { "NNN", "Nx" }, { "NN", "Nxyz" }
};

// Every variant kind at every position across the gene, on both strands, with the intron spliced or retained, alone and
// after a deletion earlier in the gene (so the haplotype offsets are not 0)
START_TEST(check_isoform_position__classification_follows_the_coding_effect) {
    const int strands[2] = { STRAND_FWD, STRAND_REV };

    for (int si = 0; si < 2; si++) {
        for (int intron_used = 0; intron_used <= 1; intron_used++) {
            NonCodingRegions ncrs = iso_regions(strands[si], intron_used);

            for (int vi = 0; vi < N_ISO_VARIANTS; vi++) {
                for (int pos = 98; pos <= 201; pos++) {
                    for (int with_earlier = 0; with_earlier <= 1; with_earlier++) {
                        if (with_earlier && pos <= 116) continue;

                        Records rs; kv_init(rs);
                        if (with_earlier) kv_push(Record *, rs, Record_init(mock_bcf1_t(112, "NNNN", "N"), NULL));
                        Record *r = Record_init(mock_bcf1_t(pos, iso_variants[vi][0], iso_variants[vi][1]), NULL);
                        kv_push(Record *, rs, r);

                        char what[96];
                        snprintf(what, sizeof(what), "%s>%s at %d, strand %d, intron %s%s", iso_variants[vi][0], iso_variants[vi][1], pos, strands[si], intron_used ? "spliced" : "retained", with_earlier ? ", after a deletion" : "");

                        // The classification and the offset follow the coding effect of the variant (the pairing rules)
                        int64_t position = -1;
                        int offset = 0;
                        const int noncoding = get_range_in_isoform(&r->alleles, &ncrs, &rs, strands[si], iso_gene, &position, &offset);
                        const CodingEffect ce = variant_coding_effect(&r->alleles, &ncrs, iso_gene);
                        if (variant_in_region(&r->alleles, iso_gene)) {
                            ck_assert_msg(noncoding == !ce.coding, "%s: non-coding %d, coding effect %d", what, noncoding, ce.coding);
                            if (!noncoding) ck_assert_msg(offset == ce.deleted - ce.inserted, "%s: offset %d, expected %d", what, offset, ce.deleted - ce.inserted);
                        }

                        // Deliberate difference: the old version flagged a variant that only touched a UTR (an SNV on the base
                        // next to it) and counted the attached base of an indel. A variant with REF bases is flagged now only
                        // when a changed base lies in the UTR, so it is flagged less often than before, never differently
                        // (insertions are covered by their own test). The old version also moved the variant by the deletions
                        // before it, so it is only comparable when there is none.
                        const enum UTRType new_utr = is_utr_variant(&r->alleles, &ncrs);
                        const enum UTRType old_utr = old_is_utr_variant(r, &ncrs, &rs, iso_gene);
                        if (r->alleles.ref_len > 0 && !with_earlier) ck_assert_msg(new_utr == 0 || new_utr == old_utr, "%s: UTR %d, old %d", what, new_utr, old_utr);

                        for (size_t k = 0; k < kv_size(rs); k++) Record_destroy(kv_A(rs, k));
                        kv_destroy(rs);
                    }
                }
            }

            kv_destroy(ncrs);
        }
    }
} END_TEST

// An insertion exactly at the border of a UTR does not split the codon that bounds the CDS, so it belongs to the UTR: it is
// non-coding on both strands. On the reverse strand it used to fall out of the region loop and count as coding, both at the CDS
// side of a 3' UTR (after the stop codon) and at the far end of the 5' UTR. Deletions and substitutions that reach the CDS
// stay coding, and an insertion at an exon/intron border stays coding (the splice site lies in the intron).
static int range_in_isoform_of(const int strand, const int pos, const char *ref, const char *alt, int64_t *position) {
    NonCodingRegions ncrs = iso_regions(strand, 1);
    Records rs; kv_init(rs);
    Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
    kv_push(Record *, rs, r);

    int offset = 0;
    const int noncoding = get_range_in_isoform(&r->alleles, &ncrs, &rs, strand, iso_gene, position, &offset);

    kv_destroy(ncrs);
    Record_destroy(r); kv_destroy(rs);
    return noncoding;
}

START_TEST(check_get_range_in_isoform__an_insertion_at_the_edge_of_a_utr_is_noncoding) {
    int64_t position;
    const int strands[2] = { STRAND_FWD, STRAND_REV };

    for (int si = 0; si < 2; si++) {
        // After 109: next to the CDS (low UTR); after 199: the far end of the high UTR
        for (int edge = 0; edge < 2; edge++) {
            const int pos = edge ? 199 : 109;
            ck_assert_msg(range_in_isoform_of(strands[si], pos, "N", "Nx", &position) == 1, "1 base at %d, strand %d", pos, strands[si]);
            ck_assert_msg(range_in_isoform_of(strands[si], pos, "N", "Nxyz", &position) == 1, "3 bases at %d, strand %d", pos, strands[si]);
        }

        // After 179: between the CDS and the high UTR (already non-coding on both strands)
        ck_assert_int_eq(range_in_isoform_of(strands[si], 179, "N", "Nx", &position), 1);
    }

    // The position of a non-coding variant is its reference position on the reverse strand
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 109, "N", "Nx", &position), 1);
    ck_assert_int_eq(position, 109);
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 199, "N", "Nx", &position), 1);
    ck_assert_int_eq(position, 199);
} END_TEST

START_TEST(check_get_range_in_isoform__only_pure_insertions_change_at_the_edge_of_a_utr) {
    int64_t position;

    // Deleting base 110, the CDS base next to the low UTR, or replacing it, still changes the CDS
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 109, "NN", "N", &position), 0); // Deletion
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 109, "NN", "Nxyz", &position), 0); // Replaces base 110 by 3 bases
    ck_assert_int_eq(range_in_isoform_of(STRAND_FWD, 109, "NN", "N", &position), 0);

    // An insertion at an exon/intron border lies in the exon, on both strands
    ck_assert_int_eq(range_in_isoform_of(STRAND_FWD, 129, "N", "Nx", &position), 0);
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 129, "N", "Nx", &position), 0);
    ck_assert_int_eq(range_in_isoform_of(STRAND_FWD, 149, "N", "Nx", &position), 0);
    ck_assert_int_eq(range_in_isoform_of(STRAND_REV, 149, "N", "Nx", &position), 0);
} END_TEST

// is_utr_variant flags a variant when a base it changes lies in a UTR. A variant that only touches the UTR (the SNV on the first
// base after it) is not flagged, and neither is the attached base of an indel. An insertion counts as in the UTR at either edge,
// as get_range_in_isoform treats it. UTRs of the layout: [100, 110) and [180, 200); forward strand: 5' then 3'.
static enum UTRType utr_of(const int strand, const int pos, const char *ref, const char *alt) {
    NonCodingRegions ncrs = iso_regions(strand, 1);
    Records rs; kv_init(rs);
    Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
    kv_push(Record *, rs, r);

    const enum UTRType utr = is_utr_variant(&r->alleles, &ncrs);

    kv_destroy(ncrs);
    Record_destroy(r); kv_destroy(rs);
    return utr;
}

// The UTRs are reference ranges, so a variant is compared with them at its reference position, whatever else the haplotype
// holds. It used to be moved into haplotype coordinates first (minus the bases deleted before it, counted up to the base before
// it), so a variant after a long deletion near the UTR could land in it: a deletion of 106 to 120 moved an SNV at 121, a base of
// the CDS, to 107, in the 5' UTR.
START_TEST(check_isoform_base_is_coding__follows_the_used_regions) {
    NonCodingRegions spliced = iso_regions(STRAND_FWD, 1), retained = iso_regions(STRAND_FWD, 0);

    ck_assert(!isoform_base_is_coding(&spliced, 105));   // 5' UTR
    ck_assert(isoform_base_is_coding(&spliced, 110));
    ck_assert(isoform_base_is_coding(&spliced, 129));
    ck_assert(!isoform_base_is_coding(&spliced, 130));   // Intron
    ck_assert(!isoform_base_is_coding(&spliced, 149));
    ck_assert(isoform_base_is_coding(&spliced, 150));
    ck_assert(!isoform_base_is_coding(&spliced, 180));   // 3' UTR
    ck_assert(isoform_base_is_coding(&retained, 140));   // A retained intron is part of the coding sequence

    kv_destroy(spliced); kv_destroy(retained);
} END_TEST

START_TEST(check_is_utr_variant__only_the_reference_position_counts) {
    const int strands[2] = { STRAND_FWD, STRAND_REV };

    for (int si = 0; si < 2; si++) {
        NonCodingRegions ncrs = iso_regions(strands[si], 1);

        Record *snv = Record_init(mock_bcf1_t(121, "N", "x"), NULL);
        ck_assert_int_eq(is_utr_variant(&snv->alleles, &ncrs), 0);

        Record *in_utr = Record_init(mock_bcf1_t(105, "N", "x"), NULL);
        ck_assert_int_ne(is_utr_variant(&in_utr->alleles, &ncrs), 0);

        Record_destroy(snv); Record_destroy(in_utr);
        kv_destroy(ncrs);
    }
} END_TEST

START_TEST(check_is_utr_variant__substitutions_need_a_base_in_the_utr) {
    ck_assert_int_eq(utr_of(STRAND_FWD, 100, "N", "x"), UTR5);
    ck_assert_int_eq(utr_of(STRAND_FWD, 109, "N", "x"), UTR5);
    ck_assert_int_eq(utr_of(STRAND_FWD, 110, "N", "x"), 0);   // The first CDS base
    ck_assert_int_eq(utr_of(STRAND_FWD, 179, "N", "x"), 0);   // The last CDS base
    ck_assert_int_eq(utr_of(STRAND_FWD, 180, "N", "x"), UTR3);
    ck_assert_int_eq(utr_of(STRAND_FWD, 199, "N", "x"), UTR3);

    ck_assert_int_eq(utr_of(STRAND_FWD, 108, "NN", "xy"), UTR5);
    ck_assert_int_eq(utr_of(STRAND_FWD, 109, "NN", "xy"), UTR5); // 109 is in the UTR, 110 is not
    ck_assert_int_eq(utr_of(STRAND_FWD, 110, "NN", "xy"), 0);
    ck_assert_int_eq(utr_of(STRAND_FWD, 178, "NN", "xy"), 0);
    ck_assert_int_eq(utr_of(STRAND_FWD, 179, "NN", "xy"), UTR3); // 180 is in the UTR

    // The low UTR is the 3' UTR on the reverse strand
    ck_assert_int_eq(utr_of(STRAND_REV, 109, "N", "x"), UTR3);
    ck_assert_int_eq(utr_of(STRAND_REV, 110, "N", "x"), 0);
    ck_assert_int_eq(utr_of(STRAND_REV, 180, "N", "x"), UTR5);
    ck_assert_int_eq(utr_of(STRAND_REV, 179, "N", "x"), 0);
} END_TEST

START_TEST(check_is_utr_variant__indels_count_by_the_changed_bases) {
    // NN>N at p deletes p + 1
    ck_assert_int_eq(utr_of(STRAND_FWD, 108, "NN", "N"), UTR5);  // Deletes 109
    ck_assert_int_eq(utr_of(STRAND_FWD, 109, "NN", "N"), 0);     // Deletes 110, the first CDS base; the anchor 109 is in the UTR
    ck_assert_int_eq(utr_of(STRAND_FWD, 178, "NN", "N"), 0);     // Deletes 179
    ck_assert_int_eq(utr_of(STRAND_FWD, 179, "NN", "N"), UTR3);  // Deletes 180; the anchor 179 is a CDS base
    ck_assert_int_eq(utr_of(STRAND_FWD, 107, "NNNN", "N"), UTR5); // Deletes 108 to 110
    ck_assert_int_eq(utr_of(STRAND_FWD, 110, "NNNN", "N"), 0);

    // A complex indel counts by the REF bases it replaces
    ck_assert_int_eq(utr_of(STRAND_FWD, 109, "NN", "Nxyz"), 0);    // Replaces 110
    ck_assert_int_eq(utr_of(STRAND_FWD, 108, "NN", "Nxyz"), UTR5); // Replaces 109
} END_TEST

START_TEST(check_is_utr_variant__an_insertion_is_in_the_utr_at_either_edge) {
    // N>Nx at p inserts between p and p + 1
    ck_assert_int_eq(utr_of(STRAND_FWD, 100, "N", "Nx"), UTR5);
    ck_assert_int_eq(utr_of(STRAND_FWD, 109, "N", "Nx"), UTR5);  // Between the last UTR base and the first CDS base
    ck_assert_int_eq(utr_of(STRAND_FWD, 110, "N", "Nx"), 0);     // Inside the CDS
    ck_assert_int_eq(utr_of(STRAND_FWD, 178, "N", "Nx"), 0);
    ck_assert_int_eq(utr_of(STRAND_FWD, 179, "N", "Nx"), UTR3);  // Between the last CDS base and the first UTR base
    ck_assert_int_eq(utr_of(STRAND_FWD, 199, "N", "Nx"), UTR3);  // After the last gene base
    ck_assert_int_eq(utr_of(STRAND_FWD, 200, "N", "Nx"), 0);     // Beyond the gene: not in a UTR (and not processed)
    ck_assert_int_eq(utr_of(STRAND_REV, 109, "N", "Nx"), UTR3);
    ck_assert_int_eq(utr_of(STRAND_REV, 179, "N", "Nx"), UTR5);
} END_TEST

// variant_in_cds_exon decides whether a record counts as being in the coding sequence of an exon, which is what the
// NON_CONSTITUTIVE and INFRAME flags are about. A variant with REF bases counts at its first changed base. An insertion lies
// between two bases: inside the CDS it is in; at an edge of the CDS it follows the border rule. Where the CDS ends at a UTR
// (the exon goes on past it) the insertion is not in; where the CDS ends at the end of the exon, at an intron, it is.
//   exon A: exon [100, 150), CDS [110, 150): a UTR before the CDS, an intron after it
//   exon B: exon [200, 260), CDS [200, 250): an intron before the CDS, a UTR after it
static bool in_cds_exon(const Range exon_bounds, const Range cds, const int pos, const char *ref, const char *alt) {
    Range bounds = exon_bounds, coding = cds;
    const Exon exon = { .cds = &coding, .ref_bounds = &bounds, .splice_sites = NULL, .constitutive = true, .in_frame = true };
    const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
    return variant_in_cds_exon(&a, &exon);
}

static bool in_exon_a(const int pos, const char *ref, const char *alt) { return in_cds_exon((Range) { 100, 150 }, (Range) { 110, 150 }, pos, ref, alt); }
static bool in_exon_b(const int pos, const char *ref, const char *alt) { return in_cds_exon((Range) { 200, 260 }, (Range) { 200, 250 }, pos, ref, alt); }

START_TEST(check_variant_in_cds_exon__substitutions_and_deletions_count_at_their_first_changed_base) {
    ck_assert(!in_exon_a(109, "N", "x"));
    ck_assert(in_exon_a(110, "N", "x"));
    ck_assert(in_exon_a(149, "N", "x"));
    ck_assert(!in_exon_a(150, "N", "x"));

    // NN>N at p deletes p + 1
    ck_assert(!in_exon_a(108, "NN", "N")); // Deletes 109, a UTR base
    ck_assert(in_exon_a(109, "NN", "N"));  // Deletes 110, the first CDS base; the attached base 109 is not in the CDS
    ck_assert(in_exon_a(148, "NN", "N"));  // Deletes 149
    ck_assert(!in_exon_a(149, "NN", "N")); // Deletes 150, past the CDS; the attached base 149 is in it

    // Unchanged: a deletion that starts in the UTR and runs into the CDS counts at its first base, so it is not in
    ck_assert(!in_exon_a(107, "NNNN", "N")); // Deletes 108 to 110
} END_TEST

START_TEST(check_variant_in_cds_exon__an_insertion_inside_the_cds_is_in) {
    ck_assert(in_exon_a(110, "N", "Nx"));  // Between 110 and 111
    ck_assert(in_exon_a(148, "N", "Nx"));  // Between 148 and 149
    ck_assert(in_exon_b(210, "N", "Nxyz"));
} END_TEST

START_TEST(check_variant_in_cds_exon__an_insertion_at_the_border_with_a_utr_is_not_in) {
    ck_assert(!in_exon_a(109, "N", "Nx")); // Between 109 (UTR) and 110 (CDS): the start codon stays whole
    ck_assert(!in_exon_b(249, "N", "Nx")); // Between 249 (CDS) and 250 (UTR): the stop codon stays whole
    ck_assert(!in_exon_b(249, "N", "Nxyz"));
    ck_assert(!in_exon_a(99, "N", "Nx"));  // Outside the exon
    ck_assert(!in_exon_b(260, "N", "Nx"));
} END_TEST

START_TEST(check_variant_in_cds_exon__an_insertion_at_the_border_with_an_intron_is_in) {
    ck_assert(in_exon_a(149, "N", "Nx"));  // Between 149 (the last base of the exon) and 150 (intron): the splice site is in the intron
    ck_assert(in_exon_b(199, "N", "Nx"));  // Between 199 (intron) and 200 (the first base of the exon)
    ck_assert(!in_exon_a(150, "N", "Nx")); // Inside the intron
    ck_assert(!in_exon_b(198, "N", "Nx"));
} END_TEST

// variant_coding_effect splits a variant by what it does to the coding sequence, with the rules agreed with the user:
//   - the REF and ALT bases are paired one to one from the left in genomic coordinates (whatever the strand); the paired bases
//     are substitutions, the REF bases left over are a deletion at the right-hand end, the ALT bases left over are an
//     insertion after the last paired base
//   - a base is coding when it is in the gene and in no used non-coding region; the coding part of a variant is applied and
//     the non-coding part is not (bases outside the gene are non-coding)
//   - an insertion point is coding strictly inside the CDS and at an intron border; not at a UTR border or outside the gene
// Layout: gene [100, 200), UTRs [100, 110) and [180, 200), CDS [110, 130) and [150, 180), intron [130, 150).
static void fill_n(char *buf, const int n) {
    memset(buf, 'N', n);
    buf[n] = '\0';
}

static void assert_effect(const int pos, const char *ref, const char *alt, const bool coding, const int substituted, const int deleted, const int inserted) {
    const int strands[2] = { STRAND_FWD, STRAND_REV };

    for (int si = 0; si < 2; si++) {
        NonCodingRegions ncrs = iso_regions(strands[si], 1);
        const Alleles a = Alleles_from_variant(mock_bcf1_t(pos, ref, alt));
        const CodingEffect e = variant_coding_effect(&a, &ncrs, iso_gene);

        ck_assert_msg(e.coding == coding && e.substituted == substituted && e.deleted == deleted && e.inserted == inserted,
            "%.8s>%.8s at %d, strand %d: coding %d, substituted %d, deleted %d, inserted %d; expected %d, %d, %d, %d",
            ref, alt, pos, strands[si], e.coding, e.substituted, e.deleted, e.inserted, coding, substituted, deleted, inserted);

        kv_destroy(ncrs);
    }
}

START_TEST(check_variant_coding_effect__substitutions) {
    assert_effect(109, "N", "x", false, 0, 0, 0);   // UTR
    assert_effect(110, "N", "x", true, 1, 0, 0);
    assert_effect(129, "N", "x", true, 1, 0, 0);
    assert_effect(130, "N", "x", false, 0, 0, 0);   // Intron
    assert_effect(179, "N", "x", true, 1, 0, 0);
    assert_effect(180, "N", "x", false, 0, 0, 0);   // UTR
    assert_effect(109, "NN", "xy", true, 1, 0, 0);  // 109 is in the UTR, 110 is coding
} END_TEST

START_TEST(check_variant_coding_effect__a_deletion_applies_only_its_coding_part) {
    char ref[32];

    // NNNN>N at p deletes p + 1 to p + 3
    assert_effect(106, "NNNN", "N", false, 0, 0, 0);  // 107 to 109: the UTR
    assert_effect(107, "NNNN", "N", true, 0, 1, 0);   // 108, 109 in the UTR, 110 coding
    assert_effect(108, "NNNN", "N", true, 0, 2, 0);
    assert_effect(109, "NNNN", "N", true, 0, 3, 0);   // 110 to 112
    assert_effect(127, "NNNN", "N", true, 0, 2, 0);   // 128, 129 coding, 130 in the intron
    assert_effect(128, "NNNN", "N", true, 0, 1, 0);
    assert_effect(129, "NNNN", "N", false, 0, 0, 0);  // 130 to 132: the intron
    assert_effect(147, "NNNN", "N", true, 0, 1, 0);   // 148, 149 in the intron, 150 coding
    assert_effect(148, "NNNN", "N", true, 0, 2, 0);
    assert_effect(149, "NNNN", "N", true, 0, 3, 0);   // 150 to 152
    assert_effect(177, "NNNN", "N", true, 0, 2, 0);   // 178, 179 coding, 180 in the UTR
    assert_effect(178, "NNNN", "N", true, 0, 1, 0);
    assert_effect(179, "NNNN", "N", false, 0, 0, 0);  // 180 to 182: the UTR

    // Deletions that run past the gene: the bases outside the gene are not coding, the coding part still applies
    fill_n(ref, 15);
    assert_effect(97, ref, "N", true, 0, 2, 0);       // 98 to 111: 98, 99 outside the gene, 100 to 109 the UTR, 110 and 111 coding
    fill_n(ref, 25);
    assert_effect(177, ref, "N", true, 0, 2, 0);      // 178 to 201: 178, 179 coding, the UTR and past the gene end do not count

    // A deletion over several regions: the coding bases of every exon it reaches
    assert_effect(127, ref, "N", true, 0, 4, 0);      // 128 to 151: 128, 129 and 150, 151 coding, the intron in between does not count
} END_TEST

START_TEST(check_variant_coding_effect__complex_indels_pair_from_the_left) {
    // NNN>Nxy at 108 replaces 109 (UTR) and 110 (coding): x pairs with the UTR base, only y is coding
    assert_effect(108, "NNN", "Nxy", true, 1, 0, 0);

    // NNN>Nx at 108: x pairs with 109, an SNV in the UTR, and 110 is deleted
    assert_effect(108, "NNN", "Nx", true, 0, 1, 0);

    // NNN>Nxyz at 108: x with 109, y with 110, and z is inserted after 110, in the CDS
    assert_effect(108, "NNN", "Nxyz", true, 1, 0, 1);

    // Mirrored, 179 coding and 180 in the UTR: pairing is fixed to genomic left to right, so x pairs with the coding base and
    // the UTR base is deleted
    assert_effect(178, "NNN", "Nx", true, 1, 0, 0);

    // NN>Nxy at 108: x pairs with 109, y is inserted after it, at the UTR border
    assert_effect(108, "NN", "Nxy", false, 0, 0, 0);

    // ... and one base on, in the CDS
    assert_effect(110, "NN", "Nxy", true, 1, 0, 1);
} END_TEST

START_TEST(check_variant_coding_effect__insertions_follow_the_border_rule) {
    // N>Nx at p inserts between p and p + 1
    assert_effect(109, "N", "Nx", false, 0, 0, 0);    // Between the UTR and the CDS
    assert_effect(110, "N", "Nx", true, 0, 0, 1);
    assert_effect(129, "N", "Nx", true, 0, 0, 1);     // Between the last exon base and the intron: the splice site is in the intron
    assert_effect(130, "N", "Nx", false, 0, 0, 0);    // Inside the intron
    assert_effect(148, "N", "Nx", false, 0, 0, 0);
    assert_effect(149, "N", "Nx", true, 0, 0, 1);     // Between the intron and the first exon base
    assert_effect(178, "N", "Nx", true, 0, 0, 1);
    assert_effect(179, "N", "Nx", false, 0, 0, 0);    // Between the CDS and the UTR
    assert_effect(199, "N", "Nx", false, 0, 0, 0);    // At the gene end
    assert_effect(110, "N", "Nxyz", true, 0, 0, 3);
} END_TEST

// isoform_coding_index counts the coding bases before a base in reading direction. Layout: exon 1 [110, 130) and exon 2
// [150, 180), 50 coding bases; the forward index of a base is the number of coding bases below it, the reverse index the number
// above it. Only the coding part of a variant counts.
static int64_t coding_index_of(const int strand, const int used, const int64_t x, const char *ref, const char *alt, const int pos) {
    NonCodingRegions ncrs = iso_regions(strand, used);
    Records rs; kv_init(rs);
    Record *r = NULL;
    if (ref != NULL) {
        r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
        kv_push(Record *, rs, r);
    }

    const int64_t index = isoform_coding_index(&rs, &ncrs, iso_gene, strand, x);

    kv_destroy(ncrs);
    if (r != NULL) Record_destroy(r);
    kv_destroy(rs);
    return index;
}

START_TEST(check_isoform_coding_index__without_variants) {
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 100, NULL, NULL, 0), 0);   // In the UTR: no coding base before it
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 110, NULL, NULL, 0), 0);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 129, NULL, NULL, 0), 19);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 135, NULL, NULL, 0), 20);  // In the intron: the next base after it is 150
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 150, NULL, NULL, 0), 20);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 179, NULL, NULL, 0), 49);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 199, NULL, NULL, 0), 50);

    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 179, NULL, NULL, 0), 0);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 150, NULL, NULL, 0), 29);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 135, NULL, NULL, 0), 30);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 129, NULL, NULL, 0), 30);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 110, NULL, NULL, 0), 49);

    // A retained intron is part of the coding sequence
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 0, 150, NULL, NULL, 0), 40);
} END_TEST

START_TEST(check_isoform_coding_index__a_variant_in_an_intron_changes_nothing) {
    // N>Nx at 135 inserts inside the intron; NNNN>N at 133 deletes 134 to 136
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 150, "N", "Nx", 135), 20);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 129, "N", "Nx", 135), 30);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 150, "NNNN", "N", 133), 20);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 129, "NNNN", "N", 133), 30);
} END_TEST

START_TEST(check_isoform_coding_index__a_coding_variant_moves_what_comes_after_it) {
    // NN>N at 120 deletes 121
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 150, "NN", "N", 120), 19);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 121, "NN", "N", 120), 11);  // The deleted base: the index of the next one
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 122, "NN", "N", 120), 11);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 110, "NN", "N", 120), 48);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 121, "NN", "N", 120), 38);

    // N>Nxy at 115 inserts two bases between 115 and 116
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 115, "N", "Nxy", 115), 5);
    ck_assert_int_eq(coding_index_of(STRAND_FWD, 1, 116, "N", "Nxy", 115), 8);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 116, "N", "Nxy", 115), 43);
    ck_assert_int_eq(coding_index_of(STRAND_REV, 1, 110, "N", "Nxy", 115), 51);
} END_TEST

/*
* Boundary table: not a test of expected values but a dump of what the code does for every variant kind around every region
* border, to decide the boundary rules from. Writes $BOUNDARY_TABLE (default boundary_table.txt, in the working directory of
* the test run). Layout of the gene as in iso_regions: UTR [100, 110), intron [130, 150), UTR [180, 200), gene [100, 200).
*
* Columns, per variant and region:
*   A  the footprint (REF bases plus the attachment base) overlaps the region   (membership today: variant_overlaps)
*   B  a changed REF base lies in the region                                    (-  for an insertion: it has none)
*   E  insertions only: both neighbours of the insertion point lie in the region (strictly inside)
*   C  the attachment base (footprint start) lies in the region                 (classification before option B)
*   D  variant_position lies in the region                                      (classification with option B)
*   '!' marks rows where A and B differ, '*' rows where C and D differ.
* Then get_range_in_isoform (nc = 1 non-coding; p = isoform position; off = coding offset) and is_utr_variant, for both
* strands with the intron spliced out.
*
* Then, for each strand, the position of a coding indel as evaluate_isoform gets it: "seq_pos->idx/len", where seq_pos is p and
* idx the same value (the position is the index into the spliced coding sequence that the stop codon scan and the pairing
* window start from), and len is the length of that sequence (built from a gene of N bases, so only its length matters). An
* idx below 0 or at or above len lies outside the coding sequence. "-" when the variant is not in the gene (it is not
* processed), or is non-coding.
* "len check f r": for every variant in the gene (also non-coding ones), the length of the spliced coding sequence on each
* strand compared with 50 - deleted + inserted from variant_coding_effect (the layout has 50 coding bases): "ok" when they
* agree, else the difference (actual - expected). A difference means the haplotype coordinates of the regions do not put the
* variant where the coding classification says it is. "-" when the variant is not in the gene.
* "ce: c s d i": variant_coding_effect (coding, substituted, deleted, inserted), which follows the rules agreed with the
* user and is not used by the pipeline yet, to compare with off.
*/
static char dump_reference[101];

// Mirrors what evaluate_isoform does with the indel shift of a record that is in the gene
static bool spliced_index(Record *r, Records *rs, const NonCodingRegions *ncrs, const int strand, int64_t *seq_pos, int64_t *index, int64_t *spliced_len) {
    if (!variant_in_region(&r->alleles, iso_gene)) return false;

    NonCodingRegions hap_ncrs = to_haplotype(ncrs, rs);
    int64_t position = 0;
    int offset = 0;
    const int noncoding = get_range_in_isoform(&r->alleles, ncrs, rs, strand, iso_gene, &position, &offset);
    if (noncoding) {
        kv_destroy(hap_ncrs);
        return false;
    }

    // The position is the index into the coding sequence
    *seq_pos = position;
    *index = position;

    memset(dump_reference, 'N', 100);
    dump_reference[100] = '\0';
    const kstring_t reference = { 100, 101, dump_reference };

    kstring_t hap = {0}, spliced = {0};
    get_haplotype_sequence(reference, rs, iso_gene, &hap);
    splice_haplotype_sequence(hap.s != NULL ? hap.s : "", &hap_ncrs, iso_gene.start, strand, &spliced);
    *spliced_len = spliced.l;

    free(hap.s);
    free(spliced.s);
    kv_destroy(hap_ncrs);
    return true;
}

// The length of the spliced coding sequence of a record that is in the gene, or -1 when it is not (and is not processed)
static int64_t spliced_length(Record *r, Records *rs, const NonCodingRegions *ncrs, const int strand) {
    if (!variant_in_region(&r->alleles, iso_gene)) return -1;

    NonCodingRegions hap_ncrs = to_haplotype(ncrs, rs);

    memset(dump_reference, 'N', 100);
    dump_reference[100] = '\0';
    const kstring_t reference = { 100, 101, dump_reference };

    kstring_t hap = {0}, spliced = {0};
    get_haplotype_sequence(reference, rs, iso_gene, &hap);
    splice_haplotype_sequence(hap.s != NULL ? hap.s : "", &hap_ncrs, iso_gene.start, strand, &spliced);
    const int64_t length = spliced.l;

    free(hap.s);
    free(spliced.s);
    kv_destroy(hap_ncrs);
    return length;
}

// The length check of a row: the spliced length minus (coding bases of the reference - deleted + inserted), "ok" when 0
static void length_check_cell(char *cell, const size_t size, const int64_t length, const CodingEffect ce) {
    const int64_t expected = 50 - ce.deleted + ce.inserted; // The layout has 50 coding bases
    if (length < 0) snprintf(cell, size, "-");
    else if (length == expected) snprintf(cell, size, "ok");
    else snprintf(cell, size, "%+ld", (long) (length - expected));
}

static void index_cell(char *cell, const size_t size, const bool computed, const int64_t seq_pos, const int64_t index, const int64_t spliced_len) {
    if (computed) snprintf(cell, size, "%ld->%ld/%ld", (long) seq_pos, (long) index, (long) spliced_len);
    else snprintf(cell, size, "-");
}

// The length of the spliced coding sequence follows the coding effect for every variant in the gene, wherever it is, on both
// strands: 50 coding bases, minus the coding bases deleted, plus the bases inserted into the coding sequence. This is the check
// the dump shows in its "len check" column.
// The coding bases of the isoform before any variant: the gene without the used non-coding regions (a retained intron is coding)
static int64_t coding_total(const NonCodingRegions *ncrs) {
    int64_t n = range_len(iso_gene);
    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        if (kv_A(*ncrs, i).is_used) n -= range_len(kv_A(*ncrs, i).coords);
    }
    return n;
}

// Gene layouts other than the default: a zero-length UTR (the CDS reaches the gene edge), two introns, an intron of length 0
// (touching exons), and each of them with the introns spliced out or retained
static const RegionLayout other_layouts[] = {
    { "zero-length UTRs", { 100, 100 }, { 200, 200 }, { { 130, 150 } }, 1 },
    { "two introns", { 100, 110 }, { 190, 200 }, { { 125, 135 }, { 155, 165 } }, 2 },
    { "touching exons", { 100, 110 }, { 180, 200 }, { { 130, 130 } }, 1 },
    { "one UTR only", { 100, 100 }, { 180, 200 }, { { 120, 125 }, { 150, 170 } }, 2 },
    // A UTR that spans exons is one region from the gene edge to the CDS (Utr_range), so an intron inside it is part of the
    // range: these UTRs [100, 140) and [170, 200) can hold introns of their own
    { "UTRs that span exons", { 100, 140 }, { 170, 200 }, { { 150, 160 } }, 1 },
};

// The spliced length follows the coding effect and the position stays near the coding sequence: from -1 (just before it) to
// its length (just after it), for every layout, strand and spliced or retained introns
START_TEST(check_spliced_length__follows_the_coding_effect_in_other_layouts) {
    const int strands[2] = { STRAND_FWD, STRAND_REV };
    const int n_layouts = sizeof(other_layouts) / sizeof(other_layouts[0]);

    for (int li = 0; li < n_layouts; li++) {
        for (int si = 0; si < 2; si++) {
            for (int used = 0; used <= 1; used++) {
                NonCodingRegions ncrs = layout_regions(&other_layouts[li], strands[si], used);
                const int64_t total = coding_total(&ncrs);

                for (int vi = 0; vi < N_ISO_VARIANTS; vi++) {
                    for (int pos = 96; pos <= 203; pos++) {
                        Records rs; kv_init(rs);
                        Record *r = Record_init(mock_bcf1_t(pos, iso_variants[vi][0], iso_variants[vi][1]), NULL);
                        kv_push(Record *, rs, r);

                        char what[128];
                        snprintf(what, sizeof(what), "%s, strand %d, introns %s: %s>%s at %d", other_layouts[li].name, strands[si], used ? "spliced" : "retained", iso_variants[vi][0], iso_variants[vi][1], pos);

                        const int64_t length = spliced_length(r, &rs, &ncrs, strands[si]);
                        if (length >= 0) {
                            const CodingEffect ce = variant_coding_effect(&r->alleles, &ncrs, iso_gene);
                            ck_assert_msg(length == total - ce.deleted + ce.inserted, "%s: length %ld, expected %ld", what, (long) length, (long) (total - ce.deleted + ce.inserted));
                        }

                        int64_t seq_pos = 0, idx = 0, len = 0;
                        if (spliced_index(r, &rs, &ncrs, strands[si], &seq_pos, &idx, &len)) {
                            ck_assert_msg(idx >= -1 && idx <= len, "%s: index %ld of %ld", what, (long) idx, (long) len);
                        }

                        Record_destroy(r); kv_destroy(rs);
                    }
                }

                kv_destroy(ncrs);
            }
        }
    }
} END_TEST

START_TEST(check_spliced_length__follows_the_coding_effect) {
    const int strands[2] = { STRAND_FWD, STRAND_REV };

    for (int si = 0; si < 2; si++) {
        for (int used = 0; used <= 1; used++) {
            NonCodingRegions ncrs = layout_regions(&default_layout, strands[si], used);
            const int64_t total = coding_total(&ncrs);

            for (int vi = 0; vi < N_ISO_VARIANTS; vi++) {
                for (int pos = 96; pos <= 203; pos++) {
                    Records rs; kv_init(rs);
                    Record *r = Record_init(mock_bcf1_t(pos, iso_variants[vi][0], iso_variants[vi][1]), NULL);
                    kv_push(Record *, rs, r);

                    const int64_t length = spliced_length(r, &rs, &ncrs, strands[si]);
                    if (length >= 0) {
                        const CodingEffect ce = variant_coding_effect(&r->alleles, &ncrs, iso_gene);
                        ck_assert_msg(length == total - ce.deleted + ce.inserted, "%s>%s at %d, strand %d, introns %s: length %ld, expected %ld",
                            iso_variants[vi][0], iso_variants[vi][1], pos, strands[si], used ? "spliced" : "retained", (long) length, (long) (total - ce.deleted + ce.inserted));
                    }

                    Record_destroy(r); kv_destroy(rs);
                }
            }

            kv_destroy(ncrs);
        }
    }
} END_TEST

// The position of a variant in the spliced coding sequence is its first changed coding base, in genomic order, the same
// convention as variant_position elsewhere: a substitution or a deletion at its first changed coding base (the first
// deleted base, which after the deletion is where the following base sits), an insertion at the coding base it is
// attached to (the nearest coding base to its left when the attached base is not coding). The index is the index in the
// coding sequence of the haplotype, counted in the reading direction of the strand: forward, 0 is the base 110; reverse,
// 0 is the base 179. Layout: exon 1 [110, 130) is 20 bases, exon 2 [150, 180) is 30 bases.
static void assert_index(const int strand, const int pos, const char *ref, const char *alt, const int64_t expected_index, const int64_t expected_length) {
    NonCodingRegions ncrs = iso_regions(strand, 1);
    Records rs; kv_init(rs);
    Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
    kv_push(Record *, rs, r);

    int64_t seq_pos = 0, index = 0, length = 0;
    const bool computed = spliced_index(r, &rs, &ncrs, strand, &seq_pos, &index, &length);

    ck_assert_msg(computed, "%s>%s at %d, strand %d: not coding", ref, alt, pos, strand);
    ck_assert_msg(index == expected_index && length == expected_length, "%s>%s at %d, strand %d: index %ld of %ld, expected %ld of %ld",
        ref, alt, pos, strand, (long) index, (long) length, (long) expected_index, (long) expected_length);

    kv_destroy(ncrs); Record_destroy(r); kv_destroy(rs);
}

START_TEST(check_spliced_index__forward_strand) {
    assert_index(STRAND_FWD, 115, "N", "x", 5, 50);
    assert_index(STRAND_FWD, 110, "NN", "N", 1, 49);          // Deletes 111: the index of the base after it, 112
    assert_index(STRAND_FWD, 107, "NNNN", "N", 0, 49);        // Deletes 108, 109 (UTR) and 110: the first coding base is the first changed one
    assert_index(STRAND_FWD, 109, "NNNN", "N", 0, 47);        // Deletes 110 to 112, anchored in the UTR
    assert_index(STRAND_FWD, 147, "NNNN", "N", 20, 49);       // Deletes 148, 149 (intron) and 150: the first base of exon 2 is gone
    assert_index(STRAND_FWD, 177, "NNNN", "N", 48, 48);       // Deletes 178, 179, 180: the end of the coding sequence, one past the last base
    assert_index(STRAND_FWD, 115, "N", "Nx", 5, 51);          // Insertion: the base it is attached to
    assert_index(STRAND_FWD, 149, "N", "Nx", 19, 51);         // Between the intron and exon 2: the nearest coding base to its left, 129
    assert_index(STRAND_FWD, 108, "NNN", "Nxy", 0, 50);       // x pairs with 109 (UTR), y with 110
    assert_index(STRAND_FWD, 178, "NNN", "Nx", 49, 50);       // x pairs with 179, 180 (UTR) is deleted
} END_TEST

START_TEST(check_spliced_index__reverse_strand) {
    assert_index(STRAND_REV, 115, "N", "x", 44, 50);
    assert_index(STRAND_REV, 110, "NN", "N", 47, 49);         // Deletes 111: the index of the base on its genomic right, 112, as on the forward strand
    assert_index(STRAND_REV, 107, "NNNN", "N", 48, 49);       // Deletes 108 to 110, the last coding base in reading direction is 110: the base on its right, 111, is the last, 48
    assert_index(STRAND_REV, 115, "N", "Nx", 45, 51);         // Insertion: the base it is attached to, pushed up by the inserted base
    assert_index(STRAND_REV, 149, "N", "Nx", 31, 51);         // Between the intron and exon 2: the nearest coding base to its left, 129
    assert_index(STRAND_REV, 178, "NN", "N", -1, 49);         // Deletes 179, the first coding base in reading direction: its genomic right neighbour 180 is the UTR, one before the coding sequence
    assert_index(STRAND_REV, 178, "NNN", "Nx", 0, 50);        // x pairs with 179, 180 (UTR) is deleted
} END_TEST

/*
* Position truth. Where a base really lands in the spliced coding sequence, against what the position code reports.
*
* The gene [100, 200) is a sequence of distinct bytes (base x is the byte 1 + x - 100, so every base has its own label), a
* single variant is applied with the real get_haplotype_sequence (its inserted bases get labels 102 and up), the regions are
* converted with the real region_to_haplotype and spliced with the real splice_haplotype_sequence. A label is then found in the
* spliced sequence at one place, and its index in the reading direction of the strand (the index counted from the end on the
* reverse strand, which is what reverse_complement does to the sequence) is the truth.
*
* Two functions are compared with it, over every layout, strand, introns spliced or retained, variant kind and position:
*   C  the index of every coding base that is still there, isoform_coding_index (what the splice shifts use)
*   A  the position of the variant itself, get_range_in_isoform (what the indel shifts use). Expected: the first changed coding
*      base in genomic order, where a substituted base is the ALT base that replaces it, a deleted base is the first coding
*      base after the deletion (no such base: the length on the forward strand, -1 on the reverse strand), and an insertion is
*      the coding base it is attached to (the nearest one to its left when that base is not coding)
* Both have to agree with the truth everywhere. A summary and the first differences are written to $POSITION_TRUTH (default
* position_truth.txt). The function the splice shifts used before isoform_coding_index differed from the truth for every
* variant that changes the length of the sequence and lies in an intron or next to the base.
*/
typedef struct TruthVariant { const char *name; int ref_n, alt_n, shared; } TruthVariant;

static const TruthVariant truth_variants[] = {
    { "SNV", 1, 1, 0 }, { "MNV", 2, 2, 0 }, { "ins 1", 1, 2, 1 }, { "ins 3", 1, 4, 1 },
    { "del 1", 2, 1, 1 }, { "del 3", 4, 1, 1 }, { "complex 2>1", 3, 2, 1 }, { "complex 1>3", 2, 4, 1 },
};
#define N_TRUTH_VARIANTS ((int) (sizeof(truth_variants) / sizeof(truth_variants[0])))

static char truth_label(const int64_t x) {
    return x >= iso_gene.start && x < iso_gene.end ? (char) (1 + (x - iso_gene.start)) : (char) 126;
}

// The alleles of a variant of the given shape at pos: REF bases are the labels of the bases there, the ALT starts with the
// shared bases and goes on with new labels
static void truth_alleles(const TruthVariant *tv, const int pos, char *ref, char *alt) {
    for (int i = 0; i < tv->ref_n; i++) ref[i] = truth_label(pos + i);
    ref[tv->ref_n] = '\0';

    for (int i = 0; i < tv->shared; i++) alt[i] = ref[i];
    for (int i = tv->shared; i < tv->alt_n; i++) alt[i] = (char) (101 + i);
    alt[tv->alt_n] = '\0';
}

static bool truth_coding_base(const int64_t b, const NonCodingRegions *ncrs) {
    if (b < iso_gene.start || b >= iso_gene.end) return false;

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        if (ncr.is_used && b >= ncr.coords.start && b < ncr.coords.end) return false;
    }

    return true;
}

typedef struct Truth {
    int64_t where[256]; // The index of every label in the spliced sequence in genomic order, -1 when it is not there
    int64_t length;
    NonCodingRegions hap_ncrs;
} Truth;

static char truth_reference[101];

// Builds the spliced sequence of a haplotype with the real functions. Returns the number of labels that occur twice.
static int truth_build(const Records *rs, const NonCodingRegions *ncrs, const int strand, Truth *t) {
    for (int i = 0; i < 100; i++) truth_reference[i] = truth_label(iso_gene.start + i);
    truth_reference[100] = '\0';
    const kstring_t reference = { 100, 101, truth_reference };

    kstring_t hap = {0}, spliced = {0};
    get_haplotype_sequence(reference, rs, iso_gene, &hap);
    t->hap_ncrs = to_haplotype(ncrs, rs);
    splice_haplotype_sequence(hap.s != NULL ? hap.s : "", &t->hap_ncrs, iso_gene.start, strand, &spliced);

    int duplicates = 0;
    for (int i = 0; i < 256; i++) t->where[i] = -1;
    for (size_t i = 0; i < spliced.l; i++) {
        const unsigned char label = (unsigned char) spliced.s[i];
        if (t->where[label] >= 0) duplicates++;
        t->where[label] = (int64_t) i;
    }
    t->length = (int64_t) spliced.l;

    free(hap.s);
    free(spliced.s);
    return duplicates;
}

static bool truth_index(const Truth *t, const int strand, const unsigned char label, int64_t *index) {
    if (t->where[label] < 0) return false;

    *index = strand == STRAND_FWD ? t->where[label] : t->length - 1 - t->where[label];
    return true;
}

// Where the variant should be placed, by the first changed base convention. False when there is nothing to compare.
static bool truth_expected_variant_index(const Alleles *a, const NonCodingRegions *ncrs, const Truth *t, const int strand, int64_t *expected) {
    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
    unsigned char label;

    int64_t first = -1;
    for (int64_t b = a->pos; b < a->pos + a->ref_len; b++) {
        if (truth_coding_base(b, ncrs)) { first = b; break; }
    }

    if (first >= 0 && first - a->pos < paired) {
        label = (unsigned char) a->alt[first - a->pos]; // The ALT base that replaces it
    } else if (first >= 0) {
        // Deleted: the first coding base after the deletion takes its place
        int64_t after = -1;
        for (int64_t b = a->pos + a->ref_len; b < iso_gene.end; b++) {
            if (truth_coding_base(b, ncrs)) { after = b; break; }
        }
        if (after < 0) {
            *expected = strand == STRAND_FWD ? t->length : -1;
            return true;
        }
        label = (unsigned char) truth_label(after);
    } else {
        // Only inserted bases are coding: the nearest coding base to the left of the insertion point
        int64_t left = -1;
        for (int64_t b = a->pos + paired - 1; b >= iso_gene.start; b--) {
            if (truth_coding_base(b, ncrs)) { left = b; break; }
        }
        if (left < 0) return false;
        label = (unsigned char) truth_label(left);
    }

    return truth_index(t, strand, label, expected);
}

START_TEST(check_position_truth__positions_match_the_spliced_sequence) {
    const char *path = getenv("POSITION_TRUTH");
    if (path == NULL) path = "position_truth.txt";

    const RegionLayout *layouts[1 + sizeof(other_layouts) / sizeof(other_layouts[0])];
    int n_layouts = 0;
    layouts[n_layouts++] = &default_layout;
    for (size_t i = 0; i < sizeof(other_layouts) / sizeof(other_layouts[0]); i++) layouts[n_layouts++] = &other_layouts[i];

    const int strands[2] = { STRAND_FWD, STRAND_REV };
    long a_total = 0, a_bad = 0, c_total = 0, c_bad = 0, duplicates_total = 0;
    long a_by_kind[N_TRUTH_VARIANTS][2] = {{0}};
    int a_examples = 0;
    int a_shown[N_TRUTH_VARIANTS][8][2][2] = {{{{0}}}}; // Examples per bucket
    kstring_t a_lines = {0}, c_lines = {0};

    for (int li = 0; li < n_layouts; li++) {
        for (int si = 0; si < 2; si++) {
            for (int used = 0; used <= 1; used++) {
                NonCodingRegions ncrs = layout_regions(layouts[li], strands[si], used);

                for (int vi = 0; vi < N_TRUTH_VARIANTS; vi++) {
                    for (int pos = 96; pos <= 203; pos++) {
                        char ref[16], alt[16];
                        truth_alleles(&truth_variants[vi], pos, ref, alt);

                        Records rs; kv_init(rs);
                        Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
                        kv_push(Record *, rs, r);

                        if (variant_in_region(&r->alleles, iso_gene)) {
                            Truth t;
                            duplicates_total += truth_build(&rs, &ncrs, strands[si], &t);

                            // C: every coding base that is still there, by isoform_coding_index
                            for (int64_t x = iso_gene.start; x < iso_gene.end; x++) {
                                int64_t want = 0;
                                if (!truth_coding_base(x, &ncrs) || !truth_index(&t, strands[si], (unsigned char) truth_label(x), &want)) continue;

                                const int64_t by_map = isoform_coding_index(&rs, &ncrs, iso_gene, strands[si], x);
                                c_total++;
                                if (by_map != want) {
                                    c_bad++;
                                    if (c_bad <= 20) ksprintf(&c_lines, "%-20s strand %d %-8s %-11s at %3d: base %3ld: index %ld, truth %ld\n", layouts[li]->name, strands[si], used ? "spliced" : "retained", truth_variants[vi].name, pos, (long) x, (long) by_map, (long) want);
                                }
                            }

                            // A: the variant itself
                            int64_t seq_pos = 0, expected = 0;
                            int offset = 0;
                            if (get_range_in_isoform(&r->alleles, &ncrs, &rs, strands[si], iso_gene, &seq_pos, &offset) == 0
                                && truth_expected_variant_index(&r->alleles, &ncrs, &t, strands[si], &expected)) {
                                const int64_t got = seq_pos; // The position is the index into the coding sequence

                                a_total++;
                                if (got != expected) {
                                    a_bad++; a_by_kind[vi][si]++;
                                    int *shown = &a_shown[vi][li][si][used];
                                    if (*shown < 2) {
                                        (*shown)++; a_examples++;
                                        ksprintf(&a_lines, "%-20s strand %d %-8s %-11s at %3d: index %ld of %ld, expected %ld\n", layouts[li]->name, strands[si], used ? "spliced" : "retained", truth_variants[vi].name, pos, (long) got, (long) t.length, (long) expected);
                                    }
                                }
                            }

                            kv_destroy(t.hap_ncrs);
                        }

                        Record_destroy(r);
                        kv_destroy(rs);
                    }
                }

                kv_destroy(ncrs);
            }
        }
    }

    FILE *out = fopen(path, "w");
    ck_assert_msg(out != NULL, "cannot write %s", path);
    fprintf(out, "A (the position of the variant, get_range_in_isoform): %ld comparisons, %ld differ\n", a_total, a_bad);
    fprintf(out, "C (every coding base, isoform_coding_index): %ld comparisons, %ld differ\n", c_total, c_bad);
    fprintf(out, "labels that occur twice in a spliced sequence (the harness assumes none): %ld\n\n", duplicates_total);
    fprintf(out, "A differences by variant kind (forward, reverse)\n");
    for (int vi = 0; vi < N_TRUTH_VARIANTS; vi++) fprintf(out, "%-12s %8ld %8ld\n", truth_variants[vi].name, a_by_kind[vi][0], a_by_kind[vi][1]);
    fprintf(out, "\nA, %d examples (two for every kind, layout, strand and spliced or retained introns that differs)\n%s", a_examples, a_lines.s != NULL ? a_lines.s : "");
    fprintf(out, "\nC, the first differences\n%s", c_lines.s != NULL ? c_lines.s : "");
    fclose(out);
    free(a_lines.s);
    free(c_lines.s);

    ck_assert_int_gt(a_total, 0);
    ck_assert_int_gt(c_total, 0);
    ck_assert_int_eq(duplicates_total, 0);
    ck_assert_msg(a_bad == 0, "get_range_in_isoform differs from the truth in %ld of %ld cases (see %s)", a_bad, a_total, path);
    ck_assert_msg(c_bad == 0, "isoform_coding_index differs from the truth in %ld of %ld cases (see %s)", c_bad, c_total, path);
} END_TEST

static int in_range(const int64_t x, const Range r) {
    return x >= r.start && x < r.end;
}

static int ref_overlaps(const Alleles *a, const Range r) {
    return a->ref_len > 0 && a->pos < r.end && a->pos + a->ref_len > r.start;
}

START_TEST(check_boundary_table__dump) {
    const char *path = getenv("BOUNDARY_TABLE");
    if (path == NULL) path = "boundary_table.txt";
    FILE *out = fopen(path, "w");
    ck_assert_msg(out != NULL, "cannot write %s", path);

    const struct { const char *label; Range region; int64_t edge; } borders[] = {
        { "gene start", { 100, 200 }, 100 }, { "UTR end", { 100, 110 }, 110 },
        { "intron start", { 130, 150 }, 130 }, { "intron end", { 130, 150 }, 150 },
        { "UTR start", { 180, 200 }, 180 }, { "gene end", { 100, 200 }, 200 },
    };
    const char *kinds[N_ISO_VARIANTS] = { "SNV", "MNV", "ins 1", "ins 3", "del 1", "del 3", "complex 2>1", "complex 1>3" };

    NonCodingRegions fwd = iso_regions(STRAND_FWD, 1);
    NonCodingRegions rev = iso_regions(STRAND_REV, 1);

    fprintf(out, "kind         ref>alt        VCF pos  border          region      A B E C D | fwd: nc  p  off | rev: nc  p  off | utr f r | fwd seq_pos->idx/len | rev seq_pos->idx/len | ce: c s d i | len check f r\n");
    int differing_membership = 0, differing_classification = 0, rows = 0;

    for (size_t bi = 0; bi < sizeof(borders) / sizeof(borders[0]); bi++) {
        for (int vi = 0; vi < N_ISO_VARIANTS; vi++) {
            for (int pos = borders[bi].edge - 4; pos <= borders[bi].edge + 3; pos++) {
                Records rs; kv_init(rs);
                Record *r = Record_init(mock_bcf1_t(pos, iso_variants[vi][0], iso_variants[vi][1]), NULL);
                kv_push(Record *, rs, r);
                const Alleles *a = &r->alleles;
                const Range region = borders[bi].region;

                const int A = variant_overlaps(a, region);
                const int B = ref_overlaps(a, region);
                const int is_insertion = a->ref_len == 0;
                const int E = is_insertion && in_range(a->pos - 1, region) && in_range(a->pos, region);
                const int C = in_range(variant_footprint(a).start, region);
                const int D = in_range(variant_position(a), region);

                int64_t fwd_pos = -1, rev_pos = -1;
                int fwd_off = 0, rev_off = 0;
                const int fwd_nc = get_range_in_isoform(a, &fwd, &rs, STRAND_FWD, iso_gene, &fwd_pos, &fwd_off);
                const int rev_nc = get_range_in_isoform(a, &rev, &rs, STRAND_REV, iso_gene, &rev_pos, &rev_off);
                const int utr_f = is_utr_variant(a, &fwd);
                const int utr_r = is_utr_variant(a, &rev);

                int64_t fwd_seq_pos = 0, fwd_index = 0, fwd_len = 0, rev_seq_pos = 0, rev_index = 0, rev_len = 0;
                const bool fwd_has = spliced_index(r, &rs, &fwd, STRAND_FWD, &fwd_seq_pos, &fwd_index, &fwd_len);
                const bool rev_has = spliced_index(r, &rs, &rev, STRAND_REV, &rev_seq_pos, &rev_index, &rev_len);
                char fwd_cell[48], rev_cell[48];
                index_cell(fwd_cell, sizeof(fwd_cell), fwd_has, fwd_seq_pos, fwd_index, fwd_len);
                index_cell(rev_cell, sizeof(rev_cell), rev_has, rev_seq_pos, rev_index, rev_len);

                const CodingEffect ce = variant_coding_effect(a, &fwd, iso_gene);

                char fwd_check[16], rev_check[16];
                length_check_cell(fwd_check, sizeof(fwd_check), spliced_length(r, &rs, &fwd, STRAND_FWD), ce);
                length_check_cell(rev_check, sizeof(rev_check), spliced_length(r, &rs, &rev, STRAND_REV), ce);

                const int membership_differs = !is_insertion && A != B;
                const int classification_differs = C != D;
                differing_membership += membership_differs;
                differing_classification += classification_differs;
                rows++;

                fprintf(out, "%-12s %-5s>%-5s %7d  %-14s [%3ld,%3ld)  %d %c %c %d %d%s%s | %d %4ld %3d | %d %4ld %3d | %d %d | %-18s | %-18s | %d %d %d %d | %-3s %-3s\n",
                    kinds[vi], iso_variants[vi][0], iso_variants[vi][1], pos, borders[bi].label, (long) region.start, (long) region.end,
                    A, is_insertion ? '-' : '0' + B, is_insertion ? '0' + E : '-', C, D,
                    membership_differs ? " !" : "", classification_differs ? " *" : "",
                    fwd_nc, (long) fwd_pos, fwd_nc ? 0 : fwd_off, rev_nc, (long) rev_pos, rev_nc ? 0 : rev_off, utr_f, utr_r,
                    fwd_cell, rev_cell, ce.coding, ce.substituted, ce.deleted, ce.inserted, fwd_check, rev_check);

                Record_destroy(r);
                kv_destroy(rs);
            }
        }
    }

    fprintf(out, "\n%d rows, %d where footprint membership and changed-base membership differ (!), %d where the attachment base and variant_position classify differently (*)\n",
        rows, differing_membership, differing_classification);
    fclose(out);
    kv_destroy(fwd); kv_destroy(rev);

    ck_assert_int_gt(rows, 0);
} END_TEST

/*
* unique_isoform_shifts: the shifts that differ in type, position or modified bases, in order. Two shifts that agree on those
* three are the same shift, whichever record they come from or what they shift by: the first one is kept.
*/

static IsoformVariant shift_at(const enum ShiftType type, const int64_t seq_pos, const int shift, const int modified_bases, Record *r) {
    return (IsoformVariant) { 1, 0, type, seq_pos, shift, modified_bases, r };
}

START_TEST(check_unique_isoform_shifts__empty) {
    IsoformVariants shifts; kv_init(shifts);
    IsoformVariants unique; kv_init(unique);

    unique_isoform_shifts(&shifts, &unique);

    ck_assert_int_eq(kv_size(unique), 0);

    kv_destroy(shifts); kv_destroy(unique);
} END_TEST

START_TEST(check_unique_isoform_shifts__different_shifts_are_all_kept_in_order) {
    Record *r = Record_init(mock_bcf1_t(10, "N", "x"), NULL);
    IsoformVariants shifts; kv_init(shifts);
    kv_push(IsoformVariant, shifts, shift_at(INDEL, 5, 1, 1, r));
    kv_push(IsoformVariant, shifts, shift_at(SPLICE, 5, 1, 1, r));  // Other type
    kv_push(IsoformVariant, shifts, shift_at(INDEL, 6, 1, 1, r));   // Other position
    kv_push(IsoformVariant, shifts, shift_at(INDEL, 5, 1, 2, r));   // Other modified bases
    IsoformVariants unique; kv_init(unique);

    unique_isoform_shifts(&shifts, &unique);

    ck_assert_int_eq(kv_size(unique), 4);
    for (size_t i = 0; i < 4; i++) {
        ck_assert_int_eq(kv_A(unique, i).seq_pos, kv_A(shifts, i).seq_pos);
        ck_assert_int_eq(kv_A(unique, i).shift_type, kv_A(shifts, i).shift_type);
        ck_assert_int_eq(kv_A(unique, i).modified_bases, kv_A(shifts, i).modified_bases);
    }

    Record_destroy(r);
    kv_destroy(shifts); kv_destroy(unique);
} END_TEST

START_TEST(check_unique_isoform_shifts__a_repeated_shift_is_kept_once_at_its_first_place) {
    Record *r = Record_init(mock_bcf1_t(10, "N", "x"), NULL);
    IsoformVariants shifts; kv_init(shifts);
    kv_push(IsoformVariant, shifts, shift_at(SPLICE, 5, -3, -3, r));
    kv_push(IsoformVariant, shifts, shift_at(INDEL, 8, 1, 1, r));
    kv_push(IsoformVariant, shifts, shift_at(SPLICE, 5, -3, -3, r)); // Repeat
    kv_push(IsoformVariant, shifts, shift_at(INDEL, 20, -1, -1, r));
    IsoformVariants unique; kv_init(unique);

    unique_isoform_shifts(&shifts, &unique);

    ck_assert_int_eq(kv_size(unique), 3);
    ck_assert_int_eq(kv_A(unique, 0).seq_pos, 5);
    ck_assert_int_eq(kv_A(unique, 1).seq_pos, 8);
    ck_assert_int_eq(kv_A(unique, 2).seq_pos, 20); // The last shift is still there, which the unfixed loop dropped

    Record_destroy(r);
    kv_destroy(shifts); kv_destroy(unique);
} END_TEST

// The key leaves out the record: of two records that changed the same site the first one is kept. This pins what the key does,
// it is not a statement that it is right (the focal record can be the one dropped).
START_TEST(check_unique_isoform_shifts__the_record_is_not_part_of_the_key) {
    Record *first = Record_init(mock_bcf1_t(10, "N", "x"), NULL);
    Record *second = Record_init(mock_bcf1_t(11, "N", "x"), NULL);
    IsoformVariants shifts; kv_init(shifts);
    kv_push(IsoformVariant, shifts, shift_at(SPLICE, 5, -3, -3, first));
    kv_push(IsoformVariant, shifts, shift_at(SPLICE, 5, -3, -3, second));
    IsoformVariants unique; kv_init(unique);

    unique_isoform_shifts(&shifts, &unique);

    ck_assert_int_eq(kv_size(unique), 1);
    ck_assert_ptr_eq(kv_A(unique, 0).r, first);

    Record_destroy(first); Record_destroy(second);
    kv_destroy(shifts); kv_destroy(unique);
} END_TEST

Suite *isoform_position_suite() {
    Suite *s = suite_create("Isoform position");

    TCase *tc = mocked_tcase("isoform_position");
    suite_add_tcase(s, tc);
    tcase_set_timeout(tc, 60);

    tcase_add_test(tc, check_isoform_position__classification_follows_the_coding_effect);
    tcase_add_test(tc, check_isoform_coding_index__without_variants);
    tcase_add_test(tc, check_isoform_coding_index__a_variant_in_an_intron_changes_nothing);
    tcase_add_test(tc, check_isoform_coding_index__a_coding_variant_moves_what_comes_after_it);
    tcase_add_test(tc, check_get_range_in_isoform__an_insertion_at_the_edge_of_a_utr_is_noncoding);
    tcase_add_test(tc, check_get_range_in_isoform__only_pure_insertions_change_at_the_edge_of_a_utr);
    tcase_add_test(tc, check_isoform_base_is_coding__follows_the_used_regions);
    tcase_add_test(tc, check_is_utr_variant__only_the_reference_position_counts);
    tcase_add_test(tc, check_is_utr_variant__substitutions_need_a_base_in_the_utr);
    tcase_add_test(tc, check_is_utr_variant__indels_count_by_the_changed_bases);
    tcase_add_test(tc, check_is_utr_variant__an_insertion_is_in_the_utr_at_either_edge);
    tcase_add_test(tc, check_variant_in_cds_exon__substitutions_and_deletions_count_at_their_first_changed_base);
    tcase_add_test(tc, check_variant_in_cds_exon__an_insertion_inside_the_cds_is_in);
    tcase_add_test(tc, check_variant_in_cds_exon__an_insertion_at_the_border_with_a_utr_is_not_in);
    tcase_add_test(tc, check_variant_in_cds_exon__an_insertion_at_the_border_with_an_intron_is_in);
    tcase_add_test(tc, check_variant_coding_effect__substitutions);
    tcase_add_test(tc, check_variant_coding_effect__a_deletion_applies_only_its_coding_part);
    tcase_add_test(tc, check_variant_coding_effect__complex_indels_pair_from_the_left);
    tcase_add_test(tc, check_variant_coding_effect__insertions_follow_the_border_rule);
    tcase_add_test(tc, check_spliced_length__follows_the_coding_effect);
    tcase_add_test(tc, check_spliced_length__follows_the_coding_effect_in_other_layouts);
    tcase_add_test(tc, check_spliced_index__forward_strand);
    tcase_add_test(tc, check_spliced_index__reverse_strand);
    tcase_add_test(tc, check_position_truth__positions_match_the_spliced_sequence);
    tcase_add_test(tc, check_boundary_table__dump);
    tcase_add_test(tc, check_unique_isoform_shifts__empty);
    tcase_add_test(tc, check_unique_isoform_shifts__different_shifts_are_all_kept_in_order);
    tcase_add_test(tc, check_unique_isoform_shifts__a_repeated_shift_is_kept_once_at_its_first_place);
    tcase_add_test(tc, check_unique_isoform_shifts__the_record_is_not_part_of_the_key);

    return s;
}
