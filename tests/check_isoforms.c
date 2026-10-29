#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <cpliceai/range.h>

#include "bcftools/gff.h"
#include "haplotype.h"
#include "isoform_position.h"
#include "isoforms.h"
#include "klib/kvec.h"
#include "mocks.h"
#include "record.h"
#include "session.h"
#include "splice_prediction.h"
#include "triggers.h"

#include "check_isoforms.h"

/*
* calculate_rescue_probabilities and calculate_ptv_probabilities: the totals of an Isoforms struct from its isoforms, whose
* posterior probabilities are already set.
*/

typedef struct IsoformSpec {
    uint8_t rescue;
    enum RescueType rescue_type;
    uint32_t affected_bases;
    float posterior;
    uint8_t ptv;
    enum PTVType ptv_type;
} IsoformSpec;

// The state in which build_isoforms_for_phase hands a struct to the calculations. Freed with Isoforms_destroy.
static Isoforms isoforms_of(const IsoformSpec *specs, const int n) {
    Isoforms isoforms = { .valid = 1, .m_rescue_types = 1, .ptv_type = NO_PTV };
    kv_init(isoforms.isoforms);

    for (int i = 0; i < n; i++) {
        Isoform isoform = {
            .posterior_probability = specs[i].posterior,
            .rescue = specs[i].rescue,
            .rescue_type = specs[i].rescue_type,
            .affected_bases = specs[i].affected_bases,
            .ptv = specs[i].ptv,
            .ptv_type = specs[i].ptv_type
        };
        kv_push(Isoform, isoforms.isoforms, isoform);
    }

    return isoforms;
}

START_TEST(check_calculate_rescue_probabilities__two_types) {
    const IsoformSpec specs[] = {
        { 1, EXON_SKIPPING, 3, 0.6f, 0, NO_PTV },
        { 0, NO_RESCUE, 0, 0.4f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_rescue_probabilities(&isoforms);

    ck_assert_float_eq_tol(isoforms.rescue_probability, 0.6f, 1e-6);
    ck_assert_int_eq(isoforms.rescue, 1);
    ck_assert_int_eq(isoforms.n_rescue_types, 2);

    const RescueProbabilities *probs = isoforms.rescue_type_probabilities;
    ck_assert_int_eq(probs->rescue_types[0], EXON_SKIPPING);
    ck_assert_int_eq(probs->affected_bases[0], 3);
    ck_assert_int_eq(probs->rescues[0], 1);
    ck_assert_int_eq(probs->n_isoforms[0], 1);
    ck_assert_float_eq_tol(probs->rescue_probabilities[0], 0.6f, 1e-6);

    ck_assert_int_eq(probs->rescue_types[1], NO_RESCUE);
    ck_assert_int_eq(probs->affected_bases[1], 0);
    ck_assert_int_eq(probs->rescues[1], 0);
    ck_assert_int_eq(probs->n_isoforms[1], 1);
    ck_assert_float_eq_tol(probs->rescue_probabilities[1], 0.4f, 1e-6);

    Isoforms_destroy(&isoforms);
} END_TEST

// Isoforms with the same rescue type and the same affected bases are one entry; another number of affected bases is another one
START_TEST(check_calculate_rescue_probabilities__same_type_and_bases_are_added_up) {
    const IsoformSpec specs[] = {
        { 1, INTRON_RETENTION, 21, 0.3f, 0, NO_PTV },
        { 1, INTRON_RETENTION, 21, 0.2f, 0, NO_PTV },
        { 1, INTRON_RETENTION, 24, 0.5f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 3);

    calculate_rescue_probabilities(&isoforms);

    ck_assert_int_eq(isoforms.n_rescue_types, 2);
    const RescueProbabilities *probs = isoforms.rescue_type_probabilities;
    ck_assert_int_eq(probs->affected_bases[0], 21);
    ck_assert_int_eq(probs->n_isoforms[0], 2);
    ck_assert_float_eq_tol(probs->rescue_probabilities[0], 0.5f, 1e-6);
    ck_assert_int_eq(probs->affected_bases[1], 24);
    ck_assert_int_eq(probs->n_isoforms[1], 1);
    ck_assert_float_eq_tol(probs->rescue_probabilities[1], 0.5f, 1e-6);
    ck_assert_float_eq_tol(isoforms.rescue_probability, 1.0f, 1e-6);

    Isoforms_destroy(&isoforms);
} END_TEST

START_TEST(check_calculate_rescue_probabilities__many_types_grow_the_arrays) {
    const IsoformSpec specs[] = {
        { 1, COMPLETE_RESCUE, 0, 0.1f, 0, NO_PTV },
        { 1, EXON_SKIPPING, 6, 0.2f, 0, NO_PTV },
        { 0, INTRON_RETENTION, 9, 0.3f, 0, NO_PTV },
        { 0, FRAME_RESTORING_INDEL, 15, 0.4f, 0, NO_PTV },
        { 0, PARTIAL_EXON_SKIPPING, 12, 0.0f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 5);

    calculate_rescue_probabilities(&isoforms);

    ck_assert_int_eq(isoforms.n_rescue_types, 5);
    ck_assert_int_ge(isoforms.m_rescue_types, 5);
    const RescueProbabilities *probs = isoforms.rescue_type_probabilities;
    for (int i = 0; i < 5; i++) {
        ck_assert_int_eq(probs->rescue_types[i], specs[i].rescue_type);
        ck_assert_int_eq(probs->affected_bases[i], specs[i].affected_bases);
        ck_assert_float_eq_tol(probs->rescue_probabilities[i], specs[i].posterior, 1e-6);
        ck_assert_int_eq(probs->rescues[i], specs[i].rescue);
    }
    ck_assert_float_eq_tol(isoforms.rescue_probability, 0.3f, 1e-6); // 0.1 + 0.2
    ck_assert_int_eq(isoforms.rescue, 0);

    Isoforms_destroy(&isoforms);
} END_TEST

// The rescue call needs more than half of the probability
START_TEST(check_calculate_rescue_probabilities__half_is_not_a_rescue) {
    const IsoformSpec specs[] = {
        { 1, COMPLETE_RESCUE, 0, 0.5f, 0, NO_PTV },
        { 0, NO_RESCUE, 0, 0.5f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_rescue_probabilities(&isoforms);

    ck_assert_float_eq(isoforms.rescue_probability, 0.5f);
    ck_assert_int_eq(isoforms.rescue, 0);

    Isoforms_destroy(&isoforms);
} END_TEST

START_TEST(check_calculate_ptv_probabilities__the_largest_type_wins) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.6f, 1, STOP_GAINED },
        { 0, NO_RESCUE, 0, 0.1f, 1, FRAMESHIFT },
        { 0, NO_RESCUE, 0, 0.3f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 3);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_float_eq_tol(isoforms.ptv_probability, 0.7f, 1e-6);
    ck_assert_int_eq(isoforms.ptv, 1);
    ck_assert_int_eq(isoforms.ptv_type, STOP_GAINED);

    kv_destroy(isoforms.isoforms);
} END_TEST

START_TEST(check_calculate_ptv_probabilities__frameshift_can_win) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.5f, 1, FRAMESHIFT },
        { 0, NO_RESCUE, 0, 0.2f, 1, STOP_GAINED },
        { 0, NO_RESCUE, 0, 0.3f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 3);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_int_eq(isoforms.ptv, 1);
    ck_assert_int_eq(isoforms.ptv_type, FRAMESHIFT);

    kv_destroy(isoforms.isoforms);
} END_TEST

// Equal probability for stop gained and frameshift gives no type
START_TEST(check_calculate_ptv_probabilities__a_tie_gives_no_type) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.5f, 1, STOP_GAINED },
        { 0, NO_RESCUE, 0, 0.5f, 1, FRAMESHIFT },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_int_eq(isoforms.ptv, 1);
    ck_assert_int_eq(isoforms.ptv_type, NO_PTV);

    kv_destroy(isoforms.isoforms);
} END_TEST

START_TEST(check_calculate_ptv_probabilities__a_long_inframe_indel_alone) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.9f, 1, LONG_INFRAME_INDEL },
        { 0, NO_RESCUE, 0, 0.1f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_float_eq_tol(isoforms.ptv_probability, 0.9f, 1e-6);
    ck_assert_int_eq(isoforms.ptv, 1);
    ck_assert_int_eq(isoforms.ptv_type, LONG_INFRAME_INDEL);

    kv_destroy(isoforms.isoforms);
} END_TEST

// An isoform with ptv 0 does not count, whatever its type says
START_TEST(check_calculate_ptv_probabilities__isoforms_without_a_ptv_are_ignored) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.8f, 0, STOP_GAINED },
        { 0, NO_RESCUE, 0, 0.2f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_float_eq(isoforms.ptv_probability, 0.0f);
    ck_assert_int_eq(isoforms.ptv, 0);
    ck_assert_int_eq(isoforms.ptv_type, NO_PTV);

    kv_destroy(isoforms.isoforms);
} END_TEST

// A ptv needs more than half of the probability
START_TEST(check_calculate_ptv_probabilities__half_is_not_a_ptv) {
    const IsoformSpec specs[] = {
        { 0, NO_RESCUE, 0, 0.5f, 1, FRAMESHIFT },
        { 0, NO_RESCUE, 0, 0.5f, 0, NO_PTV },
    };
    Isoforms isoforms = isoforms_of(specs, 2);

    calculate_ptv_probabilities(&isoforms);

    ck_assert_int_eq(isoforms.ptv, 0);
    ck_assert_int_eq(isoforms.ptv_type, FRAMESHIFT); // The type still follows the largest probability

    kv_destroy(isoforms.isoforms);
} END_TEST

/*
* shift_rescue_type
*/

START_TEST(check_shift_rescue_type__sign_and_partial) {
    IsoformVariant shift = { .shift = -21, .partial = 0 };
    ck_assert_int_eq(shift_rescue_type(&shift), INTRON_RETENTION);
    shift.partial = 1;
    ck_assert_int_eq(shift_rescue_type(&shift), PARTIAL_INTRON_RETENTION);

    shift.shift = 30;
    ck_assert_int_eq(shift_rescue_type(&shift), PARTIAL_EXON_SKIPPING);
    shift.partial = 0;
    ck_assert_int_eq(shift_rescue_type(&shift), EXON_SKIPPING);
} END_TEST

/*
* Fixtures for a gene [100, 200) on the forward strand with two coding exons: [100, 130) with CDS [110, 130), and [151, 200)
* with CDS [151, 190). The 21 bases [130, 151) are the intron, the 5' UTR is [100, 110) and the 3' UTR [190, 200). The bases of
* the gene follow the pattern ACG, which holds no stop codon, so a stop is only there where a test puts one.
*/

static Gene *new_gene(const int64_t start, const int64_t end, const uint32_t strand) {
    Gene *gene = calloc(1, sizeof(Gene));
    gene->gff_id = 1;
    gene->name = "GENE";
    gene->coords = (Range) { start, end };
    gene->strand = strand;

    const int64_t length = end - start;
    gene->seq.s = malloc(length + 1);
    for (int64_t i = 0; i < length; i++) gene->seq.s[i] = "ACG"[i % 3];
    gene->seq.s[length] = '\0';
    gene->seq.l = length;
    gene->seq.m = length + 1;

    return gene;
}

static void set_bases(Gene *gene, const int64_t pos, const char *bases) {
    memcpy(gene->seq.s + (pos - gene->coords.start), bases, strlen(bases));
}

static void free_gene(Gene *gene) {
    free(gene->seq.s);
    free(gene);
}

static Transcript *two_exon_transcript(const uint32_t strand) {
    Range cds[] = { {110, 129}, {151, 189} };
    int constitutive[] = { 1, 1 };
    Range exons[] = { {100, 129}, {151, 199} };
    TranscriptData data = { .id = 1, .start = 100, .end = 199, .strand = strand,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 2 };

    Transcript *t = Transcript_build(&data);
    for (size_t ei = 0; ei < kv_size(t->exons); ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    return t;
}

// One coding exon [100, 160) with CDS [110, 149): the 5' UTR is [100, 110) and the 3' UTR [149, 160)
static Transcript *one_exon_transcript(void) {
    Range cds[] = { {110, 148} };
    int constitutive[] = { 1 };
    Range exons[] = { {100, 159} };
    TranscriptData data = { .id = 1, .start = 100, .end = 159, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 1, .exons = exons, .nexons = 1 };

    return Transcript_build(&data);
}

static void destroy_transcript_with_sites(Transcript *t) {
    for (size_t i = 0; i < kv_size(t->exons); i++) {
        Exon *e = kv_A(t->exons, i);
        if (e->splice_sites != NULL) SpliceSites_destroy(e->splice_sites);
        e->splice_sites = NULL;
    }
    Transcript_destroy(t);
}

// A record with one flag slot per phase for its single transcript
static Record *new_record(const int64_t pos, const char *ref, const char *alt) {
    Record *r = Record_init(mock_bcf1_t(pos, ref, alt), NULL);
    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));

    return r;
}

/*
* build_noncoding_regions
*/

START_TEST(check_build_noncoding_regions__a_transcript_without_predictions) {
    Transcript *t = two_exon_transcript(STRAND_FWD);
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Record *r = new_record(120, "N", "x");

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));

    // 5' UTR, the intron and 3' UTR
    ck_assert_int_eq(kv_size(groups), 3);

    const NonCodingRegions *utr5 = &kv_A(groups, 0);
    ck_assert_int_eq(kv_size(*utr5), 1);
    ck_assert_int_eq(kv_A(*utr5, 0).is_used, 1);
    ck_assert_int_eq(kv_A(*utr5, 0).is_utr, UTR5);
    ck_assert_int_eq(kv_A(*utr5, 0).coords.start, 100);
    ck_assert_int_eq(kv_A(*utr5, 0).coords.end, 110);
    ck_assert_int_eq(kv_size(kv_A(*utr5, 0).records), 1);
    ck_assert_ptr_eq(kv_A(kv_A(*utr5, 0).records, 0), r);

    const NonCodingRegions *introns = &kv_A(groups, 1);
    ck_assert_int_eq(kv_size(*introns), 1);
    const NonCodingRegion intron = kv_A(*introns, 0);
    ck_assert_int_eq(intron.is_used, 1);
    ck_assert_int_eq(intron.is_utr, 0);
    ck_assert_int_eq(intron.coords.start, 130);
    ck_assert_int_eq(intron.coords.end, 151);
    ck_assert_int_eq(intron.canonical_coords.start, 130);
    ck_assert_int_eq(intron.canonical_coords.end, 151);
    ck_assert_int_eq(intron.uses_competitor, 0);
    ck_assert_ptr_eq(intron.donor, &kv_A(t->exons, 0)->splice_sites->canonical_donors[0]);
    ck_assert_ptr_eq(intron.acceptor, &kv_A(t->exons, 1)->splice_sites->canonical_acceptors[0]);

    const NonCodingRegions *utr3 = &kv_A(groups, 2);
    ck_assert_int_eq(kv_size(*utr3), 1);
    ck_assert_int_eq(kv_A(*utr3, 0).is_utr, UTR3);
    ck_assert_int_eq(kv_A(*utr3, 0).coords.start, 190);
    ck_assert_int_eq(kv_A(*utr3, 0).coords.end, 200);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    destroy_transcript_with_sites(t);
} END_TEST

// Both canonical sites of the intron are lost: the intron is retained, which makes it a region that is not used
START_TEST(check_build_noncoding_regions__an_intron_with_both_sites_lost_is_retained) {
    Transcript *t = two_exon_transcript(STRAND_FWD);
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Record *r = new_record(129, "N", "x");

    SpliceSite *donor = &kv_A(t->exons, 0)->splice_sites->canonical_donors[0];
    SpliceSite *acceptor = &kv_A(t->exons, 1)->splice_sites->canonical_acceptors[0];
    donor->predicted = 1; donor->posterior_probability = 0;
    acceptor->predicted = 1; acceptor->posterior_probability = 0;

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));
    ck_assert_int_eq(kv_size(groups), 3);

    const NonCodingRegions *introns = &kv_A(groups, 1);
    ck_assert_int_eq(kv_size(*introns), 1);
    ck_assert_int_eq(kv_A(*introns, 0).is_used, 0);
    ck_assert_int_eq(kv_A(*introns, 0).coords.start, 130);
    ck_assert_int_eq(kv_A(*introns, 0).coords.end, 151);
    ck_assert_ptr_eq(kv_A(*introns, 0).donor, donor);
    ck_assert_ptr_eq(kv_A(*introns, 0).acceptor, acceptor);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    destroy_transcript_with_sites(t);
} END_TEST

// Only the donor is lost: the intron is retained once for every used acceptor
START_TEST(check_build_noncoding_regions__a_lost_donor_retains_the_intron) {
    Transcript *t = two_exon_transcript(STRAND_FWD);
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Record *r = new_record(129, "N", "x");

    SpliceSite *donor = &kv_A(t->exons, 0)->splice_sites->canonical_donors[0];
    donor->predicted = 1; donor->posterior_probability = 0;

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));

    const NonCodingRegions *introns = &kv_A(groups, 1);
    ck_assert_int_eq(kv_size(*introns), 1);
    ck_assert_int_eq(kv_A(*introns, 0).is_used, 0);
    ck_assert_ptr_eq(kv_A(*introns, 0).donor, donor);
    ck_assert_ptr_eq(kv_A(*introns, 0).acceptor, &kv_A(t->exons, 1)->splice_sites->canonical_acceptors[0]);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    destroy_transcript_with_sites(t);
} END_TEST

// A cryptic donor that is used gives a second way to splice the intron, shorter than the canonical one
START_TEST(check_build_noncoding_regions__a_used_cryptic_donor_adds_a_region) {
    Transcript *t = two_exon_transcript(STRAND_FWD);
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Record *r = new_record(129, "N", "x");

    SpliceSites *sites = kv_A(t->exons, 0)->splice_sites;
    Records *cryptic_records = malloc(sizeof(Records)); kv_init(*cryptic_records);
    // A cryptic donor at 125, for the canonical donor at 129
    sites->cryptic_donors[sites->n_cryptic_donors++] = (SpliceSite) { cryptic_records, 125, 129, 1, 0, 0, 0.0, 0.5, 0.5 };

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));

    const NonCodingRegions *introns = &kv_A(groups, 1);
    ck_assert_int_eq(kv_size(*introns), 2);

    ck_assert_int_eq(kv_A(*introns, 0).coords.start, 130); // The canonical donor first
    ck_assert_int_eq(kv_A(*introns, 0).coords.end, 151);
    ck_assert_int_eq(kv_A(*introns, 1).is_used, 1);
    ck_assert_int_eq(kv_A(*introns, 1).coords.start, 126); // The cryptic donor splices from base 125
    ck_assert_int_eq(kv_A(*introns, 1).coords.end, 151);
    ck_assert_int_eq(kv_A(*introns, 1).canonical_coords.start, 130); // Compared with the canonical intron
    ck_assert_int_eq(kv_A(*introns, 1).canonical_coords.end, 151);
    ck_assert_ptr_eq(kv_A(*introns, 1).donor, &sites->cryptic_donors[0]);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    destroy_transcript_with_sites(t);
} END_TEST

// On the reverse strand the groups are in transcription order: the 5' UTR is the high one
START_TEST(check_build_noncoding_regions__reverse_strand_runs_in_transcription_order) {
    Transcript *t = two_exon_transcript(STRAND_REV);
    Gene *gene = new_gene(100, 200, STRAND_REV);
    Record *r = new_record(120, "N", "x");

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));
    ck_assert_int_eq(kv_size(groups), 3);

    const NonCodingRegion utr5 = kv_A(kv_A(groups, 0), 0);
    ck_assert_int_eq(utr5.is_utr, UTR5);
    ck_assert_int_eq(utr5.coords.start, 190);
    ck_assert_int_eq(utr5.coords.end, 200);

    const NonCodingRegion intron = kv_A(kv_A(groups, 1), 0);
    ck_assert_int_eq(intron.coords.start, 130);
    ck_assert_int_eq(intron.coords.end, 151);
    ck_assert_ptr_eq(intron.donor, &kv_A(t->exons, 1)->splice_sites->canonical_donors[0]);
    ck_assert_ptr_eq(intron.acceptor, &kv_A(t->exons, 0)->splice_sites->canonical_acceptors[0]);

    const NonCodingRegion utr3 = kv_A(kv_A(groups, 2), 0);
    ck_assert_int_eq(utr3.is_utr, UTR3);
    ck_assert_int_eq(utr3.coords.start, 100);
    ck_assert_int_eq(utr3.coords.end, 110);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    destroy_transcript_with_sites(t);
} END_TEST

START_TEST(check_build_noncoding_regions__one_coding_exon_has_only_the_utrs) {
    Transcript *t = one_exon_transcript();
    Gene *gene = new_gene(100, 160, STRAND_FWD);
    Record *r = new_record(120, "N", "x");

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));
    ck_assert_int_eq(kv_size(groups), 2);
    ck_assert_int_eq(kv_A(kv_A(groups, 0), 0).coords.end, 110);
    ck_assert_int_eq(kv_A(kv_A(groups, 1), 0).coords.start, 149);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    Transcript_destroy(t);
} END_TEST

START_TEST(check_build_noncoding_regions__a_transcript_without_a_coding_exon) {
    Range exons[] = { {100, 159} };
    TranscriptData data = { .id = 1, .start = 100, .end = 159, .strand = STRAND_FWD,
                            .cds = NULL, .constitutive = NULL, .ncds = 0, .exons = exons, .nexons = 1 };
    Transcript *t = Transcript_build(&data);
    Gene *gene = new_gene(100, 160, STRAND_FWD);
    Record *r = new_record(120, "N", "x");

    NonCodingRegionSets groups;
    ck_assert(!build_noncoding_regions(t, gene, r, &groups));
    ck_assert_int_eq(kv_size(groups), 0);

    free_noncoding_region_groups(&groups);
    Record_destroy(r);
    free_gene(gene);
    Transcript_destroy(t);
} END_TEST

/*
* collect_splice_shifts: the shifts of the splice sites of the regions of an isoform. The regions here: the 5' UTR [100, 110),
* an intron [130, 150) and the 3' UTR [190, 200) of a gene [100, 200), so that the bases 110 to 129 are the first coding bases.
*/

static SpliceSite site_of(Records *records, const uint32_t position, const uint32_t canonical, const float posterior) {
    return (SpliceSite) { records, position, canonical, 1, 1, 0, 1.0f, 1.0f, posterior };
}

static NonCodingRegions regions_with_intron(const int is_used, const Range intron, const Range canonical, SpliceSite *donor, SpliceSite *acceptor, const int uses_competitor) {
    NonCodingRegions ncrs; kv_init(ncrs);

    NonCodingRegion utr5 = { 1, { 100, 110 }, { 100, 110 }, NULL, NULL, UTR5, 0 };
    NonCodingRegion region = { is_used, intron, canonical, donor, acceptor, 0, uses_competitor };
    NonCodingRegion utr3 = { 1, { 190, 200 }, { 190, 200 }, NULL, NULL, UTR3, 0 };
    kv_push(NonCodingRegion, ncrs, utr5);
    kv_push(NonCodingRegion, ncrs, region);
    kv_push(NonCodingRegion, ncrs, utr3);

    return ncrs;
}

START_TEST(check_collect_splice_shifts__a_lost_donor_is_an_intron_retention_at_the_donor_base) {
    Transcript t = { .strand = STRAND_FWD };
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Records hap_in; kv_init(hap_in);
    Record *r = new_record(129, "N", "x");

    Records donor_records; kv_init(donor_records); kv_push(Record *, donor_records, r);
    Records acceptor_records; kv_init(acceptor_records);
    SpliceSite donor = site_of(&donor_records, 129, 129, 0.0f);
    SpliceSite acceptor = site_of(&acceptor_records, 150, 150, 1.0f);
    NonCodingRegions ncrs = regions_with_intron(0, (Range) { 130, 150 }, (Range) { 130, 150 }, &donor, &acceptor, 0);

    IsoformVariants shifts; kv_init(shifts);
    collect_splice_shifts(&ncrs, &hap_in, &t, gene, &shifts);

    ck_assert_int_eq(kv_size(shifts), 1);
    const IsoformVariant shift = kv_A(shifts, 0);
    ck_assert_int_eq(shift.shift_type, SPLICE);
    ck_assert_int_eq(shift.coding, 1);
    ck_assert_int_eq(shift.partial, 0);
    ck_assert_int_eq(shift.seq_pos, 19); // The donor base 129 comes after the 19 coding bases 110 to 128
    ck_assert_int_eq(shift.shift, -20);  // The 20 bases of the intron are added
    ck_assert_int_eq(shift.modified_bases, -20);
    ck_assert_ptr_eq(shift.r, r);

    kv_destroy(shifts); kv_destroy(ncrs); kv_destroy(donor_records); kv_destroy(acceptor_records); kv_destroy(hap_in);
    Record_destroy(r);
    free_gene(gene);
} END_TEST

START_TEST(check_collect_splice_shifts__both_sites_lost_gives_a_shift_for_each_record) {
    Transcript t = { .strand = STRAND_FWD };
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Records hap_in; kv_init(hap_in);
    Record *donor_record = new_record(129, "N", "x");
    Record *acceptor_record = new_record(150, "N", "x");

    Records donor_records; kv_init(donor_records); kv_push(Record *, donor_records, donor_record);
    Records acceptor_records; kv_init(acceptor_records); kv_push(Record *, acceptor_records, acceptor_record);
    SpliceSite donor = site_of(&donor_records, 129, 129, 0.0f);
    SpliceSite acceptor = site_of(&acceptor_records, 150, 150, 0.0f);
    NonCodingRegions ncrs = regions_with_intron(0, (Range) { 130, 150 }, (Range) { 130, 150 }, &donor, &acceptor, 0);

    IsoformVariants shifts; kv_init(shifts);
    collect_splice_shifts(&ncrs, &hap_in, &t, gene, &shifts);

    ck_assert_int_eq(kv_size(shifts), 2);
    ck_assert_ptr_eq(kv_A(shifts, 0).r, donor_record);
    ck_assert_ptr_eq(kv_A(shifts, 1).r, acceptor_record);
    for (int i = 0; i < 2; i++) {
        ck_assert_int_eq(kv_A(shifts, i).seq_pos, 19); // Both are placed at the donor base
        ck_assert_int_eq(kv_A(shifts, i).shift, -20);
    }

    kv_destroy(shifts); kv_destroy(ncrs); kv_destroy(donor_records); kv_destroy(acceptor_records); kv_destroy(hap_in);
    Record_destroy(donor_record); Record_destroy(acceptor_record);
    free_gene(gene);
} END_TEST

// A donor that moved 3 bases down splices a shorter intron: the exon gets 3 bases
START_TEST(check_collect_splice_shifts__a_moved_donor_shifts_by_the_distance) {
    Transcript t = { .strand = STRAND_FWD };
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Records hap_in; kv_init(hap_in);
    Record *r = new_record(129, "N", "x");

    Records donor_records; kv_init(donor_records); kv_push(Record *, donor_records, r);
    Records acceptor_records; kv_init(acceptor_records);
    SpliceSite donor = site_of(&donor_records, 132, 129, 0.7f);
    SpliceSite acceptor = site_of(&acceptor_records, 150, 150, 1.0f);
    NonCodingRegions ncrs = regions_with_intron(1, (Range) { 133, 150 }, (Range) { 130, 150 }, &donor, &acceptor, 0);

    IsoformVariants shifts; kv_init(shifts);
    collect_splice_shifts(&ncrs, &hap_in, &t, gene, &shifts);

    ck_assert_int_eq(kv_size(shifts), 1); // The acceptor has no records
    const IsoformVariant shift = kv_A(shifts, 0);
    ck_assert_int_eq(shift.shift, 3);
    ck_assert_int_eq(shift.modified_bases, 3);
    ck_assert_int_eq(shift.partial, 1); // The site is not a competitor of the canonical one
    ck_assert_int_eq(shift.seq_pos, 22); // The new donor base 132: the coding bases 110 to 131 come before it
    ck_assert_ptr_eq(shift.r, r);

    kv_destroy(shifts); kv_destroy(ncrs); kv_destroy(donor_records); kv_destroy(acceptor_records); kv_destroy(hap_in);
    Record_destroy(r);
    free_gene(gene);
} END_TEST

START_TEST(check_collect_splice_shifts__a_used_site_of_a_competitor_is_not_partial) {
    Transcript t = { .strand = STRAND_FWD };
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Records hap_in; kv_init(hap_in);
    Record *r = new_record(129, "N", "x");

    Records donor_records; kv_init(donor_records); kv_push(Record *, donor_records, r);
    Records acceptor_records; kv_init(acceptor_records);
    SpliceSite donor = site_of(&donor_records, 132, 129, 0.7f);
    SpliceSite acceptor = site_of(&acceptor_records, 150, 150, 1.0f);
    NonCodingRegions ncrs = regions_with_intron(1, (Range) { 133, 150 }, (Range) { 130, 150 }, &donor, &acceptor, 1);

    IsoformVariants shifts; kv_init(shifts);
    collect_splice_shifts(&ncrs, &hap_in, &t, gene, &shifts);

    ck_assert_int_eq(kv_size(shifts), 1);
    ck_assert_int_eq(kv_A(shifts, 0).partial, 0);

    kv_destroy(shifts); kv_destroy(ncrs); kv_destroy(donor_records); kv_destroy(acceptor_records); kv_destroy(hap_in);
    Record_destroy(r);
    free_gene(gene);
} END_TEST

START_TEST(check_collect_splice_shifts__regions_without_sites_give_nothing) {
    Transcript t = { .strand = STRAND_FWD };
    Gene *gene = new_gene(100, 200, STRAND_FWD);
    Records hap_in; kv_init(hap_in);

    NonCodingRegions ncrs = regions_with_intron(1, (Range) { 130, 150 }, (Range) { 130, 150 }, NULL, NULL, 0);

    IsoformVariants shifts; kv_init(shifts);
    collect_splice_shifts(&ncrs, &hap_in, &t, gene, &shifts);

    ck_assert_int_eq(kv_size(shifts), 0);

    kv_destroy(shifts); kv_destroy(ncrs); kv_destroy(hap_in);
    free_gene(gene);
} END_TEST

/*
* evaluate_isoform: the rescue and the PTV of one isoform of a record.
*/

// The pieces evaluate_isoform needs, for a gene [100, 160) with one coding exon
typedef struct EvaluateOneExon {
    bcf_hdr_t *hdr;
    Transcript *t;
    Gene *gene;
    NonCodingRegionSets groups;
    NonCodingRegionSets sets;
    Records hap_in;
    kstring_t hap_seq;
} EvaluateOneExon;

static void evaluate_setup(EvaluateOneExon *e, Record *first, Record *second) {
    e->hdr = mock_hdr();
    e->t = one_exon_transcript();
    e->gene = new_gene(100, 160, STRAND_FWD);

    ck_assert(build_noncoding_regions(e->t, e->gene, first, &e->groups));
    enumerate_region_combinations(e->groups.a, kv_size(e->groups), &e->sets);
    ck_assert_int_eq(kv_size(e->sets), 1);

    kv_init(e->hap_in);
    kv_push(Record *, e->hap_in, first);
    if (second != NULL) kv_push(Record *, e->hap_in, second);
    get_haplotype_sequence(e->gene->seq, &e->hap_in, e->gene->coords, &e->hap_seq);
}

static void evaluate_teardown(EvaluateOneExon *e) {
    NonCodingRegionSets_destroy(&e->sets);
    free_noncoding_region_groups(&e->groups);
    kv_destroy(e->hap_in);
    free(e->hap_seq.s);
    free_gene(e->gene);
    Transcript_destroy(e->t);
    bcf_hdr_destroy(e->hdr);
}

// A single deletion of one base in the CDS shifts the frame, and nothing restores it
START_TEST(check_evaluate_isoform__a_frameshift_deletion_alone_is_not_rescued) {
    Record *r = new_record(120, "NN", "N"); // Deletes the base 121
    EvaluateOneExon e;
    evaluate_setup(&e, r, NULL);

    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;
    const IsoformContext ctx = { e.hdr, 0, r, 0, e.t, e.gene, &csq, &e.hap_in, &e.hap_seq };

    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(e.sets, 0));

    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.rescue_type, NO_RESCUE);
    ck_assert_int_eq(isoform.affected_bases, 0);
    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);
    ck_assert_float_eq(isoform.prior_probability, 1.0f);
    ck_assert_float_eq(isoform.posterior_probability, 0.0f);

    // The CDS has 39 bases, so the exon is in frame, and it is constitutive. The frameshift is unrescued and the sequence
    // ("CGA" repeated, no stop in any frame) has no stop to read on to, so the isoform has a PTV for NMD but no escape region
    ck_assert_uint_eq(r->flags[0], INFRAME);
    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    evaluate_teardown(&e);
    Record_destroy(r);
} END_TEST

// A deletion of one base and an insertion of one base 5 bases further: the frame is back
START_TEST(check_evaluate_isoform__an_insertion_after_a_deletion_restores_the_frame) {
    Record *r = new_record(120, "NN", "N");      // Deletes the base 121
    Record *other = new_record(125, "N", "Nx");  // Inserts a base after 125
    EvaluateOneExon e;
    evaluate_setup(&e, r, other);

    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;
    const IsoformContext ctx = { e.hdr, 0, r, 0, e.t, e.gene, &csq, &e.hap_in, &e.hap_seq };

    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(e.sets, 0));

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.rescue_type, FRAME_RESTORING_INDEL);
    ck_assert_int_eq(isoform.affected_bases, 15); // From the start of the sequence up to the insertion, and the inserted base
    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);

    evaluate_teardown(&e);
    Record_destroy(r); Record_destroy(other);
} END_TEST

// Two SNVs two bases apart, 113 C>T and 114 C>A, that make the stop codon TAA out of CCA only together. The reading frame conflict
// is found for phase 0 (reading_frame[0] of the consequence), and the isoform is evaluated for the given phase.
static Isoform evaluate_reading_frame_pair(const int phase) {
    Record *r = new_record(113, "C", "T");
    Record *other = new_record(114, "C", "A");
    other->tscripts = malloc(sizeof(Transcripts)); kv_init(*other->tscripts);
    kv_push(Transcript *, *other->tscripts, NULL);

    EvaluateOneExon e;
    evaluate_setup(&e, r, other);
    set_bases(e.gene, 113, "CCA");
    free(e.hap_seq.s);
    get_haplotype_sequence(e.gene->seq, &e.hap_in, e.gene->coords, &e.hap_seq);

    TscriptCsqs csq = {0};
    csq.reading_frame[0] = 1;
    const IsoformContext ctx = { e.hdr, phase, r, 0, e.t, e.gene, &csq, &e.hap_in, &e.hap_seq };

    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(e.sets, 0));

    evaluate_teardown(&e);
    kv_destroy(*other->tscripts); free(other->tscripts); other->tscripts = NULL;
    Record_destroy(r); Record_destroy(other);

    return isoform;
}

START_TEST(check_evaluate_isoform__a_reading_frame_conflict_that_makes_a_stop_gives_a_stop_gained) {
    const Isoform isoform = evaluate_reading_frame_pair(0);

    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_int_eq(isoform.ptv_type, STOP_GAINED);
    ck_assert_int_eq(isoform.rescue, 0);
} END_TEST

// In the phase without the conflict the trigger is not there, so the stop is not looked for
START_TEST(check_evaluate_isoform__the_reading_frame_trigger_of_another_phase_is_not_used) {
    const Isoform isoform = evaluate_reading_frame_pair(1);

    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);
} END_TEST

/*
* A stop that VEP found (STOP_RESCUE_VALUE) is looked for in the codons the change touches. The gene is [100, 160) with the CDS
* [110, 149): the coding index of a base is its position minus 110, and the codons are the indices 0 to 2, 3 to 5, and so on. If
* the stop is not found, the isoform does not have it and the variant counts as rescued; if it is, it is not rescued.
*
* A change that begins late in a codon and runs into the next one has to have that codon looked at too: the window runs from the
* start of the codon of the first changed base up to the last changed base, not for the length of the change from the codon start.
*/
static Isoform evaluate_vep_stop(Record *r, const int64_t bases_at, const char *bases) {
    EvaluateOneExon e;
    evaluate_setup(&e, r, NULL);
    set_bases(e.gene, bases_at, bases);
    free(e.hap_seq.s);
    get_haplotype_sequence(e.gene->seq, &e.hap_in, e.gene->coords, &e.hap_seq);

    TscriptCsqs csq = {0};
    csq.triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;
    const IsoformContext ctx = { e.hdr, 0, r, 0, e.t, e.gene, &csq, &e.hap_in, &e.hap_seq };

    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(e.sets, 0));

    evaluate_teardown(&e);
    return isoform;
}

// A substitution of one base never leaves its codon
START_TEST(check_evaluate_isoform__a_vep_stop_from_a_substitution_is_not_rescued) {
    Record *r = new_record(113, "C", "T");   // The codon CAA at 3 to 5 becomes TAA

    const Isoform isoform = evaluate_vep_stop(r, 110, "CGACAA");

    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.rescue_type, NO_RESCUE);

    Record_destroy(r);
} END_TEST

// The two bases at the indices 2 and 3, the last of the codon CGA and the first of CAA, become G and T: CGG TAA. The change has
// 2 bases, and a window of that length from the start of the codon at 0 holds only CGG
START_TEST(check_evaluate_isoform__a_vep_stop_made_by_an_mnv_that_runs_into_the_next_codon_is_not_rescued) {
    Record *r = new_record(112, "AC", "GT");

    const Isoform isoform = evaluate_vep_stop(r, 110, "CGACAA");

    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.rescue_type, NO_RESCUE);
    ck_assert_int_eq(isoform.ptv, 0);

    Record_destroy(r);
} END_TEST

// An insertion of TA after the index 2, in front of the A at the index 3: CGA TAA A. The attached base and the 2 inserted ones are
// 3 bases, which from the start of the codon at 0 holds only CGA
START_TEST(check_evaluate_isoform__a_vep_stop_made_by_an_insertion_that_runs_into_the_next_codon_is_not_rescued) {
    Record *r = new_record(112, "A", "ATA");

    const Isoform isoform = evaluate_vep_stop(r, 110, "CGAAAA");

    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.rescue_type, NO_RESCUE);

    Record_destroy(r);
} END_TEST

// Without a stop in the codons of the change the isoform does not have the stop, and the variant is rescued
START_TEST(check_evaluate_isoform__a_vep_stop_without_a_stop_in_the_codons_of_the_change_is_rescued) {
    Record *r = new_record(112, "AC", "GT");

    const Isoform isoform = evaluate_vep_stop(r, 110, "CGACGC");   // CGG TGC: no stop

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.rescue_type, COMPLETE_RESCUE);

    Record_destroy(r);
} END_TEST

/*
* A retained intron in an isoform, from a lost donor. The intron is the 21 bases [130, 151), which are in the isoform once and
* directly after the first coding exon: the coding sequence is the 20 bases of the first CDS, the 21 of the intron, and the 39
* of the second CDS. The retention is placed at the donor base 129, which is index 19 of the coding sequence, and the check
* for a stop codon covers the codons from the one holding index 19 up to the last intron base. A TAG in frame in the intron
* (the intron starts at index 20, so its codons start at index 21) makes a stop; the same TAG one base further is out of frame.
*/

typedef struct Retention {
    Transcript *t;
    Gene *gene;
    Record *r;
    bcf_hdr_t *hdr;
} Retention;

static Retention retention_setup(const int64_t stop_at) {
    Retention x;
    x.hdr = mock_hdr();
    x.t = two_exon_transcript(STRAND_FWD);
    x.gene = new_gene(100, 200, STRAND_FWD);
    if (stop_at > 0) set_bases(x.gene, stop_at, "TAG");
    x.r = new_record(129, "N", "x");

    // The donor is lost by this record
    SpliceSite *donor = &kv_A(x.t->exons, 0)->splice_sites->canonical_donors[0];
    donor->predicted = 1;
    donor->posterior_probability = 0;
    kv_push(Record *, *(Records *) donor->records_ptr, x.r);

    return x;
}

static void retention_teardown(Retention *x) {
    destroy_transcript_with_sites(x->t);
    free_gene(x->gene);
    Record_destroy(x->r);
    bcf_hdr_destroy(x->hdr);
}

static Isoform evaluate_retention(Retention *x) {
    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(x->t, x->gene, x->r, &groups));

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups.a, kv_size(groups), &sets);
    ck_assert_int_eq(kv_size(sets), 1);

    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, x->r);
    kstring_t hap_seq;
    get_haplotype_sequence(x->gene->seq, &hap_in, x->gene->coords, &hap_seq);

    TscriptCsqs csq = {0};
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
    const IsoformContext ctx = { x->hdr, 0, x->r, 0, x->t, x->gene, &csq, &hap_in, &hap_seq };

    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(sets, 0));

    NonCodingRegionSets_destroy(&sets);
    free_noncoding_region_groups(&groups);
    kv_destroy(hap_in);
    free(hap_seq.s);

    return isoform;
}

START_TEST(check_evaluate_isoform__a_retained_intron_without_a_stop_is_a_rescue) {
    Retention x = retention_setup(0);

    const Isoform isoform = evaluate_retention(&x);

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.rescue_type, INTRON_RETENTION);
    ck_assert_int_eq(isoform.affected_bases, 21);
    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);

    retention_teardown(&x);
} END_TEST

// The TAG at 131 to 133 is the codon at index 21 of the isoform: in frame
START_TEST(check_evaluate_isoform__a_retained_intron_with_an_in_frame_stop_gives_a_stop_gained) {
    Retention x = retention_setup(131);

    const Isoform isoform = evaluate_retention(&x);

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.rescue_type, INTRON_RETENTION);
    ck_assert_int_eq(isoform.affected_bases, 21);
    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_int_eq(isoform.ptv_type, STOP_GAINED);

    retention_teardown(&x);
} END_TEST

// The TAG at 132 to 134 starts one base after a codon starts: out of frame, so no stop
START_TEST(check_evaluate_isoform__a_retained_intron_with_an_out_of_frame_stop_has_no_stop) {
    Retention x = retention_setup(132);

    const Isoform isoform = evaluate_retention(&x);

    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);

    retention_teardown(&x);
} END_TEST

// The codons checked run from the one that holds index 19 up to the last intron base: the indexes 18, 21, ..., 39. The last one
// is the bases 149 and 150 of the intron and the base 151, the first base of the second CDS. A TAG there is found.
START_TEST(check_evaluate_isoform__a_stop_in_the_last_codon_of_the_window_is_found) {
    Retention x = retention_setup(149);

    const Isoform isoform = evaluate_retention(&x);

    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_int_eq(isoform.ptv_type, STOP_GAINED);

    retention_teardown(&x);
} END_TEST

// The next codon (index 42, the bases 152 to 154) is past the window
START_TEST(check_evaluate_isoform__a_stop_after_the_window_is_not_found) {
    Retention x = retention_setup(152);

    const Isoform isoform = evaluate_retention(&x);

    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.ptv_type, NO_PTV);

    retention_teardown(&x);
} END_TEST

/*
* The NMD flag of an unrescued or induced PTV. The gene [100, 951) has three coding exons [100, 400), [500, 700) and [800, 951) with
* no UTR, so 651 coding bases (217 codons): the indices 0 to 299, 300 to 499 and 500 to 650. The coding sequence is GCC repeated
* with TAA at the end, so no frame has a stop before its own, and the introns are ACG repeated, which has none either. A test
* plants the stops it needs. The escape regions are the coding indices 0 to 149 (start-proximal) and 450 to 650 (the last 50 bases
* of the second exon and all of the last, distal), and a retained intron between the last two exons (the indices after 499 when it
* is kept, distal too); no exon is long enough for a long exon region, and the rest is the NMD target region. The isoform keeps the
* flag; build_isoforms_for_phase sets it on the record.
*/
#define NMD_BITS (NMD_START_PROXIMAL | NMD_LONG_EXON | NMD_DISTAL)

static int64_t coding_to_genomic(const int64_t c) {
    return c < 300 ? 100 + c : c < 500 ? 500 + (c - 300) : 800 + (c - 500);
}

static void set_coding_bases(Gene *gene, const int64_t index, const char *bases) {
    for (size_t i = 0; i < strlen(bases); i++) gene->seq.s[coding_to_genomic(index + (int64_t) i) - gene->coords.start] = bases[i];
}

static Transcript *three_exon_transcript(void) {
    Range cds[] = { {100, 399}, {500, 699}, {800, 950} };
    int constitutive[] = { 1, 1, 1 };
    TranscriptData data = { .id = 1, .start = 100, .end = 950, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = cds, .nexons = 3 };

    Transcript *t = Transcript_build(&data);
    for (size_t ei = 0; ei < kv_size(t->exons); ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    return t;
}

static Gene *three_exon_gene(void) {
    Gene *gene = new_gene(100, 951, STRAND_FWD);

    char cds[652];
    for (int c = 0; c < 648; c++) cds[c] = "GCC"[c % 3];
    memcpy(cds + 648, "TAA", 4);
    set_coding_bases(gene, 0, cds);

    return gene;
}

// Evaluates the isoform of r for the phase. other is a second record of the haplotype or NULL, and lost_donor_exon the exon whose
// canonical donor r takes away (so that the intron after it is retained), or -1
static Isoform evaluate_three_exons(Gene *gene, Transcript *t, Record *r, Record *other, TscriptCsqs *csq, const int phase, const int lost_donor_exon) {
    if (lost_donor_exon >= 0) {
        SpliceSite *donor = &kv_A(t->exons, lost_donor_exon)->splice_sites->canonical_donors[0];
        donor->predicted = 1;
        donor->posterior_probability = 0;
        kv_push(Record *, *(Records *) donor->records_ptr, r);
    }

    NonCodingRegionSets groups;
    ck_assert(build_noncoding_regions(t, gene, r, &groups));

    NonCodingRegionSets sets;
    enumerate_region_combinations(groups.a, kv_size(groups), &sets);
    ck_assert_int_eq(kv_size(sets), 1);

    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, r);
    if (other != NULL) kv_push(Record *, hap_in, other);
    kstring_t hap_seq;
    get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_seq);

    bcf_hdr_t *hdr = mock_hdr();
    const IsoformContext ctx = { hdr, phase, r, 0, t, gene, csq, &hap_in, &hap_seq };
    const Isoform isoform = evaluate_isoform(&ctx, &kv_A(sets, 0));

    NonCodingRegionSets_destroy(&sets);
    free_noncoding_region_groups(&groups);
    kv_destroy(hap_in);
    free(hap_seq.s);
    bcf_hdr_destroy(hdr);

    return isoform;
}

// A stop that VEP found: the codon codon_index, planted as CAA, has its C turned into a T. Gives the NMD flag of the isoform
static uint32_t nmd_of_vep_stop(const int codon_index) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_coding_bases(gene, 3 * codon_index, "CAA");
    Record *r = new_record(coding_to_genomic(3 * codon_index), "C", "T");

    TscriptCsqs csq = {0};
    csq.triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;
    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);
    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.nmd_ptv, 1);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);

    return isoform.nmd_flag;
}

START_TEST(check_nmd_flag__a_vep_stop_in_the_start_proximal_region_is_start_proximal) {
    ck_assert_uint_eq(nmd_of_vep_stop(20), NMD_START_PROXIMAL);
    ck_assert_uint_eq(nmd_of_vep_stop(1), NMD_START_PROXIMAL);
} END_TEST

START_TEST(check_nmd_flag__a_vep_stop_in_the_target_region_has_no_flag) {
    ck_assert_uint_eq(nmd_of_vep_stop(70), 0);
    ck_assert_uint_eq(nmd_of_vep_stop(120), 0);   // In the second exon, before its last 50 bases
} END_TEST

START_TEST(check_nmd_flag__a_vep_stop_in_the_last_50_bases_of_the_penultimate_exon_or_the_last_exon_is_distal) {
    ck_assert_uint_eq(nmd_of_vep_stop(160), NMD_DISTAL);
    ck_assert_uint_eq(nmd_of_vep_stop(200), NMD_DISTAL);
} END_TEST

// The stop is placed by the first base of its codon: the codon 49 (147 to 149) begins before 150, and the codon 50 at 150. The
// codon 149 (447 to 449) begins before the last 50 bases of the second exon (450 on), and the codon 150 at 450
START_TEST(check_nmd_flag__the_stop_is_placed_by_the_first_base_of_its_codon) {
    ck_assert_uint_eq(nmd_of_vep_stop(49), NMD_START_PROXIMAL);
    ck_assert_uint_eq(nmd_of_vep_stop(50), 0);
    ck_assert_uint_eq(nmd_of_vep_stop(149), 0);
    ck_assert_uint_eq(nmd_of_vep_stop(150), NMD_DISTAL);
} END_TEST

// evaluate_isoform keeps the flag on the isoform: the record gets it from build_isoforms_for_phase, once for all isoforms
START_TEST(check_nmd_flag__is_kept_on_the_isoform_and_not_set_on_the_record) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_coding_bases(gene, 30, "CAA");
    Record *r = new_record(coding_to_genomic(30), "C", "T");
    TscriptCsqs csq = {0};
    csq.triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 1, -1);

    ck_assert_uint_eq(isoform.nmd_flag, NMD_START_PROXIMAL);
    ck_assert_uint_eq(r->flags[0] & NMD_BITS, 0);
    ck_assert_uint_eq(r->flags[1] & NMD_BITS, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// A stop VEP found that the isoform does not have (the codon is not a stop in it) is rescued: no PTV, no flag
START_TEST(check_nmd_flag__a_rescued_vep_stop_has_no_flag) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    Record *r = new_record(coding_to_genomic(210), "G", "A");   // GCC to ACC
    TscriptCsqs csq = {0};
    csq.triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.nmd_ptv, 0);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// Two SNVs that make a stop only together (the reading frame trigger): an induced PTV, with its stop at the codon of the pair
START_TEST(check_nmd_flag__a_stop_made_by_two_variants_together_is_located_at_the_codon_of_the_pair) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_coding_bases(gene, 300, "CCA");          // The first codon of the second exon
    Record *r = new_record(coding_to_genomic(300), "C", "T");
    Record *other = new_record(coding_to_genomic(301), "C", "A");
    other->tscripts = malloc(sizeof(Transcripts)); kv_init(*other->tscripts);
    kv_push(Transcript *, *other->tscripts, NULL);

    TscriptCsqs csq = {0};
    csq.reading_frame[0] = 1;
    const Isoform isoform = evaluate_three_exons(gene, t, r, other, &csq, 0, -1);

    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_int_eq(isoform.ptv_type, STOP_GAINED);
    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);   // The stop is at the coding index 300: in the target region

    destroy_transcript_with_sites(t); free_gene(gene);
    kv_destroy(*other->tscripts); free(other->tscripts); other->tscripts = NULL;
    Record_destroy(r); Record_destroy(other);
} END_TEST

// A VEP frameshift that nothing restores is unrescued, though the isoform does not call it a PTV. The deletion of the base 210
// (genomic 310) of GCC CTA AGC leaves CCC TAA GCC: the stop is in the next codon, at 213, which the codon of the change does not
// have. The base it came from is 214: the NMD target region
START_TEST(check_nmd_flag__an_unrescued_frameshift_is_located_by_reading_on_to_the_next_stop) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_coding_bases(gene, 213, "CTAAGC");
    Record *r = new_record(coding_to_genomic(210) - 1, "CG", "C");
    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);

    ck_assert_int_eq(isoform.rescue, 0);
    ck_assert_int_eq(isoform.ptv, 0);
    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// The same frameshift near the start of the coding sequence puts the stop in the start-proximal region
START_TEST(check_nmd_flag__an_unrescued_frameshift_with_a_stop_near_the_start_is_start_proximal) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_coding_bases(gene, 33, "CTAAGC");   // The codons 11 and 12
    Record *r = new_record(coding_to_genomic(30) - 1, "CG", "C");   // Deletes the base 30, the first of the codon 10
    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);

    ck_assert_uint_eq(isoform.nmd_flag, NMD_START_PROXIMAL);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// A frameshift with no stop before the end of the sequence: no stop can be found, and the original stop codon is not one. The
// isoform still has a PTV for NMD, without an escape region
START_TEST(check_nmd_flag__a_frameshift_without_a_new_stop_has_no_flag) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    Record *r = new_record(coding_to_genomic(609) - 1, "CG", "C");
    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);

    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// A deletion of three bases is in frame and rescued: not a PTV
START_TEST(check_nmd_flag__an_in_frame_deletion_has_no_flag) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    Record *r = new_record(coding_to_genomic(210) - 1, "CGCC", "C");
    TscriptCsqs csq = {0};
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, -1);

    ck_assert_int_eq(isoform.rescue, 1);
    ck_assert_int_eq(isoform.nmd_ptv, 0);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// A lost donor retains the 100 bases of the intron after the first exon: a frameshift of the isoform (100 is not a multiple of 3),
// the PTV is induced. Reading on from the retention with no stop planted finds none
START_TEST(check_nmd_flag__a_retention_that_shifts_the_frame_without_a_stop_has_no_flag) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    Record *r = new_record(399, "C", "T");   // The last base of the first exon
    TscriptCsqs csq = {0};
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, 0);

    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_int_eq(isoform.ptv_type, FRAMESHIFT);
    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// A stop in the retained intron between the first two exons, in the isoform frame (the intron starts at the index 300 and the TAG at
// 303): a retained intron that is not between the last two exons is not an escape region
START_TEST(check_nmd_flag__a_stop_in_a_retained_intron_that_is_not_the_last_has_no_flag) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_bases(gene, 403, "TAG");
    Record *r = new_record(399, "C", "T");
    TscriptCsqs csq = {0};
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, 0);

    ck_assert_int_eq(isoform.nmd_ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, 0);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

// The retained intron between the last two exons is distal: the intron after the second exon starts at the index 500, and the TAG
// at 704 to 706 is in frame at 504
START_TEST(check_nmd_flag__a_stop_in_the_retained_intron_between_the_last_two_exons_is_distal) {
    Gene *gene = three_exon_gene();
    Transcript *t = three_exon_transcript();
    set_bases(gene, 704, "TAG");
    Record *r = new_record(699, "C", "T");   // The last base of the second exon
    TscriptCsqs csq = {0};
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;

    const Isoform isoform = evaluate_three_exons(gene, t, r, NULL, &csq, 0, 1);

    ck_assert_int_eq(isoform.ptv, 1);
    ck_assert_uint_eq(isoform.nmd_flag, NMD_DISTAL);

    destroy_transcript_with_sites(t); free_gene(gene); Record_destroy(r);
} END_TEST

/*
* nmd_flag_of_isoforms: the one NMD flag of a transcript and phase, from the most probable isoform with a PTV for NMD.
*/

// Isoforms with the given posteriors, PTVs for NMD and NMD flags. Freed with kv_destroy of the isoforms vector
static Isoforms nmd_isoforms(const float *posteriors, const uint8_t *nmd_ptvs, const uint32_t *nmd_flags, const int n) {
    Isoforms isoforms = { .valid = 1 };
    kv_init(isoforms.isoforms);
    for (int i = 0; i < n; i++) {
        const Isoform isoform = { .posterior_probability = posteriors[i], .nmd_ptv = nmd_ptvs[i], .nmd_flag = nmd_flags[i] };
        kv_push(Isoform, isoforms.isoforms, isoform);
    }

    return isoforms;
}

START_TEST(check_nmd_flag_of_isoforms__the_most_probable_isoform_with_a_ptv_decides) {
    const float posteriors[] = { 0.2f, 0.5f, 0.3f };
    const uint8_t nmd_ptvs[] = { 1, 1, 1 };
    const uint32_t nmd_flags[] = { NMD_START_PROXIMAL, NMD_DISTAL, NMD_LONG_EXON };
    Isoforms isoforms = nmd_isoforms(posteriors, nmd_ptvs, nmd_flags, 3);

    ck_assert_uint_eq(nmd_flag_of_isoforms(&isoforms), NMD_DISTAL);

    kv_destroy(isoforms.isoforms);
} END_TEST

// The most probable isoform with a PTV has its stop outside every escape region: no flag, though a less probable one escapes
START_TEST(check_nmd_flag_of_isoforms__no_escape_on_the_most_probable_isoform_gives_no_flag) {
    const float posteriors[] = { 0.3f, 0.7f };
    const uint8_t nmd_ptvs[] = { 1, 1 };
    const uint32_t nmd_flags[] = { NMD_DISTAL, 0 };
    Isoforms isoforms = nmd_isoforms(posteriors, nmd_ptvs, nmd_flags, 2);

    ck_assert_uint_eq(nmd_flag_of_isoforms(&isoforms), 0);

    kv_destroy(isoforms.isoforms);
} END_TEST

// An isoform without a PTV (a rescue) has nothing to say about NMD, however probable it is
START_TEST(check_nmd_flag_of_isoforms__isoforms_without_a_ptv_are_ignored) {
    const float posteriors[] = { 0.9f, 0.1f };
    const uint8_t nmd_ptvs[] = { 0, 1 };
    const uint32_t nmd_flags[] = { 0, NMD_LONG_EXON };
    Isoforms isoforms = nmd_isoforms(posteriors, nmd_ptvs, nmd_flags, 2);

    ck_assert_uint_eq(nmd_flag_of_isoforms(&isoforms), NMD_LONG_EXON);

    kv_destroy(isoforms.isoforms);
} END_TEST

START_TEST(check_nmd_flag_of_isoforms__a_tie_goes_to_the_first_and_no_ptv_gives_no_flag) {
    const float posteriors[] = { 0.5f, 0.5f };
    const uint8_t nmd_ptvs[] = { 1, 1 };
    const uint32_t nmd_flags[] = { NMD_START_PROXIMAL, NMD_DISTAL };
    Isoforms tie = nmd_isoforms(posteriors, nmd_ptvs, nmd_flags, 2);
    ck_assert_uint_eq(nmd_flag_of_isoforms(&tie), NMD_START_PROXIMAL);
    kv_destroy(tie.isoforms);

    const uint8_t no_ptvs[] = { 0, 0 };
    Isoforms rescued = nmd_isoforms(posteriors, no_ptvs, nmd_flags, 2);
    ck_assert_uint_eq(nmd_flag_of_isoforms(&rescued), 0);
    kv_destroy(rescued.isoforms);

    Isoforms empty = nmd_isoforms(NULL, NULL, NULL, 0);
    ck_assert_uint_eq(nmd_flag_of_isoforms(&empty), 0);
    kv_destroy(empty.isoforms);
} END_TEST

/*
* build_isoforms_for_phase: the whole stage for one record.
*/

START_TEST(check_build_isoforms_for_phase__a_retained_intron_with_a_stop) {
    Retention x = retention_setup(131);

    // The record as the loader and the earlier stages leave it
    x.r->tscripts = malloc(sizeof(Transcripts)); kv_init(*x.r->tscripts);
    kv_push(Transcript *, *x.r->tscripts, x.t);
    x.r->transcript_genes = malloc(sizeof(Genes)); kv_init(*x.r->transcript_genes);
    kv_push(Gene *, *x.r->transcript_genes, x.gene);
    x.r->genes = malloc(sizeof(Genes)); kv_init(*x.r->genes);
    kv_push(Gene *, *x.r->genes, x.gene);

    x.r->isoforms = malloc(NUM_PHASES * sizeof(Isoforms));
    for (int k = 0; k < NUM_PHASES; k++) x.r->isoforms[k].valid = 0;

    TscriptCsqs csq = {0};
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
    x.r->vep = malloc(sizeof(VepCsqs)); kv_init(*x.r->vep);
    kv_push(TscriptCsqs *, *x.r->vep, &csq);

    Session s = {0};
    s.hdr = x.hdr;
    s.prediction_radius = 50;

    Records records; kv_init(records);
    kv_push(Record *, records, x.r);

    build_isoforms_for_phase(&s, 0, &records);

    const Isoforms *isoforms = &x.r->isoforms[0];
    ck_assert_int_eq(isoforms->valid, 1);
    ck_assert_int_eq(kv_size(isoforms->isoforms), 1);
    ck_assert_float_eq(kv_A(isoforms->isoforms, 0).posterior_probability, 1.0f);
    ck_assert_int_eq(isoforms->rescue, 1);
    ck_assert_float_eq(isoforms->rescue_probability, 1.0f);
    ck_assert_int_eq(isoforms->n_rescue_types, 1);
    ck_assert_int_eq(isoforms->rescue_type_probabilities->rescue_types[0], INTRON_RETENTION);
    ck_assert_int_eq(isoforms->rescue_type_probabilities->affected_bases[0], 21);
    ck_assert_int_eq(isoforms->ptv, 1);
    ck_assert_int_eq(isoforms->ptv_type, STOP_GAINED);
    ck_assert_float_eq(isoforms->ptv_probability, 1.0f);

    // The stop is in the retained intron between the two exons, which is distal. One NMD flag, on the phase that was built
    ck_assert_uint_eq(x.r->flags[0] & NMD_BITS, NMD_DISTAL);
    ck_assert_uint_eq(x.r->flags[1] & NMD_BITS, 0);

    ck_assert_int_eq(x.r->isoforms[1].valid, 0); // The other phase was not built

    kv_destroy(records);
    kv_destroy(*x.r->vep); free(x.r->vep); x.r->vep = NULL;
    destroy_transcript_with_sites(x.t);
    free_gene(x.gene);
    Record_destroy(x.r); // Also frees the isoforms, the transcript and gene vectors, and the flags
    bcf_hdr_destroy(x.hdr);
} END_TEST

// A record that is too long for the prediction radius, or a transcript without a trigger, gets no isoforms
START_TEST(check_build_isoforms_for_phase__without_a_trigger_nothing_is_built) {
    Retention x = retention_setup(0);

    x.r->tscripts = malloc(sizeof(Transcripts)); kv_init(*x.r->tscripts);
    kv_push(Transcript *, *x.r->tscripts, x.t);
    x.r->transcript_genes = malloc(sizeof(Genes)); kv_init(*x.r->transcript_genes);
    kv_push(Gene *, *x.r->transcript_genes, x.gene);
    x.r->genes = malloc(sizeof(Genes)); kv_init(*x.r->genes);
    kv_push(Gene *, *x.r->genes, x.gene);
    x.r->isoforms = malloc(NUM_PHASES * sizeof(Isoforms));
    for (int k = 0; k < NUM_PHASES; k++) x.r->isoforms[k].valid = 0;

    TscriptCsqs csq = {0}; // No trigger
    x.r->vep = malloc(sizeof(VepCsqs)); kv_init(*x.r->vep);
    kv_push(TscriptCsqs *, *x.r->vep, &csq);

    Session s = {0};
    s.hdr = x.hdr;
    s.prediction_radius = 50;

    Records records; kv_init(records);
    kv_push(Record *, records, x.r);

    build_isoforms_for_phase(&s, 0, &records);

    ck_assert_int_eq(x.r->isoforms[0].valid, 0);

    kv_destroy(records);
    kv_destroy(*x.r->vep); free(x.r->vep); x.r->vep = NULL;
    destroy_transcript_with_sites(x.t);
    free_gene(x.gene);
    Record_destroy(x.r);
    bcf_hdr_destroy(x.hdr);
} END_TEST

Suite *isoforms_suite() {
    Suite *s = suite_create("Isoforms");

    TCase *tc_rescue = mocked_tcase("calculate_rescue_probabilities");
    tcase_add_test(tc_rescue, check_calculate_rescue_probabilities__two_types);
    tcase_add_test(tc_rescue, check_calculate_rescue_probabilities__same_type_and_bases_are_added_up);
    tcase_add_test(tc_rescue, check_calculate_rescue_probabilities__many_types_grow_the_arrays);
    tcase_add_test(tc_rescue, check_calculate_rescue_probabilities__half_is_not_a_rescue);
    suite_add_tcase(s, tc_rescue);

    TCase *tc_ptv = mocked_tcase("calculate_ptv_probabilities");
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__the_largest_type_wins);
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__frameshift_can_win);
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__a_tie_gives_no_type);
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__a_long_inframe_indel_alone);
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__isoforms_without_a_ptv_are_ignored);
    tcase_add_test(tc_ptv, check_calculate_ptv_probabilities__half_is_not_a_ptv);
    suite_add_tcase(s, tc_ptv);

    TCase *tc_type = mocked_tcase("shift_rescue_type");
    tcase_add_test(tc_type, check_shift_rescue_type__sign_and_partial);
    suite_add_tcase(s, tc_type);

    TCase *tc_regions = mocked_tcase("build_noncoding_regions");
    tcase_add_test(tc_regions, check_build_noncoding_regions__a_transcript_without_predictions);
    tcase_add_test(tc_regions, check_build_noncoding_regions__an_intron_with_both_sites_lost_is_retained);
    tcase_add_test(tc_regions, check_build_noncoding_regions__a_lost_donor_retains_the_intron);
    tcase_add_test(tc_regions, check_build_noncoding_regions__a_used_cryptic_donor_adds_a_region);
    tcase_add_test(tc_regions, check_build_noncoding_regions__reverse_strand_runs_in_transcription_order);
    tcase_add_test(tc_regions, check_build_noncoding_regions__one_coding_exon_has_only_the_utrs);
    tcase_add_test(tc_regions, check_build_noncoding_regions__a_transcript_without_a_coding_exon);
    suite_add_tcase(s, tc_regions);

    TCase *tc_shifts = mocked_tcase("collect_splice_shifts");
    tcase_add_test(tc_shifts, check_collect_splice_shifts__a_lost_donor_is_an_intron_retention_at_the_donor_base);
    tcase_add_test(tc_shifts, check_collect_splice_shifts__both_sites_lost_gives_a_shift_for_each_record);
    tcase_add_test(tc_shifts, check_collect_splice_shifts__a_moved_donor_shifts_by_the_distance);
    tcase_add_test(tc_shifts, check_collect_splice_shifts__a_used_site_of_a_competitor_is_not_partial);
    tcase_add_test(tc_shifts, check_collect_splice_shifts__regions_without_sites_give_nothing);
    suite_add_tcase(s, tc_shifts);

    TCase *tc_evaluate = mocked_tcase("evaluate_isoform");
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_frameshift_deletion_alone_is_not_rescued);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__an_insertion_after_a_deletion_restores_the_frame);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_reading_frame_conflict_that_makes_a_stop_gives_a_stop_gained);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__the_reading_frame_trigger_of_another_phase_is_not_used);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_vep_stop_from_a_substitution_is_not_rescued);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_vep_stop_made_by_an_mnv_that_runs_into_the_next_codon_is_not_rescued);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_vep_stop_made_by_an_insertion_that_runs_into_the_next_codon_is_not_rescued);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_vep_stop_without_a_stop_in_the_codons_of_the_change_is_rescued);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_vep_stop_in_the_start_proximal_region_is_start_proximal);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_vep_stop_in_the_target_region_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_vep_stop_in_the_last_50_bases_of_the_penultimate_exon_or_the_last_exon_is_distal);
    tcase_add_test(tc_evaluate, check_nmd_flag__the_stop_is_placed_by_the_first_base_of_its_codon);
    tcase_add_test(tc_evaluate, check_nmd_flag__is_kept_on_the_isoform_and_not_set_on_the_record);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_rescued_vep_stop_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_stop_made_by_two_variants_together_is_located_at_the_codon_of_the_pair);
    tcase_add_test(tc_evaluate, check_nmd_flag__an_unrescued_frameshift_is_located_by_reading_on_to_the_next_stop);
    tcase_add_test(tc_evaluate, check_nmd_flag__an_unrescued_frameshift_with_a_stop_near_the_start_is_start_proximal);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_frameshift_without_a_new_stop_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__an_in_frame_deletion_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_retention_that_shifts_the_frame_without_a_stop_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_stop_in_a_retained_intron_that_is_not_the_last_has_no_flag);
    tcase_add_test(tc_evaluate, check_nmd_flag__a_stop_in_the_retained_intron_between_the_last_two_exons_is_distal);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_retained_intron_without_a_stop_is_a_rescue);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_retained_intron_with_an_in_frame_stop_gives_a_stop_gained);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_retained_intron_with_an_out_of_frame_stop_has_no_stop);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_stop_in_the_last_codon_of_the_window_is_found);
    tcase_add_test(tc_evaluate, check_evaluate_isoform__a_stop_after_the_window_is_not_found);
    suite_add_tcase(s, tc_evaluate);

    TCase *tc_nmd = mocked_tcase("nmd_flag_of_isoforms");
    tcase_add_test(tc_nmd, check_nmd_flag_of_isoforms__the_most_probable_isoform_with_a_ptv_decides);
    tcase_add_test(tc_nmd, check_nmd_flag_of_isoforms__no_escape_on_the_most_probable_isoform_gives_no_flag);
    tcase_add_test(tc_nmd, check_nmd_flag_of_isoforms__isoforms_without_a_ptv_are_ignored);
    tcase_add_test(tc_nmd, check_nmd_flag_of_isoforms__a_tie_goes_to_the_first_and_no_ptv_gives_no_flag);
    suite_add_tcase(s, tc_nmd);

    TCase *tc_stage = mocked_tcase("build_isoforms_for_phase");
    tcase_add_test(tc_stage, check_build_isoforms_for_phase__a_retained_intron_with_a_stop);
    tcase_add_test(tc_stage, check_build_isoforms_for_phase__without_a_trigger_nothing_is_built);
    suite_add_tcase(s, tc_stage);

    return s;
}
