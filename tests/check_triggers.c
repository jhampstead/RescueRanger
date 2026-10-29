#include <check.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>

#include "klib/kvec.h"
#include "record.h"
#include "triggers.h"

#include "check_triggers.h"
#include "mocks.h"

// A shared transcript with two coding exons, [100, 110) and [200, 300), and an intron of 90 bases between them
static Transcript *make_shared(const uint32_t id) {
    Range cds[] = { {100, 109}, {200, 299} };
    int constitutive[] = { 1, 1 };
    Range exons[] = { {100, 109}, {200, 299} };
    TranscriptData data = { .id = id, .start = 100, .end = 299, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 2 };

    return Transcript_build(&data);
}

// A record on chr1 without transcript entries, add them with add_entry. af is a negative value for a record without INFO/AF.
static Record *make_record_with_alleles(bcf_hdr_t *hdr, const int64_t pos, const char *alleles, const int gt0, const int gt1, const bool phased, const float af) {
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, pos, alleles, gt0, gt1, phased, af);

    VepCsqs *vep = malloc(sizeof(VepCsqs));
    kv_init(*vep);

    Record *r = Record_init(v, vep); // Keeps a copy of the record
    bcf_destroy(v);

    r->tscripts = malloc(sizeof(Transcripts));
    kv_init(*r->tscripts);

    return r;
}

static Record *make_record(bcf_hdr_t *hdr, const int64_t pos, const int gt0, const int gt1, const bool phased, const float af) {
    return make_record_with_alleles(hdr, pos, "A,G", gt0, gt1, phased, af);
}

// Adds the VEP consequences for one transcript to a record, together with the record's own copy of the shared transcript.
// A NULL shared transcript is a consequence without a coding transcript. Returns the copy.
static Transcript *add_entry(Record *r, const Transcript *shared, const size_t n, const char *consequences[]) {
    TscriptCsqs *csq = calloc(1, sizeof(TscriptCsqs));
    csq->tid = strdup("transcript");
    csq->num_consequences = n;
    csq->consequences = malloc((n ? n : 1) * sizeof(char *));
    for (size_t i = 0; i < n; i++) csq->consequences[i] = strdup(consequences[i]);
    kv_push(TscriptCsqs *, *r->vep, csq);

    Transcript *copy = shared == NULL ? NULL : Transcript_phase_copy(shared);
    kv_push(Transcript *, *r->tscripts, copy);

    return copy;
}

static void free_record(Record *r) {
    if (r->tscripts != NULL) {
        for (size_t i = 0; i < kv_size(*r->tscripts); i++) {
            if (kv_A(*r->tscripts, i) != NULL) Transcript_phase_destroy(kv_A(*r->tscripts, i));
        }
    }

    Record_destroy(r);
}

static void free_records(Records *rs) {
    for (size_t i = 0; i < kv_size(*rs); i++) free_record(kv_A(*rs, i));
    kv_destroy(*rs);
}

static int csq_trigger(const Record *r, const int entry, const int trigger) {
    return kv_A(*r->vep, entry)->triggers[trigger];
}

/*
* has_af_trigger
*/
START_TEST(check_has_af_trigger__at_most_the_threshold) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *rare = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, 0.001f);
    bcf1_t *common = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, 0.5f);
    bcf1_t *none = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, -1);

    ck_assert_int_eq(has_af_trigger(rare, hdr, "AF", 0.01f), 1);
    ck_assert_int_eq(has_af_trigger(common, hdr, "AF", 0.01f), 0);
    ck_assert_int_eq(has_af_trigger(common, hdr, "AF", 0.5f), 1); // The threshold itself still counts
    ck_assert_int_eq(has_af_trigger(none, hdr, "AF", 0.01f), 0);  // No frequency, no trigger

    bcf_destroy(rare);
    bcf_destroy(common);
    bcf_destroy(none);
    bcf_hdr_destroy(hdr);
}
END_TEST

/*
* set_vep_triggers
*/
START_TEST(check_set_vep_triggers__splice_donor_and_acceptor_set_the_canonical_splice_trigger) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *donor[] = { SPLICE_DONOR_VARIANT };
    const char *acceptor[] = { SPLICE_ACCEPTOR_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1); add_entry(a, shared, 1, donor);
    Record *b = make_record(hdr, 107, 1, 1, true, -1); add_entry(b, shared, 1, acceptor);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_vep_triggers(rs, 1);

    ck_assert_int_eq(csq_trigger(a, 0, SPLICE_TRIGGER), CANONICAL_SPLICE_VALUE);
    ck_assert_int_eq(csq_trigger(b, 0, SPLICE_TRIGGER), CANONICAL_SPLICE_VALUE);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_vep_triggers__frameshift_and_stop_gained) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *frameshift[] = { FRAMESHIFT_VARIANT };
    const char *stop[] = { STOP_GAINED_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1); add_entry(a, shared, 1, frameshift);
    Record *b = make_record(hdr, 107, 1, 1, true, -1); add_entry(b, shared, 1, stop);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_vep_triggers(rs, 1);

    ck_assert_int_eq(csq_trigger(a, 0, FRAMESHIFT_TRIGGER), 1);
    ck_assert_int_eq(csq_trigger(a, 0, STOP_TRIGGER), 0);
    ck_assert_int_eq(csq_trigger(b, 0, STOP_TRIGGER), STOP_RESCUE_VALUE);
    ck_assert_int_eq(csq_trigger(b, 0, FRAMESHIFT_TRIGGER), 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_vep_triggers__other_consequences_set_nothing) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant", SPLICE_REGION_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1); add_entry(a, shared, 2, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_vep_triggers(rs, 1);

    for (int t = 0; t < 3; t++) ck_assert_int_eq(csq_trigger(a, 0, t), 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_vep_triggers__skips_the_non_coding_consequence_but_not_the_others_of_the_entry) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *mixed[] = { "non_coding_transcript_variant", FRAMESHIFT_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1); add_entry(a, shared, 2, mixed);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_vep_triggers(rs, 1);

    ck_assert_int_eq(csq_trigger(a, 0, FRAMESHIFT_TRIGGER), 1);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_vep_triggers__only_the_given_transcript) {
    Transcript *first = make_shared(1);
    Transcript *second = make_shared(2);
    bcf_hdr_t *hdr = mock_hdr();
    const char *frameshift[] = { FRAMESHIFT_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    add_entry(a, first, 1, frameshift);
    add_entry(a, second, 1, frameshift);
    add_entry(a, NULL, 1, frameshift); // No coding transcript for this consequence
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_vep_triggers(rs, 1);

    ck_assert_int_eq(csq_trigger(a, 0, FRAMESHIFT_TRIGGER), 1);
    ck_assert_int_eq(csq_trigger(a, 1, FRAMESHIFT_TRIGGER), 0);
    ck_assert_int_eq(csq_trigger(a, 2, FRAMESHIFT_TRIGGER), 0);

    free_records(&rs);
    Transcript_destroy(first);
    Transcript_destroy(second);
    bcf_hdr_destroy(hdr);
}
END_TEST

/*
* set_af_triggers
*/
START_TEST(check_set_af_triggers__rare_variant_triggers_every_consequence) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, 0.001f);
    add_entry(a, shared, 1, other);
    add_entry(a, NULL, 1, other); // The trigger does not depend on the transcript
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_af_triggers(rs, hdr, "AF", 0.01f);

    ck_assert_int_eq(csq_trigger(a, 0, SPLICE_TRIGGER), AF_SPLICE_VALUE);
    ck_assert_int_eq(csq_trigger(a, 1, SPLICE_TRIGGER), AF_SPLICE_VALUE);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_af_triggers__keeps_a_canonical_splice_trigger) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *donor[] = { SPLICE_DONOR_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, 0.001f);
    add_entry(a, shared, 1, donor);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_vep_triggers(rs, 1);
    set_af_triggers(rs, hdr, "AF", 0.01f);

    ck_assert_int_eq(csq_trigger(a, 0, SPLICE_TRIGGER), CANONICAL_SPLICE_VALUE);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_set_af_triggers__common_or_unknown_frequency_triggers_nothing) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *common = make_record(hdr, 105, 1, 1, true, 0.5f); add_entry(common, shared, 1, other);
    Record *unknown = make_record(hdr, 107, 1, 1, true, -1); add_entry(unknown, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, common); kv_push(Record *, rs, unknown);

    set_af_triggers(rs, hdr, "AF", 0.01f);

    ck_assert_int_eq(csq_trigger(common, 0, SPLICE_TRIGGER), 0);
    ck_assert_int_eq(csq_trigger(unknown, 0, SPLICE_TRIGGER), 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

/*
* set_mnv_triggers: the transcript has coding positions 100..109 (0..9) and, after the 90 base intron, 200..299 (10..109)
*/

// Two records with their own copy of the transcript, at the given positions, and the reading frame trigger of the first record's
// consequence for phase 0 afterwards
static int mnv_trigger_between(const int64_t first_pos, const int64_t second_pos) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, first_pos, 1, 1, true, -1);
    Transcript *copy_a = add_entry(a, shared, 1, other);
    Record *b = make_record(hdr, second_pos, 1, 1, true, -1); add_entry(b, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_mnv_triggers(rs, copy_a, 0);
    int trigger = Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);

    return trigger;
}

// The coding position of a record is its first changed base (variant_position): the first deleted base of a deletion, the base
// an insertion is attached to (the VCF position of a normalised insertion)
static int mnv_trigger_between_alleles(const int64_t first_pos, const char *first_alleles, const int64_t second_pos, const char *second_alleles) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record_with_alleles(hdr, first_pos, first_alleles, 1, 1, true, -1);
    Transcript *copy_a = add_entry(a, shared, 1, other);
    Record *b = make_record_with_alleles(hdr, second_pos, second_alleles, 1, 1, true, -1); add_entry(b, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_mnv_triggers(rs, copy_a, 0);
    int trigger = Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);

    return trigger;
}

START_TEST(check_set_mnv_triggers__an_insertion_counts_at_the_base_it_is_attached_to) {
    ck_assert_int_eq(mnv_trigger_between_alleles(105, "A,AT", 107, "A,G"), STOP_GAIN_VALUE); // Attached to 105
    ck_assert_int_eq(mnv_trigger_between_alleles(104, "A,ATT", 107, "A,G"), 0); // Attached to 104, three bases apart
    ck_assert_int_eq(mnv_trigger_between_alleles(109, "A,AT", 200, "A,G"), STOP_GAIN_VALUE); // Coding positions 9 and 10
} END_TEST

// A deletion counts at its first deleted base, not at the base before it
START_TEST(check_set_mnv_triggers__a_deletion_counts_at_its_first_deleted_base) {
    // Deletes 106; the base before it is 105. 108 is two bases from 106 but three from 105
    ck_assert_int_eq(mnv_trigger_between_alleles(105, "AT,A", 108, "A,G"), STOP_GAIN_VALUE);
    // Deletes 110, the first intron base, so it is not coding although its attachment base 109 is. 200 is coding position 10
    ck_assert_int_eq(mnv_trigger_between_alleles(109, "AT,A", 200, "A,G"), 0);
    // Deletes 200, coding position 10, although its attachment base 199 is in the intron. 201 is coding position 11
    ck_assert_int_eq(mnv_trigger_between_alleles(199, "AT,A", 201, "A,G"), STOP_GAIN_VALUE);
} END_TEST

START_TEST(check_set_mnv_triggers__two_bases_apart_in_an_exon_triggers) {
    ck_assert_int_eq(mnv_trigger_between(105, 107), STOP_GAIN_VALUE);
    ck_assert_int_eq(mnv_trigger_between(105, 106), STOP_GAIN_VALUE);
}
END_TEST

START_TEST(check_set_mnv_triggers__three_bases_apart_does_not_trigger) {
    ck_assert_int_eq(mnv_trigger_between(105, 108), 0);
}
END_TEST

START_TEST(check_set_mnv_triggers__across_an_intron_uses_the_true_intron_length) {
    // Coding positions 8 and 10: 2 apart. The intron is 90 bases, from 110 up to 199
    ck_assert_int_eq(mnv_trigger_between(108, 200), STOP_GAIN_VALUE);

    // Coding positions 7 and 10: 3 apart. An intron length one too large would have made this 2 apart and a trigger
    ck_assert_int_eq(mnv_trigger_between(107, 200), 0);
}
END_TEST

START_TEST(check_set_mnv_triggers__positions_outside_the_coding_exons_are_ignored) {
    ck_assert_int_eq(mnv_trigger_between(150, 151), 0); // In the intron
    ck_assert_int_eq(mnv_trigger_between(50, 51), 0);   // Before the first exon
    ck_assert_int_eq(mnv_trigger_between(400, 401), 0); // After the last exon
}
END_TEST

START_TEST(check_set_mnv_triggers__a_single_coding_record_cannot_trigger) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    Transcript *copy = add_entry(a, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    set_mnv_triggers(rs, copy, 0);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), 0);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

// Both records of a pair get the trigger, on their own consequence for the transcript, for the phase that is given and no other
START_TEST(check_set_mnv_triggers__both_records_of_the_pair_get_it_for_the_phase_only) {
    Transcript *shared = make_shared(1);
    Transcript *other_transcript = make_shared(2);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    Transcript *copy_a = add_entry(a, shared, 1, other);
    add_entry(a, other_transcript, 1, other); // A second consequence, for another transcript
    Record *b = make_record(hdr, 107, 1, 1, true, -1);
    add_entry(b, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_mnv_triggers(rs, copy_a, 1);

    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 1), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 1), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), 0); // The other phase
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 0), 0);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 1), STOP_TRIGGER, 1), 0); // The other transcript

    free_records(&rs);
    Transcript_destroy(shared);
    Transcript_destroy(other_transcript);
    bcf_hdr_destroy(hdr);
}
END_TEST

// A stop that VEP already found on the consequence keeps its value: the reading frame trigger does not replace it
START_TEST(check_set_mnv_triggers__a_stop_trigger_from_vep_keeps_its_value) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };
    const char *stop_gained[] = { STOP_GAINED_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    Transcript *copy_a = add_entry(a, shared, 1, stop_gained);
    Record *b = make_record(hdr, 107, 1, 1, true, -1);
    add_entry(b, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    set_vep_triggers(rs, 1);
    set_mnv_triggers(rs, copy_a, 0);

    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), STOP_RESCUE_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

/*
* Trigger_value: the value of a trigger of a consequence for a phase
*/
START_TEST(check_trigger_value__the_reading_frame_trigger_is_per_phase_and_only_for_the_stop_trigger) {
    TscriptCsqs csq = {0};
    csq.reading_frame[1] = 1;
    csq.triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
    csq.triggers[FRAMESHIFT_TRIGGER] = 1;

    ck_assert_int_eq(Trigger_value(&csq, STOP_TRIGGER, 0), 0);
    ck_assert_int_eq(Trigger_value(&csq, STOP_TRIGGER, 1), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(&csq, SPLICE_TRIGGER, 0), CANONICAL_SPLICE_VALUE);
    ck_assert_int_eq(Trigger_value(&csq, SPLICE_TRIGGER, 1), CANONICAL_SPLICE_VALUE);
    ck_assert_int_eq(Trigger_value(&csq, FRAMESHIFT_TRIGGER, 1), 1);

    csq.triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE; // Found by VEP: the same in every phase, and it comes first
    ck_assert_int_eq(Trigger_value(&csq, STOP_TRIGGER, 0), STOP_RESCUE_VALUE);
    ck_assert_int_eq(Trigger_value(&csq, STOP_TRIGGER, 1), STOP_RESCUE_VALUE);
}
END_TEST

/*
* Triggers_compute_for_phase
*/
START_TEST(check_triggers_compute_for_phase__both_records_of_a_reading_frame_conflict_get_the_trigger) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 0, true, -1); // Both on the first haplotype
    add_entry(a, shared, 1, other);
    Record *b = make_record(hdr, 107, 1, 0, true, -1);
    add_entry(b, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.01f, 0);

    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 1), 0); // Nothing was computed for the other phase

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

// The conflict depends on who is on the haplotype: in the phase where only one of the two is, there is none, and computing that
// phase does not touch the result of the other one
START_TEST(check_triggers_compute_for_phase__the_reading_frame_trigger_follows_the_records_of_the_phase) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, -1); // Homozygous: on both haplotypes
    add_entry(a, shared, 1, other);
    Record *b = make_record(hdr, 107, 1, 0, true, -1); // Only on the first
    add_entry(b, shared, 1, other);

    Records first; kv_init(first);
    kv_push(Record *, first, a); kv_push(Record *, first, b);
    Records second; kv_init(second);
    kv_push(Record *, second, a);

    Triggers_compute_for_phase(&first, hdr, "AF", 0.01f, 0);
    Triggers_compute_for_phase(&second, hdr, "AF", 0.01f, 1);

    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 1), 0);
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
    ck_assert_int_eq(Trigger_value(kv_A(*b->vep, 0), STOP_TRIGGER, 1), 0);

    kv_destroy(first);
    kv_destroy(second);
    free_record(a);
    free_record(b);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

// Computing a phase again with other records starts again: a conflict that is gone is not kept
START_TEST(check_triggers_compute_for_phase__computing_a_phase_again_clears_its_reading_frame_trigger) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    add_entry(a, shared, 1, other);
    Record *b = make_record(hdr, 107, 1, 1, true, -1);
    add_entry(b, shared, 1, other);
    Records both; kv_init(both);
    kv_push(Record *, both, a); kv_push(Record *, both, b);
    Records only_a; kv_init(only_a);
    kv_push(Record *, only_a, a);

    Triggers_compute_for_phase(&both, hdr, "AF", 0.01f, 0);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), STOP_GAIN_VALUE);

    Triggers_compute_for_phase(&only_a, hdr, "AF", 0.01f, 0);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), 0);

    kv_destroy(both);
    kv_destroy(only_a);
    free_record(a);
    free_record(b);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_triggers_compute_for_phase__only_the_given_records_get_triggers) {
    // The records are the ones selected for the phase, so a record that is not in the list gets nothing and does not
    // take part in a reading frame conflict, like a record that was removed because it conflicts with an earlier variant
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *frameshift[] = { FRAMESHIFT_VARIANT };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    add_entry(a, shared, 1, frameshift);
    Record *removed = make_record(hdr, 106, 1, 1, true, 0.001f);
    add_entry(removed, shared, 1, frameshift);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.01f, 0);

    ck_assert_int_eq(csq_trigger(a, 0, FRAMESHIFT_TRIGGER), 1);
    ck_assert_int_eq(csq_trigger(removed, 0, FRAMESHIFT_TRIGGER), 0);
    ck_assert_int_eq(csq_trigger(removed, 0, SPLICE_TRIGGER), 0);
    ck_assert_int_eq(Trigger_value(kv_A(*a->vep, 0), STOP_TRIGGER, 0), 0);       // One record on its own cannot conflict, and 106 is not in the list
    ck_assert_int_eq(Trigger_value(kv_A(*removed->vep, 0), STOP_TRIGGER, 0), 0);

    free_records(&rs);
    free_record(removed);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_triggers_compute_for_phase__allele_frequency_threshold_is_passed_on) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, 0.001f);
    add_entry(a, shared, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.0001f, 0);
    ck_assert_int_eq(csq_trigger(a, 0, SPLICE_TRIGGER), 0);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.01f, 0);
    ck_assert_int_eq(csq_trigger(a, 0, SPLICE_TRIGGER), AF_SPLICE_VALUE);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_triggers_compute_for_phase__records_without_transcripts_are_skipped) {
    Transcript *shared = make_shared(1);
    bcf_hdr_t *hdr = mock_hdr();
    const char *frameshift[] = { FRAMESHIFT_VARIANT };

    // The loader leaves tscripts NULL for a record without coding transcripts
    Record *without = make_record(hdr, 105, 1, 1, true, 0.001f);
    kv_destroy(*without->tscripts);
    free(without->tscripts);
    without->tscripts = NULL;

    Record *with = make_record(hdr, 107, 1, 1, true, -1);
    add_entry(with, shared, 1, frameshift);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, without); kv_push(Record *, rs, with);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.01f, 0);

    ck_assert_int_eq(csq_trigger(with, 0, FRAMESHIFT_TRIGGER), 1);

    free_records(&rs);
    Transcript_destroy(shared);
    bcf_hdr_destroy(hdr);
}
END_TEST

START_TEST(check_triggers_compute_for_phase__every_transcript_gets_its_own_triggers) {
    Transcript *first = make_shared(1);
    Transcript *second = make_shared(2);
    bcf_hdr_t *hdr = mock_hdr();
    const char *other[] = { "missense_variant" };

    Record *a = make_record(hdr, 105, 1, 1, true, -1);
    add_entry(a, first, 1, other);
    add_entry(a, second, 1, other);
    Record *b = make_record(hdr, 107, 1, 1, true, -1);
    add_entry(b, first, 1, other);
    add_entry(b, second, 1, other);
    Records rs; kv_init(rs);
    kv_push(Record *, rs, a); kv_push(Record *, rs, b);

    Triggers_compute_for_phase(&rs, hdr, "AF", 0.01f, 0);

    for (int entry = 0; entry < 2; entry++) {
        ck_assert_int_eq(Trigger_value(kv_A(*a->vep, entry), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
        ck_assert_int_eq(Trigger_value(kv_A(*b->vep, entry), STOP_TRIGGER, 0), STOP_GAIN_VALUE);
    }

    free_records(&rs);
    Transcript_destroy(first);
    Transcript_destroy(second);
    bcf_hdr_destroy(hdr);
}
END_TEST

Suite *triggers_suite() {
    Suite *s = suite_create("Triggers");

    TCase *tc_af = tcase_create("has_af_trigger");
    suite_add_tcase(s, tc_af);
    tcase_add_test(tc_af, check_has_af_trigger__at_most_the_threshold);

    TCase *tc_vep = tcase_create("set_vep_triggers");
    suite_add_tcase(s, tc_vep);
    tcase_add_test(tc_vep, check_set_vep_triggers__splice_donor_and_acceptor_set_the_canonical_splice_trigger);
    tcase_add_test(tc_vep, check_set_vep_triggers__frameshift_and_stop_gained);
    tcase_add_test(tc_vep, check_set_vep_triggers__other_consequences_set_nothing);
    tcase_add_test(tc_vep, check_set_vep_triggers__skips_the_non_coding_consequence_but_not_the_others_of_the_entry);
    tcase_add_test(tc_vep, check_set_vep_triggers__only_the_given_transcript);

    TCase *tc_af_triggers = tcase_create("set_af_triggers");
    suite_add_tcase(s, tc_af_triggers);
    tcase_add_test(tc_af_triggers, check_set_af_triggers__rare_variant_triggers_every_consequence);
    tcase_add_test(tc_af_triggers, check_set_af_triggers__keeps_a_canonical_splice_trigger);
    tcase_add_test(tc_af_triggers, check_set_af_triggers__common_or_unknown_frequency_triggers_nothing);

    TCase *tc_mnv = tcase_create("set_mnv_triggers");
    suite_add_tcase(s, tc_mnv);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__two_bases_apart_in_an_exon_triggers);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__an_insertion_counts_at_the_base_it_is_attached_to);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__a_deletion_counts_at_its_first_deleted_base);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__three_bases_apart_does_not_trigger);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__across_an_intron_uses_the_true_intron_length);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__positions_outside_the_coding_exons_are_ignored);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__a_single_coding_record_cannot_trigger);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__both_records_of_the_pair_get_it_for_the_phase_only);
    tcase_add_test(tc_mnv, check_set_mnv_triggers__a_stop_trigger_from_vep_keeps_its_value);
    tcase_add_test(tc_mnv, check_trigger_value__the_reading_frame_trigger_is_per_phase_and_only_for_the_stop_trigger);

    TCase *tc_phase = tcase_create("compute_for_phase");
    suite_add_tcase(s, tc_phase);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__both_records_of_a_reading_frame_conflict_get_the_trigger);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__the_reading_frame_trigger_follows_the_records_of_the_phase);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__computing_a_phase_again_clears_its_reading_frame_trigger);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__only_the_given_records_get_triggers);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__allele_frequency_threshold_is_passed_on);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__records_without_transcripts_are_skipped);
    tcase_add_test(tc_phase, check_triggers_compute_for_phase__every_transcript_gets_its_own_triggers);

    return s;
}
