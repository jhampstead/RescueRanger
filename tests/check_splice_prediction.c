#include <check.h>
#include <stdlib.h>
#include <cpliceai/range.h>

#include "bcftools/gff.h"
#include "klib/kvec.h"
#include "mocks.h"
#include "record.h"
#include "splice_prediction.h"
#include "utils.h"

#include "check_splice_prediction.h"

/*
* Small helpers on the predictions: they work on plain arrays of NUM_SCORES values per position.
*/

// A session without precomputed reference scores: the caller has to predict the sequence, so there are none to give
START_TEST(check_precomputed_reference_predictions__without_reference_scores_there_are_none) {
    Session session = {0};
    session.prediction_radius = 50;
    Gene gene = { .name = "G1", .coords = { 100, 160 }, .strand = STRAND_REV };

    ck_assert_ptr_null(precomputed_reference_predictions(&session, &gene, "chr1", 110));
} END_TEST

START_TEST(check_filter_to_most_likely_prediction__keeps_the_largest_of_every_position) {
    float probabilities[] = { 0.1f, 0.5f, 0.2f,   0.7f, 0.2f, 0.1f,   0.0f, 0.1f, 0.9f };

    filter_to_most_likely_prediction(9, probabilities);

    const float expected[] = { 0.0f, 0.5f, 0.0f,   0.7f, 0.0f, 0.0f,   0.0f, 0.0f, 0.9f };
    for (int i = 0; i < 9; i++) ck_assert_float_eq(probabilities[i], expected[i]);
} END_TEST

// Of equal values the first one stays; a position without any probability stays 0
START_TEST(check_filter_to_most_likely_prediction__ties_and_zeros) {
    float probabilities[] = { 0.4f, 0.4f, 0.1f,   0.0f, 0.0f, 0.0f };

    filter_to_most_likely_prediction(6, probabilities);

    ck_assert_float_eq(probabilities[0], 0.4f);
    ck_assert_float_eq(probabilities[1], 0.0f);
    ck_assert_float_eq(probabilities[2], 0.0f);
    for (int i = 3; i < 6; i++) ck_assert_float_eq(probabilities[i], 0.0f);
} END_TEST

START_TEST(check_transform_to_posterior_probabilities__scales_the_nonzero_values_to_sum_to_one) {
    const float probabilities[] = { 0.0f, 0.25f, 0.75f, 0.0f };
    float posteriors[4] = { 9, 9, 9, 9 };

    const int n = transform_to_posterior_probabilities(4, probabilities, posteriors);

    ck_assert_int_eq(n, 2);
    ck_assert_float_eq(posteriors[0], 0.0f);
    ck_assert_float_eq_tol(posteriors[1], 0.25f, 1e-6);
    ck_assert_float_eq_tol(posteriors[2], 0.75f, 1e-6);
    ck_assert_float_eq(posteriors[3], 0.0f);
} END_TEST

START_TEST(check_transform_to_posterior_probabilities__unequal_values) {
    const float probabilities[] = { 0.2f, 0.6f };
    float posteriors[2];

    ck_assert_int_eq(transform_to_posterior_probabilities(2, probabilities, posteriors), 2);
    ck_assert_float_eq_tol(posteriors[0], 0.25f, 1e-6);
    ck_assert_float_eq_tol(posteriors[1], 0.75f, 1e-6);
} END_TEST

START_TEST(check_transform_to_posterior_probabilities__all_zero) {
    const float probabilities[] = { 0.0f, 0.0f, 0.0f };
    float posteriors[3] = { 9, 9, 9 };

    ck_assert_int_eq(transform_to_posterior_probabilities(3, probabilities, posteriors), 0);
    for (int i = 0; i < 3; i++) ck_assert_float_eq(posteriors[i], 0.0f);
} END_TEST

START_TEST(check_check_proximity__within_the_threshold_on_both_sides) {
    ck_assert(check_proximity(100, 100));
    ck_assert(check_proximity(115, 100));
    ck_assert(check_proximity(85, 100));
    ck_assert(!check_proximity(116, 100));
    ck_assert(!check_proximity(84, 100));
} END_TEST

START_TEST(check_kv_push_new_record__a_record_is_added_once) {
    Records rs; kv_init(rs);
    Record *a = calloc(1, sizeof(Record));
    Record *b = calloc(1, sizeof(Record));

    kv_push_new_record(&rs, a);
    kv_push_new_record(&rs, b);
    kv_push_new_record(&rs, a);
    kv_push_new_record(&rs, b);

    ck_assert_int_eq(kv_size(rs), 2);
    ck_assert_ptr_eq(kv_A(rs, 0), a);
    ck_assert_ptr_eq(kv_A(rs, 1), b);

    free(a); free(b);
    kv_destroy(rs);
} END_TEST

/*
* store_splice_sites: what a predicted splice site does to the sites of an exon. The sites of the exon: a canonical donor at 200
* and a canonical acceptor at 100, both from the GFF.
*/

static SpliceSites *sites_with_canonicals(void) {
    SpliceSites *sites = malloc(sizeof(SpliceSites));
    init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
    init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
    init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
    init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

    Records *donor_records = malloc(sizeof(Records)); kv_init(*donor_records);
    sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { donor_records, 200, 200, 0, 1, 0, 1.0, 1.0, 1.0 };

    Records *acceptor_records = malloc(sizeof(Records)); kv_init(*acceptor_records);
    sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { acceptor_records, 100, 100, 0, 1, 0, 1.0, 1.0, 1.0 };

    return sites;
}

static SpliceProbability prediction_of(Record *r, const enum SpliceType type, const float ref_prior, const float prior, const float posterior, const uint32_t position, const uint32_t canonical) {
    return (SpliceProbability) { r, type, ref_prior, prior, posterior, position, canonical };
}

START_TEST(check_store_splice_sites__a_lost_canonical_donor) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, DONOR, 0.9f, 0.0f, 0.0f, 200, 200), &pair, &competitor, sites);

    const SpliceSite *donor = &sites->canonical_donors[0];
    ck_assert_int_eq(donor->predicted, 1);
    ck_assert_int_eq(donor->lost, 1);
    ck_assert_float_eq(donor->ref_prior_probability, 0.9f);
    ck_assert_float_eq(donor->prior_probability, 0.0f);
    ck_assert_float_eq(donor->posterior_probability, 0.0f);
    ck_assert_int_eq(kv_size(*(Records *) donor->records_ptr), 1);
    ck_assert_ptr_eq(kv_A(*(Records *) donor->records_ptr, 0), r);
    ck_assert_int_eq(pair, 1);
    ck_assert_int_eq(competitor, 1);
    ck_assert_int_eq(sites->n_canonical_donors, 1);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

START_TEST(check_store_splice_sites__a_lost_canonical_acceptor) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, ACCEPTOR, 0.8f, 0.0f, 0.0f, 100, 100), &pair, &competitor, sites);

    const SpliceSite *acceptor = &sites->canonical_acceptors[0];
    ck_assert_int_eq(acceptor->predicted, 1);
    ck_assert_int_eq(acceptor->lost, 1);
    ck_assert_float_eq(acceptor->ref_prior_probability, 0.8f);
    ck_assert_int_eq(kv_size(*(Records *) acceptor->records_ptr), 1);
    ck_assert_int_eq(pair, 1);
    ck_assert_int_eq(competitor, 1);
    ck_assert_int_eq(sites->canonical_donors[0].predicted, 0); // The other kind is left alone

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

// The canonical site is the only one in the window: it is kept as it is predicted, and not lost
START_TEST(check_store_splice_sites__the_only_site_in_the_window_is_the_canonical_one) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, DONOR, 0.9f, 0.85f, 1.0f, 200, 200), &pair, &competitor, sites);

    const SpliceSite *donor = &sites->canonical_donors[0];
    ck_assert_int_eq(donor->predicted, 1);
    ck_assert_int_eq(donor->lost, 0);
    ck_assert_float_eq(donor->prior_probability, 0.85f);
    ck_assert_float_eq(donor->posterior_probability, 1.0f);
    ck_assert_int_eq(pair, 1);
    ck_assert_int_eq(competitor, 1);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

// A site close to the canonical one that is predicted too is a canonical site of its own, and the pair is predicted for it
START_TEST(check_store_splice_sites__a_site_next_to_the_canonical_one_becomes_a_canonical_site) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, DONOR, 0.0f, 0.6f, 0.6f, 205, 200), &pair, &competitor, sites);

    ck_assert_int_eq(sites->n_canonical_donors, 2);
    const SpliceSite *added = &sites->canonical_donors[1];
    ck_assert_int_eq(added->position, 205);
    ck_assert_int_eq(added->canonical_position, 200);
    ck_assert_int_eq(added->predicted, 1);
    ck_assert_int_eq(added->gff_canonical, 0);
    ck_assert_float_eq(added->posterior_probability, 0.6f);
    ck_assert_int_eq(kv_size(*(Records *) added->records_ptr), 1);
    ck_assert_int_eq(pair, 1);
    ck_assert_int_eq(competitor, 0);
    ck_assert_int_eq(sites->canonical_donors[0].predicted, 0); // The GFF site itself was not predicted

    // The same site again changes nothing and asks for no new prediction
    pair = 0;
    store_splice_sites(prediction_of(r, DONOR, 0.0f, 0.6f, 0.6f, 205, 200), &pair, &competitor, sites);
    ck_assert_int_eq(sites->n_canonical_donors, 2);
    ck_assert_int_eq(pair, 0);
    ck_assert_int_eq(kv_size(*(Records *) sites->canonical_donors[1].records_ptr), 1);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

// Another record that predicts the same site is added to the site, which is not predicted again
START_TEST(check_store_splice_sites__another_record_joins_a_known_site) {
    SpliceSites *sites = sites_with_canonicals();
    Record *first = calloc(1, sizeof(Record));
    Record *second = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(first, ACCEPTOR, 0.0f, 0.6f, 0.6f, 95, 100), &pair, &competitor, sites);
    pair = 0;
    store_splice_sites(prediction_of(second, ACCEPTOR, 0.0f, 0.6f, 0.6f, 95, 100), &pair, &competitor, sites);

    ck_assert_int_eq(sites->n_canonical_acceptors, 2);
    ck_assert_int_eq(pair, 0);
    Records *records = (Records *) sites->canonical_acceptors[1].records_ptr;
    ck_assert_int_eq(kv_size(*records), 2);
    ck_assert_ptr_eq(kv_A(*records, 0), first);
    ck_assert_ptr_eq(kv_A(*records, 1), second);

    SpliceSites_destroy(sites);
    free(first); free(second);
} END_TEST

// A site far from the canonical one is a cryptic site, and a competitor is predicted for it
START_TEST(check_store_splice_sites__a_far_site_is_cryptic) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, DONOR, 0.0f, 0.5f, 0.5f, 260, 200), &pair, &competitor, sites);

    ck_assert_int_eq(sites->n_cryptic_donors, 1);
    ck_assert_int_eq(sites->cryptic_donors[0].position, 260);
    ck_assert_int_eq(sites->cryptic_donors[0].canonical_position, 200);
    ck_assert_int_eq(sites->cryptic_donors[0].predicted, 1);
    ck_assert_int_eq(sites->n_canonical_donors, 1);
    ck_assert_int_eq(pair, 0);
    ck_assert_int_eq(competitor, 1);

    // Again: known, so no competitor is asked for
    competitor = 0;
    store_splice_sites(prediction_of(r, DONOR, 0.0f, 0.5f, 0.5f, 260, 200), &pair, &competitor, sites);
    ck_assert_int_eq(sites->n_cryptic_donors, 1);
    ck_assert_int_eq(competitor, 0);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

START_TEST(check_store_splice_sites__many_cryptic_sites_grow_the_array) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));

    for (int i = 0; i < 6; i++) {
        int pair = 0, competitor = 0;
        store_splice_sites(prediction_of(r, ACCEPTOR, 0.0f, 0.3f, 0.3f, 20 + i, 100), &pair, &competitor, sites);
    }

    ck_assert_int_eq(sites->n_cryptic_acceptors, 6);
    for (int i = 0; i < 6; i++) ck_assert_int_eq(sites->cryptic_acceptors[i].position, 20 + i);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

// A site with no probability left that is not the canonical one is nothing to store
START_TEST(check_store_splice_sites__an_unused_other_site_is_ignored) {
    SpliceSites *sites = sites_with_canonicals();
    Record *r = calloc(1, sizeof(Record));
    int pair = 0, competitor = 0;

    store_splice_sites(prediction_of(r, DONOR, 0.5f, 0.0f, 0.0f, 230, 200), &pair, &competitor, sites);

    ck_assert_int_eq(sites->n_canonical_donors, 1);
    ck_assert_int_eq(sites->n_cryptic_donors, 0);
    ck_assert_int_eq(sites->canonical_donors[0].predicted, 0);
    ck_assert_int_eq(pair, 0);
    ck_assert_int_eq(competitor, 0);

    SpliceSites_destroy(sites);
    free(r);
} END_TEST

/*
* init_exon_splice_sites and init_splice_sites_for_phase: the canonical sites an exon has in the GFF. Three exons [100, 150),
* [200, 250) and [300, 350).
*/

static Transcript *three_exon_transcript(const uint32_t strand) {
    Range cds[] = { {100, 149}, {200, 249}, {300, 349} };
    int constitutive[] = { 1, 1, 1 };
    Range exons[] = { {100, 149}, {200, 249}, {300, 349} };
    TranscriptData data = { .id = 1, .start = 100, .end = 349, .strand = strand,
                            .cds = cds, .constitutive = constitutive, .ncds = 3, .exons = exons, .nexons = 3 };

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

START_TEST(check_init_exon_splice_sites__forward_strand) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    // The first exon has no acceptor and the last no donor
    const size_t acceptors[] = { 0, 1, 1 };
    const size_t donors[] = { 1, 1, 0 };
    for (size_t ei = 0; ei < 3; ei++) {
        const Exon *e = kv_A(t->exons, ei);
        ck_assert_int_eq(e->splice_sites->n_canonical_acceptors, acceptors[ei]);
        ck_assert_int_eq(e->splice_sites->n_canonical_donors, donors[ei]);
        ck_assert_int_eq(e->splice_sites->n_cryptic_acceptors, 0);
        ck_assert_int_eq(e->splice_sites->n_cryptic_donors, 0);

        if (donors[ei]) {
            const SpliceSite *site = &e->splice_sites->canonical_donors[0];
            ck_assert_int_eq(site->position, Exon_donor_base(e, STRAND_FWD));
            ck_assert_int_eq(site->canonical_position, site->position);
            ck_assert_int_eq(site->predicted, 0);
            ck_assert_int_eq(site->gff_canonical, 1);
            ck_assert_int_eq(site->lost, 0);
            ck_assert_float_eq(site->posterior_probability, 1.0f);
            ck_assert_int_eq(kv_size(*(Records *) site->records_ptr), 0);
        }
        if (acceptors[ei]) {
            const SpliceSite *site = &e->splice_sites->canonical_acceptors[0];
            ck_assert_int_eq(site->position, Exon_acceptor_base(e, STRAND_FWD));
            ck_assert_int_eq(site->gff_canonical, 1);
        }
    }

    ck_assert_int_eq(kv_A(t->exons, 0)->splice_sites->canonical_donors[0].position, 149);
    ck_assert_int_eq(kv_A(t->exons, 1)->splice_sites->canonical_acceptors[0].position, 200);

    destroy_transcript_with_sites(t);
} END_TEST

// On the reverse strand the first exon in transcription order is the highest one
START_TEST(check_init_exon_splice_sites__reverse_strand) {
    Transcript *t = three_exon_transcript(STRAND_REV);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    const size_t acceptors[] = { 1, 1, 0 };
    const size_t donors[] = { 0, 1, 1 };
    for (size_t ei = 0; ei < 3; ei++) {
        const Exon *e = kv_A(t->exons, ei);
        ck_assert_int_eq(e->splice_sites->n_canonical_acceptors, acceptors[ei]);
        ck_assert_int_eq(e->splice_sites->n_canonical_donors, donors[ei]);

        if (donors[ei]) ck_assert_int_eq(e->splice_sites->canonical_donors[0].position, Exon_donor_base(e, STRAND_REV));
        if (acceptors[ei]) ck_assert_int_eq(e->splice_sites->canonical_acceptors[0].position, Exon_acceptor_base(e, STRAND_REV));
    }

    ck_assert_int_eq(kv_A(t->exons, 1)->splice_sites->canonical_donors[0].position, 200);
    ck_assert_int_eq(kv_A(t->exons, 0)->splice_sites->canonical_acceptors[0].position, 149);

    destroy_transcript_with_sites(t);
} END_TEST

START_TEST(check_init_exon_splice_sites__unknown_strand_has_no_sites) {
    Transcript *t = three_exon_transcript(STRAND_UNK);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    for (size_t ei = 0; ei < 3; ei++) {
        const Exon *e = kv_A(t->exons, ei);
        ck_assert_ptr_nonnull(e->splice_sites);
        ck_assert_int_eq(e->splice_sites->n_canonical_acceptors, 0);
        ck_assert_int_eq(e->splice_sites->n_canonical_donors, 0);
    }

    destroy_transcript_with_sites(t);
} END_TEST

// A record with one transcript of a gene [100, 350)
static Record *record_on_transcript(const int64_t pos, Transcript *t, Gene *gene) {
    Record *r = Record_init(mock_bcf1_t(pos, "N", "x"), NULL);

    r->tscripts = malloc(sizeof(Transcripts)); kv_init(*r->tscripts);
    kv_push(Transcript *, *r->tscripts, t);
    r->transcript_genes = malloc(sizeof(Genes)); kv_init(*r->transcript_genes);
    kv_push(Gene *, *r->transcript_genes, gene);

    return r;
}

static void destroy_record_on_transcript(Record *r) {
    kv_destroy(*r->tscripts); free(r->tscripts); r->tscripts = NULL;
    kv_destroy(*r->transcript_genes); free(r->transcript_genes); r->transcript_genes = NULL;
    Record_destroy(r);
}

START_TEST(check_init_splice_sites_for_phase__a_record_in_the_gene_gives_its_exons_sites) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 350 } };
    Record *r = record_on_transcript(120, t, &gene);
    Records records; kv_init(records);
    kv_push(Record *, records, r);

    init_splice_sites_for_phase(&records);

    for (size_t ei = 0; ei < 3; ei++) ck_assert_ptr_nonnull(kv_A(t->exons, ei)->splice_sites);
    ck_assert_int_eq(kv_A(t->exons, 1)->splice_sites->n_canonical_donors, 1);

    kv_destroy(records);
    destroy_record_on_transcript(r);
    destroy_transcript_with_sites(t);
} END_TEST

START_TEST(check_init_splice_sites_for_phase__a_record_outside_the_gene_or_without_transcript_gives_none) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 350 } };
    Record *outside = record_on_transcript(500, t, &gene);

    Record *without = record_on_transcript(120, NULL, &gene); // No transcript for this entry
    Records records; kv_init(records);
    kv_push(Record *, records, outside);
    kv_push(Record *, records, without);

    init_splice_sites_for_phase(&records);

    for (size_t ei = 0; ei < 3; ei++) ck_assert_ptr_null(kv_A(t->exons, ei)->splice_sites);

    kv_destroy(records);
    destroy_record_on_transcript(outside);
    destroy_record_on_transcript(without);
    destroy_transcript_with_sites(t);
} END_TEST

/*
* locate_prediction_on_exon: whether a predicted site is in the splice window of an exon. If it is, its canonical position becomes
* the position of the canonical site of that exon.
*/

START_TEST(check_locate_prediction_on_exon__a_donor_in_the_window) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    Range window;
    ck_assert(Exon_splice_window(t, 0, true, &window));

    SpliceProbability inside = prediction_of(NULL, DONOR, 0, 0.5f, 0.5f, window.start, 0);
    ck_assert(locate_prediction_on_exon(t, 0, &inside));
    ck_assert_int_eq(inside.canonical_position, 149);

    SpliceProbability last = prediction_of(NULL, DONOR, 0, 0.5f, 0.5f, window.end - 1, 0);
    ck_assert(locate_prediction_on_exon(t, 0, &last));

    SpliceProbability after = prediction_of(NULL, DONOR, 0, 0.5f, 0.5f, window.end, 7);
    ck_assert(!locate_prediction_on_exon(t, 0, &after));
    ck_assert_int_eq(after.canonical_position, 7); // Not touched

    SpliceProbability before = prediction_of(NULL, DONOR, 0, 0.5f, 0.5f, window.start - 1, 7);
    ck_assert(!locate_prediction_on_exon(t, 0, &before));

    destroy_transcript_with_sites(t);
} END_TEST

START_TEST(check_locate_prediction_on_exon__an_acceptor_in_the_window) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    Range window;
    ck_assert(Exon_splice_window(t, 1, false, &window));

    SpliceProbability inside = prediction_of(NULL, ACCEPTOR, 0, 0.5f, 0.5f, window.start, 0);
    ck_assert(locate_prediction_on_exon(t, 1, &inside));
    ck_assert_int_eq(inside.canonical_position, 200);

    SpliceProbability outside = prediction_of(NULL, ACCEPTOR, 0, 0.5f, 0.5f, window.end, 0);
    ck_assert(!locate_prediction_on_exon(t, 1, &outside));

    destroy_transcript_with_sites(t);
} END_TEST

// The first exon has no acceptor window and the last has no donor window
START_TEST(check_locate_prediction_on_exon__an_exon_without_a_window_gives_false) {
    Transcript *t = three_exon_transcript(STRAND_FWD);
    for (size_t ei = 0; ei < 3; ei++) init_exon_splice_sites(kv_A(t->exons, ei), ei, t);

    SpliceProbability acceptor = prediction_of(NULL, ACCEPTOR, 0, 0.5f, 0.5f, 100, 0);
    ck_assert(!locate_prediction_on_exon(t, 0, &acceptor));

    SpliceProbability donor = prediction_of(NULL, DONOR, 0, 0.5f, 0.5f, 349, 0);
    ck_assert(!locate_prediction_on_exon(t, 2, &donor));

    destroy_transcript_with_sites(t);
} END_TEST

/*
* collect_splice_probabilities: the splice sites a record changes in a window, from the predictions of the reference and of the
* haplotype. The window here is [100, 111) around a record at 105 (radius 5), and the first exon of the transcript ends at the
* base 105, so its donor is the base 105. The predictions are plain arrays with NUM_SCORES values per position.
*/

static Transcript *donor_at_105_transcript(void) {
    Range cds[] = { {100, 105}, {200, 250} };
    int constitutive[] = { 1, 1 };
    Range exons[] = { {100, 105}, {200, 250} };
    TranscriptData data = { .id = 1, .start = 100, .end = 250, .strand = STRAND_FWD,
                            .cds = cds, .constitutive = constitutive, .ncds = 2, .exons = exons, .nexons = 2 };

    return Transcript_build(&data);
}

#define WINDOW_RADIUS 5
#define WINDOW_POSITIONS 11

START_TEST(check_collect_splice_probabilities__a_lost_canonical_donor) {
    Transcript *t = donor_at_105_transcript();
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 251 } };
    Record *r = record_on_transcript(105, t, &gene);
    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, r);

    float ref[WINDOW_POSITIONS * NUM_SCORES] = {0};
    float hap[WINDOW_POSITIONS * NUM_SCORES] = {0};
    ref[(105 - 100) * NUM_SCORES + DONOR_POS] = 0.9f; // The donor is predicted on the reference and gone on the haplotype

    int in_gene[1] = { 1 };
    SpliceProbabilities probabilities; kv_init(probabilities);
    collect_splice_probabilities(r, 0, in_gene, range_around(105, WINDOW_RADIUS), WINDOW_RADIUS, ref, hap, 0, &hap_in, &probabilities);

    ck_assert_int_eq(kv_size(probabilities), 1);
    const SpliceProbability p = kv_A(probabilities, 0);
    ck_assert_ptr_eq(p.record, r);
    ck_assert_int_eq(p.type, DONOR);
    ck_assert_float_eq(p.ref_prior_probability, 0.9f);
    ck_assert_float_eq(p.prior_probability, 0.0f);
    ck_assert_float_eq(p.posterior_probability, 0.0f);
    ck_assert_int_eq(p.position, 105);
    ck_assert_int_eq(p.canonical_position, 105);
    ck_assert_uint_eq(r->flags[0], ABLATED_CANONICAL_DONOR);
    ck_assert_uint_eq(r->flags[1], 0); // The other phase

    kv_destroy(probabilities);
    kv_destroy(hap_in);
    destroy_record_on_transcript(r);
    Transcript_destroy(t);
} END_TEST

// A donor that is only on the haplotype, away from the canonical one, is cryptic; the canonical donor is weak on the reference
START_TEST(check_collect_splice_probabilities__a_cryptic_donor_on_the_haplotype) {
    Transcript *t = donor_at_105_transcript();
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 251 } };
    Record *r = record_on_transcript(105, t, &gene);
    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, r);

    float ref[WINDOW_POSITIONS * NUM_SCORES] = {0};
    float hap[WINDOW_POSITIONS * NUM_SCORES] = {0};
    hap[(103 - 100) * NUM_SCORES + DONOR_POS] = 0.8f;

    int in_gene[1] = { 1 };
    SpliceProbabilities probabilities; kv_init(probabilities);
    collect_splice_probabilities(r, 0, in_gene, range_around(105, WINDOW_RADIUS), WINDOW_RADIUS, ref, hap, 0, &hap_in, &probabilities);

    ck_assert_int_eq(kv_size(probabilities), 1);
    const SpliceProbability p = kv_A(probabilities, 0);
    ck_assert_int_eq(p.type, DONOR);
    ck_assert_float_eq(p.ref_prior_probability, 0.0f);
    ck_assert_float_eq(p.prior_probability, 0.8f);
    ck_assert_float_eq(p.posterior_probability, 1.0f); // The only site of the haplotype
    ck_assert_int_eq(p.position, 103);
    ck_assert_int_eq(p.canonical_position, 0); // Not known yet: it is set when the site is placed on an exon
    ck_assert_uint_eq(r->flags[0], WEAK_CANONICAL_DONOR | CRYPTIC_DONOR);

    kv_destroy(probabilities);
    kv_destroy(hap_in);
    destroy_record_on_transcript(r);
    Transcript_destroy(t);
} END_TEST

// A canonical donor that is weak on the reference and viable on the haplotype is collected with every field in its place
START_TEST(check_collect_splice_probabilities__a_weak_canonical_donor_viable_on_the_haplotype) {
    Transcript *t = donor_at_105_transcript();
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 251 } };
    Record *r = record_on_transcript(105, t, &gene);
    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, r);

    float ref[WINDOW_POSITIONS * NUM_SCORES] = {0};
    float hap[WINDOW_POSITIONS * NUM_SCORES] = {0};
    ref[(105 - 100) * NUM_SCORES + DONOR_POS] = 0.005f;
    hap[(105 - 100) * NUM_SCORES + DONOR_POS] = 0.8f;

    int in_gene[1] = { 1 };
    SpliceProbabilities probabilities; kv_init(probabilities);
    collect_splice_probabilities(r, 0, in_gene, range_around(105, WINDOW_RADIUS), WINDOW_RADIUS, ref, hap, 0, &hap_in, &probabilities);

    ck_assert_int_eq(kv_size(probabilities), 1);
    const SpliceProbability p = kv_A(probabilities, 0);
    ck_assert_ptr_eq(p.record, r);
    ck_assert_int_eq(p.type, DONOR);
    ck_assert_float_eq(p.ref_prior_probability, 0.005f);
    ck_assert_float_eq(p.prior_probability, 0.8f);
    ck_assert_float_eq(p.posterior_probability, 1.0f); // The only site of the haplotype
    ck_assert_int_eq(p.position, 105);
    ck_assert_int_eq(p.canonical_position, 105);
    ck_assert_uint_eq(r->flags[0], WEAK_CANONICAL_DONOR);

    kv_destroy(probabilities);
    kv_destroy(hap_in);
    destroy_record_on_transcript(r);
    Transcript_destroy(t);
} END_TEST

// Nothing changes: nothing is collected and nothing is flagged
START_TEST(check_collect_splice_probabilities__no_predictions_gives_nothing) {
    Transcript *t = donor_at_105_transcript();
    Gene gene = { .gff_id = 1, .name = "G", .coords = { 100, 251 } };
    Record *r = record_on_transcript(105, t, &gene);
    r->flags = calloc(NUM_PHASES, sizeof(uint32_t));
    Records hap_in; kv_init(hap_in);
    kv_push(Record *, hap_in, r);

    float ref[WINDOW_POSITIONS * NUM_SCORES] = {0};
    float hap[WINDOW_POSITIONS * NUM_SCORES] = {0};

    int in_gene[1] = { 1 };
    SpliceProbabilities probabilities; kv_init(probabilities);
    collect_splice_probabilities(r, 0, in_gene, range_around(105, WINDOW_RADIUS), WINDOW_RADIUS, ref, hap, 0, &hap_in, &probabilities);

    ck_assert_int_eq(kv_size(probabilities), 0);
    // The canonical donor is in the window and has no probability on the reference
    ck_assert_uint_eq(r->flags[0], WEAK_CANONICAL_DONOR);

    kv_destroy(probabilities);
    kv_destroy(hap_in);
    destroy_record_on_transcript(r);
    Transcript_destroy(t);
} END_TEST

/*
* Where the sequence of the gene lies in the model window. The window has to hold exactly the bases of the sequence that are in
* it: a window that takes more reads past the end of the sequence, and one that takes fewer puts N over bases it has.
*/
#define SEQUENCE_WINDOW_RADIUS 50
#define SEQUENCE_WINDOW_GENE ((Range) { 79625, 80425 })
#define SEQUENCE_WINDOW_SITE 80275

// The sequence in the window around SEQUENCE_WINDOW_SITE of the gene SEQUENCE_WINDOW_GENE, which is 800 bases long on the reference and ends within
// the context after the site, with the haplotype of one indel (mock_bcf1_t(position, ref, alt)). ref and alt NULL is the
// reference. Checks that the window has the expected start and holds hap_len bases, the length of the haplotype of the gene.
static void check_sequence_in_window(const uint64_t position, const char *ref, const char *alt, const int expected_start, const int hap_len) {
    Records rs;
    kv_init(rs);
    if (ref != NULL) kv_push(Record *, rs, Record_init(mock_bcf1_t(position, ref, alt), NULL));

    // The width the window is predicted with (see predict_splice_changes_for_phase): the model width, shorter or longer by the
    // indel where it is in the prediction window
    const int width = model_width(SEQUENCE_WINDOW_RADIUS) - haplotype_offset(&rs, range_around(SEQUENCE_WINDOW_SITE, SEQUENCE_WINDOW_RADIUS));
    const Range offsets = sequence_in_window(SEQUENCE_WINDOW_SITE, SEQUENCE_WINDOW_GENE, ref != NULL ? &rs : NULL, SEQUENCE_WINDOW_RADIUS);

    ck_assert_int_eq(offsets.start, expected_start);
    ck_assert_int_eq(width - offsets.start - offsets.end, hap_len);

    for (size_t i = 0; i < kv_size(rs); i++) Record_destroy(kv_A(rs, i));
    kv_destroy(rs);
}

START_TEST(check_sequence_in_window__the_reference) {
    check_sequence_in_window(0, NULL, NULL, 4400, 800);
} END_TEST

// The window of the haplotype whose prediction varied: the deletion before the site made the window take 5 bases past the end
// of the sequence, which were whatever was in memory there
START_TEST(check_sequence_in_window__a_deletion_before_the_site) {
    check_sequence_in_window(80254, "NNNNNN", "N", 4400, 795);
} END_TEST

START_TEST(check_sequence_in_window__an_insertion_before_the_site) {
    check_sequence_in_window(80254, "N", "NNNNNN", 4400, 805);
} END_TEST

START_TEST(check_sequence_in_window__a_deletion_after_the_site) {
    check_sequence_in_window(80299, "NNNNNN", "N", 4400, 795);
} END_TEST

START_TEST(check_sequence_in_window__a_deletion_after_the_prediction_window) {
    check_sequence_in_window(80399, "NNNNNN", "N", 4400, 795);
} END_TEST

Suite *splice_prediction_suite() {
    Suite *s = suite_create("Splice prediction");

    TCase *tc_helpers = mocked_tcase("prediction_helpers");
    tcase_add_test(tc_helpers, check_precomputed_reference_predictions__without_reference_scores_there_are_none);
    tcase_add_test(tc_helpers, check_filter_to_most_likely_prediction__keeps_the_largest_of_every_position);
    tcase_add_test(tc_helpers, check_filter_to_most_likely_prediction__ties_and_zeros);
    tcase_add_test(tc_helpers, check_transform_to_posterior_probabilities__scales_the_nonzero_values_to_sum_to_one);
    tcase_add_test(tc_helpers, check_transform_to_posterior_probabilities__unequal_values);
    tcase_add_test(tc_helpers, check_transform_to_posterior_probabilities__all_zero);
    tcase_add_test(tc_helpers, check_check_proximity__within_the_threshold_on_both_sides);
    tcase_add_test(tc_helpers, check_kv_push_new_record__a_record_is_added_once);
    suite_add_tcase(s, tc_helpers);

    TCase *tc_store = mocked_tcase("store_splice_sites");
    tcase_add_test(tc_store, check_store_splice_sites__a_lost_canonical_donor);
    tcase_add_test(tc_store, check_store_splice_sites__a_lost_canonical_acceptor);
    tcase_add_test(tc_store, check_store_splice_sites__the_only_site_in_the_window_is_the_canonical_one);
    tcase_add_test(tc_store, check_store_splice_sites__a_site_next_to_the_canonical_one_becomes_a_canonical_site);
    tcase_add_test(tc_store, check_store_splice_sites__another_record_joins_a_known_site);
    tcase_add_test(tc_store, check_store_splice_sites__a_far_site_is_cryptic);
    tcase_add_test(tc_store, check_store_splice_sites__many_cryptic_sites_grow_the_array);
    tcase_add_test(tc_store, check_store_splice_sites__an_unused_other_site_is_ignored);
    suite_add_tcase(s, tc_store);

    TCase *tc_init = mocked_tcase("init_splice_sites");
    tcase_add_test(tc_init, check_init_exon_splice_sites__forward_strand);
    tcase_add_test(tc_init, check_init_exon_splice_sites__reverse_strand);
    tcase_add_test(tc_init, check_init_exon_splice_sites__unknown_strand_has_no_sites);
    tcase_add_test(tc_init, check_init_splice_sites_for_phase__a_record_in_the_gene_gives_its_exons_sites);
    tcase_add_test(tc_init, check_init_splice_sites_for_phase__a_record_outside_the_gene_or_without_transcript_gives_none);
    suite_add_tcase(s, tc_init);

    TCase *tc_collect = mocked_tcase("collect_splice_probabilities");
    tcase_add_test(tc_collect, check_collect_splice_probabilities__a_lost_canonical_donor);
    tcase_add_test(tc_collect, check_collect_splice_probabilities__a_cryptic_donor_on_the_haplotype);
    tcase_add_test(tc_collect, check_collect_splice_probabilities__a_weak_canonical_donor_viable_on_the_haplotype);
    tcase_add_test(tc_collect, check_collect_splice_probabilities__no_predictions_gives_nothing);
    suite_add_tcase(s, tc_collect);

    TCase *tc_locate = mocked_tcase("locate_prediction_on_exon");
    tcase_add_test(tc_locate, check_locate_prediction_on_exon__a_donor_in_the_window);
    tcase_add_test(tc_locate, check_locate_prediction_on_exon__an_acceptor_in_the_window);
    tcase_add_test(tc_locate, check_locate_prediction_on_exon__an_exon_without_a_window_gives_false);
    suite_add_tcase(s, tc_locate);

    TCase *tc_window = mocked_tcase("sequence_in_window");
    tcase_add_test(tc_window, check_sequence_in_window__the_reference);
    tcase_add_test(tc_window, check_sequence_in_window__a_deletion_before_the_site);
    tcase_add_test(tc_window, check_sequence_in_window__an_insertion_before_the_site);
    tcase_add_test(tc_window, check_sequence_in_window__a_deletion_after_the_site);
    tcase_add_test(tc_window, check_sequence_in_window__a_deletion_after_the_prediction_window);
    suite_add_tcase(s, tc_window);

    return s;
}
