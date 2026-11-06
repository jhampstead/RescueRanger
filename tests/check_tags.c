#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/vcf.h>

#include "klib/kvec.h"
#include "mocks.h"
#include "record.h"
#include "tags.h"

#include "check_tags.h"

/*
* add_rescue_tags: the header gets the INFO lines of every output tag, in the order of INFO_TAGS.
*/

START_TEST(check_add_rescue_tags__every_output_tag_is_in_the_header) {
    bcf_hdr_t *hdr = mock_hdr();

    ck_assert_int_eq(add_rescue_tags(hdr), 0);
    ck_assert_int_eq(bcf_hdr_sync(hdr), 0);

    const char *tags[] = INFO_TAGS;
    for (size_t i = 0; i < sizeof(tags) / sizeof(tags[0]); i++) {
        ck_assert_msg(bcf_hdr_id2int(hdr, BCF_DT_ID, tags[i]) >= 0, "tag %s is not in the header", tags[i]);
    }

    bcf_hdr_destroy(hdr);
} END_TEST

START_TEST(check_add_rescue_tags__tags_and_descriptions_match) {
    const char *tags[] = INFO_TAGS;
    const char *descriptions[] = INFO_TAG_DESCS;

    ck_assert_int_eq(sizeof(tags) / sizeof(tags[0]), sizeof(descriptions) / sizeof(descriptions[0]));

    // Every description declares the INFO field of the tag at the same place
    for (size_t i = 0; i < sizeof(tags) / sizeof(tags[0]); i++) {
        char id[64];
        snprintf(id, sizeof(id), "ID=%s,", tags[i]);
        ck_assert_msg(strstr(descriptions[i], id) != NULL, "description %zu does not declare %s", i, tags[i]);
    }
} END_TEST

/*
* get_prediction_flags: the names of the set bits of an enum PredictionFlags value, in bit order, joined with ':'.
*/

static char *flag_names(const uint32_t flags) {
    kstring_t names = {0};
    get_prediction_flags(flags, &names);
    return names.s; // Freed by the caller
}

START_TEST(check_get_prediction_flags__one_flag) {
    char *names = flag_names(UNPHASED_VARIANT_IN_PREDICTION);
    ck_assert_str_eq(names, "unphased_variant_in_prediction");
    free(names);
} END_TEST

START_TEST(check_get_prediction_flags__flags_come_in_bit_order_joined_with_colons) {
    // Passed in the other order than they come out
    char *names = flag_names(SKIPPED_VARIANT_IN_PREDICTION | UNPHASED_VARIANT_IN_PREDICTION | UTR3_VARIANT | LAST_EXON);
    ck_assert_str_eq(names, "unphased_variant_in_prediction:skipped_variant_in_prediction:3utr_variant:last_exon");
    free(names);
} END_TEST

// The name list has to follow the enum: a flag added to the enum without a name, or a name in the wrong place, shows here
START_TEST(check_get_prediction_flags__every_flag_has_its_own_name) {
    static const struct { uint32_t flag; const char *name; } expected[] = {
        { UNPHASED_VARIANT_IN_PREDICTION, "unphased_variant_in_prediction" },
        { SKIPPED_VARIANT_IN_PREDICTION, "skipped_variant_in_prediction" },
        { VARIANT_LARGER_THAN_SPLICE_WINDOW, "variant_larger_than_splice_window" },
        { RESCUED_BY_SPLICE_FRAMESHIFT, "rescued_by_splice_frameshift" },
        { POSSIBLE_REINITIATION, "possible_reinitiation" },
        { MANUFACTURED_SPLICE_RESCUE, "manufactured_splice_rescue" },
        { HAPLOTYPE_SPECIFIC_DONOR, "haplotype_specific_donor" },
        { HAPLOTYPE_SPECIFIC_ACCEPTOR, "haplotype_specific_acceptor" },
        { CRYPTIC_DONOR, "cryptic_donor" },
        { CRYPTIC_ACCEPTOR, "cryptic_acceptor" },
        { WEAK_CANONICAL_DONOR, "weak_canonical_donor" },
        { WEAK_CANONICAL_ACCEPTOR, "weak_canonical_acceptor" },
        { ABLATED_CANONICAL_DONOR, "ablated_canonical_donor" },
        { ABLATED_CANONICAL_ACCEPTOR, "ablated_canonical_acceptor" },
        { ABLATED_CRYPTIC_DONOR, "ablated_cryptic_donor" },
        { ABLATED_CRYPTIC_ACCEPTOR, "ablated_cryptic_acceptor" },
        { UPSTREAM_PTV, "upstream_ptv" },
        { UTR5_VARIANT, "5utr_variant" },
        { UTR3_VARIANT, "3utr_variant" },
        { LAST_EXON, "last_exon" },
        { SINGLE_EXON, "single_exon" },
        { NON_CONSTITUTIVE, "non_constitutive" },
        { INFRAME, "in_frame" },
        { NMD_START_PROXIMAL, "start_proximal" },
        { NMD_LONG_EXON, "long_exon" },
        { NMD_DISTAL, "distal" },
    };

    for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); i++) {
        char *names = flag_names(expected[i].flag);
        ck_assert_msg(strcmp(names, expected[i].name) == 0, "flag %u gives \"%s\", expected \"%s\"", expected[i].flag, names, expected[i].name);
        free(names);
    }
} END_TEST

/*
* generate_rescue_type_tag: rescue|type|probability|affected_bases for every rescue type of an isoforms struct, joined with ':'.
*/

// An Isoforms with the given rescue types as calculate_rescue_probabilities leaves them. The caller frees it with Isoforms_destroy.
static Isoforms isoforms_with_rescue_types(const int n, const uint8_t *rescues, const enum RescueType *types, const float *probabilities, const uint32_t *affected_bases) {
    Isoforms isoforms = {0};
    isoforms.valid = 1;
    isoforms.n_rescue_types = n;
    isoforms.m_rescue_types = n;
    kv_init(isoforms.isoforms);

    RescueProbabilities *probs = malloc(sizeof(RescueProbabilities));
    probs->n_isoforms = calloc(n, sizeof(uint32_t));
    probs->rescues = malloc(n * sizeof(uint8_t));
    probs->rescue_types = malloc(n * sizeof(enum RescueType));
    probs->rescue_probabilities = malloc(n * sizeof(float));
    probs->affected_bases = malloc(n * sizeof(uint32_t));
    for (int i = 0; i < n; i++) {
        probs->rescues[i] = rescues[i];
        probs->rescue_types[i] = types[i];
        probs->rescue_probabilities[i] = probabilities[i];
        probs->affected_bases[i] = affected_bases[i];
    }
    isoforms.rescue_type_probabilities = probs;

    return isoforms;
}

START_TEST(check_generate_rescue_type_tag__one_type) {
    const uint8_t rescues[] = { 1 };
    const enum RescueType types[] = { EXON_SKIPPING };
    const float probabilities[] = { 0.75f };
    const uint32_t affected[] = { 3 };
    Isoforms isoforms = isoforms_with_rescue_types(1, rescues, types, probabilities, affected);

    kstring_t tag = generate_rescue_type_tag(&isoforms);
    ck_assert_str_eq(tag.s, "1|exon_skipping|0.7500|3");

    free(tag.s);
    Isoforms_destroy(&isoforms);
} END_TEST

START_TEST(check_generate_rescue_type_tag__types_are_joined_with_colons) {
    const uint8_t rescues[] = { 1, 0 };
    const enum RescueType types[] = { INTRON_RETENTION, NO_RESCUE };
    const float probabilities[] = { 0.5f, 0.25f };
    const uint32_t affected[] = { 21, 0 };
    Isoforms isoforms = isoforms_with_rescue_types(2, rescues, types, probabilities, affected);

    kstring_t tag = generate_rescue_type_tag(&isoforms);
    ck_assert_str_eq(tag.s, "1|intron_retention|0.5000|21:0|no_rescue|0.2500|0");

    free(tag.s);
    Isoforms_destroy(&isoforms);
} END_TEST

START_TEST(check_generate_rescue_type_tag__every_rescue_type_has_a_name) {
    static const struct { enum RescueType type; const char *name; } expected[] = {
        { NO_RESCUE, "no_rescue" },
        { COMPLETE_RESCUE, "complete_rescue" },
        { POSSIBLE_REINITIATION_SITE, "possible_reinitiation" },
        { FRAME_RESTORING_INDEL, "frame_restoring_indel" },
        { PARTIAL_EXON_SKIPPING, "partial_exon_skipping" },
        { PARTIAL_INTRON_RETENTION, "partial_intron_retention" },
        { EXON_SKIPPING, "exon_skipping" },
        { INTRON_RETENTION, "intron_retention" },
    };

    for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); i++) {
        const uint8_t rescues[] = { 0 };
        const enum RescueType types[] = { expected[i].type };
        const float probabilities[] = { 1.0f };
        const uint32_t affected[] = { 0 };
        Isoforms isoforms = isoforms_with_rescue_types(1, rescues, types, probabilities, affected);

        kstring_t tag = generate_rescue_type_tag(&isoforms);
        char want[64];
        snprintf(want, sizeof(want), "0|%s|1.0000|0", expected[i].name);
        ck_assert_str_eq(tag.s, want);

        free(tag.s);
        Isoforms_destroy(&isoforms);
    }
} END_TEST

START_TEST(check_generate_rescue_type_tag__no_types_gives_an_empty_tag) {
    Isoforms isoforms = {0};
    isoforms.valid = 1;
    kv_init(isoforms.isoforms);

    kstring_t tag = generate_rescue_type_tag(&isoforms);
    ck_assert_int_eq(tag.l, 0);

    free(tag.s);
    kv_destroy(isoforms.isoforms);
} END_TEST

/*
* format_gene_scores: one SpliceAI-style value per gene and phase: genes separated by ',', phases by '&'.
*/

// A record on chr1 with the given genes, as the loader gives it. The caller frees it with Record_destroy; the genes stay owned by the test.
static Record *record_with_genes(Gene *genes, const int n_genes) {
    bcf_hdr_t *hdr = mock_hdr();
    bcf1_t *v = mock_genotyped_bcf1_t(hdr, 100, "A,G", 0, 1, true, -1);
    Record *r = Record_init(v, NULL);
    bcf_destroy(v);
    bcf_hdr_destroy(hdr);

    r->genes = malloc(sizeof(Genes));
    kv_init(*r->genes);
    for (int i = 0; i < n_genes; i++) kv_push(Gene *, *r->genes, &genes[i]);

    return r;
}

START_TEST(check_format_gene_scores__a_score_and_a_missing_phase) {
    Gene genes[1] = { { .name = "GENEA" } };
    Record *r = record_with_genes(genes, 1);

    Score first = { .alt = "G", .gene = "GENEA", .ag = 0.1f, .al = 0.25f, .dg = 0.5f, .dl = 0.75f, .ag_idx = 1, .al_idx = -2, .dg_idx = 3, .dl_idx = 4 };
    Score *scores[NUM_PHASES] = { &first, NULL };

    kstring_t out = {0};
    format_gene_scores(r, scores, &out);

    ck_assert_str_eq(out.s, "G|GENEA|0.10|0.25|0.50|0.75|1|-2|3|4&G|GENEA|.|.|.|.|.|.|.|.");

    free(out.s);
    kv_destroy(*r->genes); free(r->genes); r->genes = NULL;
    Record_destroy(r);
} END_TEST

START_TEST(check_format_gene_scores__genes_are_separated_by_commas) {
    Gene genes[2] = { { .name = "GENEA" }, { .name = "GENEB" } };
    Record *r = record_with_genes(genes, 2);

    Score *scores[2 * NUM_PHASES] = { NULL, NULL, NULL, NULL };

    kstring_t out = {0};
    format_gene_scores(r, scores, &out);

    ck_assert_str_eq(out.s, "G|GENEA|.|.|.|.|.|.|.|.&G|GENEA|.|.|.|.|.|.|.|.,G|GENEB|.|.|.|.|.|.|.|.&G|GENEB|.|.|.|.|.|.|.|.");

    free(out.s);
    kv_destroy(*r->genes); free(r->genes); r->genes = NULL;
    Record_destroy(r);
} END_TEST

Suite *tags_suite() {
    Suite *s = suite_create("Tags");

    TCase *tc_header = tcase_create("add_rescue_tags");
    tcase_add_test(tc_header, check_add_rescue_tags__every_output_tag_is_in_the_header);
    tcase_add_test(tc_header, check_add_rescue_tags__tags_and_descriptions_match);
    suite_add_tcase(s, tc_header);

    TCase *tc_flags = tcase_create("get_prediction_flags");
    tcase_add_test(tc_flags, check_get_prediction_flags__one_flag);
    tcase_add_test(tc_flags, check_get_prediction_flags__flags_come_in_bit_order_joined_with_colons);
    tcase_add_test(tc_flags, check_get_prediction_flags__every_flag_has_its_own_name);
    suite_add_tcase(s, tc_flags);

    TCase *tc_rescue = tcase_create("generate_rescue_type_tag");
    tcase_add_test(tc_rescue, check_generate_rescue_type_tag__one_type);
    tcase_add_test(tc_rescue, check_generate_rescue_type_tag__types_are_joined_with_colons);
    tcase_add_test(tc_rescue, check_generate_rescue_type_tag__every_rescue_type_has_a_name);
    tcase_add_test(tc_rescue, check_generate_rescue_type_tag__no_types_gives_an_empty_tag);
    suite_add_tcase(s, tc_rescue);

    TCase *tc_scores = tcase_create("format_gene_scores");
    tcase_add_test(tc_scores, check_format_gene_scores__a_score_and_a_missing_phase);
    tcase_add_test(tc_scores, check_format_gene_scores__genes_are_separated_by_commas);
    suite_add_tcase(s, tc_scores);

    return s;
}
