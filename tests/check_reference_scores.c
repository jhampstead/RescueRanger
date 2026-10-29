#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <htslib/faidx.h>

#include "reference_scores.h"

#include "check_reference_scores.h"

// This file includes the gene and reference headers of cpliceai (when this build has them) and so not record.h, whose Gene they
// clash with. The tests are about ReferenceScores_window and ReferenceScores_open only.

#ifdef CPLICEAI_HAS_REFERENCE_SCORES

#include <cpliceai/gene_regions.h>
#include <cpliceai/reference.h>
#include <cpliceai/utils.h>

/*
* A small genome and the reference scores made from it, written the way cpliceai_reference writes them: the contig chr1 of 300
* bases, and the genes G1 at [100, 160) on the forward strand and G2 at [200, 280) on the reverse strand. The scores that are not
* zero (the file only keeps the positions where the acceptor or the donor is at least 1e-3), as offsets from the start of the gene:
*   G1:  10 acceptor 0.5, donor 0.1   40 donor 0.9   59 acceptor 0.25
*   G2:   0 acceptor 0.7, donor 0.2   79 acceptor 0.3, donor 0.3
* There is no gene regions file: the fingerprint of one that the file holds is not checked.
*/
typedef struct Fixture {
    char fasta[64], reference[64];
    faidx_t *fai;
} Fixture;

static void write_text(const char *path, const char *text) {
    FILE *out = fopen(path, "w");
    ck_assert_ptr_nonnull(out);
    fputs(text, out);
    fclose(out);
}

static void temp_path(char *path, const char *prefix) {
    snprintf(path, 64, "%s_XXXXXX", prefix);
    const int fd = mkstemp(path);
    ck_assert_int_ge(fd, 0);
    close(fd);
}

static void add_gene(Reference *ref, const char *name, const char strand, const int64_t start, const int64_t end,
                     const PositionScore *scores, const int n_scores) {
    Region region = {0};
    region.chunk_size = CHUNK_SIZE;
    region.size = end - start;
    region.strand = strand;
    region.chunk_start = ref->n_chunks;
    region.name_start = ref->region_names_len;
    reference_add_region(region, name, start, end, ref);

    Chunk chunk = {0};
    chunk.scores_start = ref->n_scores;
    reference_add_chunk(chunk, ref);
    for (int i = 0; i < n_scores; i++) reference_add_score(scores[i], ref);
}

static Fixture make_fixture(void) {
    Fixture f = {0};

    temp_path(f.fasta, "check_reference_fa");
    char fasta[400] = ">chr1\n";
    memset(fasta + strlen(fasta), 'A', 300);
    strcpy(fasta + 306, "\n");
    write_text(f.fasta, fasta);
    f.fai = fai_load(f.fasta);   // Writes the index next to the FASTA
    ck_assert_ptr_nonnull(f.fai);

    Reference ref;
    ck_assert_int_eq(reference_alloc(&ref, 1, 5, 3, 9, 3, 16), EXIT_SUCCESS);

    Contig contig = {0};
    contig.region_start = ref.n_genes;
    contig.name_start = ref.contig_names_len;
    reference_add_contig(contig, "chr1", &ref);

    const PositionScore g1[] = { { 10, 0.5f, 0.1f }, { 40, 0.0f, 0.9f }, { 59, 0.25f, 0.0f } };
    const PositionScore g2[] = { { 0, 0.7f, 0.2f }, { 79, 0.3f, 0.3f } };
    const PositionScore g3[] = { { 5, 0.4f, 0.0f } };
    add_gene(&ref, "G1", '+', 100, 160, g1, 3);
    add_gene(&ref, "G2", '-', 200, 280, g2, 2);
    add_gene(&ref, "G3", '+', 60, 90, g3, 1);   // The annotation will have it one base further on

    ref.fasta_digest = fasta_digest(f.fai);
    ref.regions_digest = 0x1234;   // Whatever it is: nothing checks it

    temp_path(f.reference, "check_reference_scores");
    ck_assert_int_eq(reference_write(f.reference, &ref), EXIT_SUCCESS);
    reference_free(&ref);

    return f;
}

static void destroy_fixture(Fixture *f) {
    char fai_path[80];
    snprintf(fai_path, sizeof(fai_path), "%s.fai", f->fasta);

    fai_destroy(f->fai);
    unlink(f->fasta); unlink(fai_path); unlink(f->reference);
}

#define NEITHER 0
#define ACCEPTOR 1
#define DONOR 2

static void check_position(const float *window, const int index, const float acceptor, const float donor) {
    ck_assert_msg(fabsf(window[index * NUM_SCORES + ACCEPTOR] - acceptor) < 1e-6f, "position %d: acceptor %f, expected %f", index, window[index * NUM_SCORES + ACCEPTOR], acceptor);
    ck_assert_msg(fabsf(window[index * NUM_SCORES + DONOR] - donor) < 1e-6f, "position %d: donor %f, expected %f", index, window[index * NUM_SCORES + DONOR], donor);
    ck_assert_msg(fabsf(window[index * NUM_SCORES + NEITHER] - (1.0f - acceptor - donor)) < 1e-6f, "position %d: neither %f, expected %f", index, window[index * NUM_SCORES + NEITHER], 1.0f - acceptor - donor);
}

START_TEST(check_reference_scores__is_supported_by_this_build) {
    ck_assert_int_eq(ReferenceScores_supported(), 1);
} END_TEST

START_TEST(check_reference_scores__a_window_inside_the_gene_holds_its_scores) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);
    ck_assert_ptr_nonnull(scores);

    // The 11 bases 100 to 110 of G1, the offsets 0 to 10: only the last has scores
    float *window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 105, 5);
    ck_assert_ptr_nonnull(window);
    for (int i = 0; i < 10; i++) check_position(window, i, 0.0f, 0.0f);
    check_position(window, 10, 0.5f, 0.1f);
    free(window);

    // Around the donor at 40: the offset 40 is the position 140
    window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 140, 2);
    ck_assert_ptr_nonnull(window);
    check_position(window, 0, 0.0f, 0.0f);
    check_position(window, 2, 0.0f, 0.9f);
    check_position(window, 4, 0.0f, 0.0f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

// The window begins before the gene: the bases before it have no splice site, and the gene begins where its offset 0 is
START_TEST(check_reference_scores__a_window_over_the_start_of_the_gene_has_no_sites_before_it) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    // The positions 83 to 113: the offset 0 (the position 100) is the index 17, and the offset 10 the index 27
    float *window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 98, 15);
    ck_assert_ptr_nonnull(window);
    for (int i = 0; i < 27; i++) check_position(window, i, 0.0f, 0.0f);
    check_position(window, 27, 0.5f, 0.1f);
    check_position(window, 28, 0.0f, 0.0f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

START_TEST(check_reference_scores__a_window_over_the_end_of_the_gene_has_no_sites_after_it) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    // The positions 153 to 163: the last base of the gene, 159, is the index 6
    float *window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 158, 5);
    ck_assert_ptr_nonnull(window);
    for (int i = 0; i < 6; i++) check_position(window, i, 0.0f, 0.0f);
    check_position(window, 6, 0.25f, 0.0f);
    for (int i = 7; i < 11; i++) check_position(window, i, 0.0f, 0.0f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

// The reverse strand: the scores are in ascending order of position like the predictions of the sequence
START_TEST(check_reference_scores__a_reverse_gene_is_in_ascending_order) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    float *window = ReferenceScores_window(scores, "chr1", "G2", 200, 280, '-', 200, 3);
    ck_assert_ptr_nonnull(window);
    for (int i = 0; i < 3; i++) check_position(window, i, 0.0f, 0.0f);
    check_position(window, 3, 0.7f, 0.2f);
    check_position(window, 4, 0.0f, 0.0f);
    free(window);

    // The last base of the gene, 279
    window = ReferenceScores_window(scores, "chr1", "G2", 200, 280, '-', 279, 3);
    ck_assert_ptr_nonnull(window);
    check_position(window, 3, 0.3f, 0.3f);
    check_position(window, 4, 0.0f, 0.0f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

START_TEST(check_reference_scores__calls_for_different_genes_in_turn_give_the_scores_of_each) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    for (int round = 0; round < 3; round++) {
        float *g1 = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 110, 0);
        float *g2 = ReferenceScores_window(scores, "chr1", "G2", 200, 280, '-', 200, 0);
        ck_assert_ptr_nonnull(g1);
        ck_assert_ptr_nonnull(g2);
        check_position(g1, 0, 0.5f, 0.1f);
        check_position(g2, 0, 0.7f, 0.2f);
        free(g1); free(g2);
    }

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

// A gene that is not the annotation's gene is not used, and the next call for a gene that is works as before
START_TEST(check_reference_scores__a_gene_that_is_not_in_the_scores_gives_none) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 100, 161, '+', 110, 5));    // The end differs
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 101, 160, '+', 110, 5));    // The start differs
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 100, 160, '-', 110, 5));    // The strand differs
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G4", 100, 160, '+', 110, 5));    // No such gene
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr2", "G1", 100, 160, '+', 110, 5));    // No such contig
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G2", 100, 160, '-', 110, 5));    // A gene with the coordinates of another

    // Asking again is no different, and the genes that are there are not harmed
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 100, 161, '+', 110, 5));

    float *window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 110, 0);
    ck_assert_ptr_nonnull(window);
    check_position(window, 0, 0.5f, 0.1f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

// The case this check is for: the scores come from an annotation where the gene begins one base before. The gene is not read from
// them, and the genes that are the same are
START_TEST(check_reference_scores__scores_made_from_another_annotation_are_not_used_for_the_genes_that_differ) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);
    ck_assert_ptr_nonnull(scores);

    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G3", 59, 90, '+', 70, 5));      // The annotation's G3 begins at 59, the scores' at 60

    float *window = ReferenceScores_window(scores, "chr1", "G3", 60, 90, '+', 65, 0);          // Its coordinates in the scores
    ck_assert_ptr_nonnull(window);
    check_position(window, 0, 0.4f, 0.0f);
    free(window);

    window = ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 110, 0);
    ck_assert_ptr_nonnull(window);
    check_position(window, 0, 0.5f, 0.1f);
    free(window);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

START_TEST(check_reference_scores__counts_the_windows_and_the_genes_that_were_refused) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);
    long windows, refused;

    ReferenceScores_counts(scores, &windows, &refused);
    ck_assert_int_eq(windows, 0);
    ck_assert_int_eq(refused, 0);

    free(ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 110, 1));
    free(ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 120, 1));
    free(ReferenceScores_window(scores, "chr1", "G2", 200, 280, '-', 210, 1));

    // The same gene is one gene however often it is asked for, and a gene with other coordinates is another
    ReferenceScores_window(scores, "chr1", "G1", 100, 161, '+', 110, 1);
    ReferenceScores_window(scores, "chr1", "G1", 100, 161, '+', 120, 1);
    ReferenceScores_window(scores, "chr1", "G4", 100, 160, '+', 110, 1);

    ReferenceScores_counts(scores, &windows, &refused);
    ck_assert_int_eq(windows, 3);
    ck_assert_int_eq(refused, 2);

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

// A window that has nothing of the gene in it is a bug of the caller, and gives none
START_TEST(check_reference_scores__a_window_away_from_the_gene_gives_none) {
    Fixture f = make_fixture();
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);

    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 10, 5));
    ck_assert_ptr_null(ReferenceScores_window(scores, "chr1", "G1", 100, 160, '+', 250, 5));

    ReferenceScores_close(scores);
    destroy_fixture(&f);
} END_TEST

START_TEST(check_reference_scores__without_scores_there_is_no_window) {
    ck_assert_ptr_null(ReferenceScores_window(NULL, "chr1", "G1", 100, 160, '+', 110, 5));
    ReferenceScores_close(NULL);   // And closing nothing is fine

    long windows = 5, refused = 5;
    ReferenceScores_counts(NULL, &windows, &refused);
    ck_assert_int_eq(windows, 0);
    ck_assert_int_eq(refused, 0);
} END_TEST

// Scores for another FASTA are refused: they would give wrong predictions silently
START_TEST(check_reference_scores__scores_made_for_another_fasta_are_refused) {
    Fixture f = make_fixture();

    // A FASTA with a contig of another length
    char other_fasta[64];
    temp_path(other_fasta, "check_reference_other_fa");
    write_text(other_fasta, ">chr1\nAAAAAAAAAA\n");
    faidx_t *other = fai_load(other_fasta);
    ck_assert_ptr_nonnull(other);

    ck_assert_ptr_null(ReferenceScores_open(f.reference, other));

    char other_fai[80];
    snprintf(other_fai, sizeof(other_fai), "%s.fai", other_fasta);
    fai_destroy(other);
    unlink(other_fasta); unlink(other_fai);

    // The one they were made from is accepted
    ReferenceScores *scores = ReferenceScores_open(f.reference, f.fai);
    ck_assert_ptr_nonnull(scores);
    ReferenceScores_close(scores);

    destroy_fixture(&f);
} END_TEST

START_TEST(check_reference_scores__a_missing_file_or_one_that_is_not_scores_is_refused) {
    Fixture f = make_fixture();

    ck_assert_ptr_null(ReferenceScores_open("check_reference_no_such_file.bin", f.fai));

    // The FASTA is not reference scores, and it is shorter than the header of one, or fails on its magic
    ck_assert_ptr_null(ReferenceScores_open(f.fasta, f.fai));

    destroy_fixture(&f);
} END_TEST

#else

START_TEST(check_reference_scores__a_build_without_them_says_so) {
    ck_assert_int_eq(ReferenceScores_supported(), 0);
    ck_assert_ptr_null(ReferenceScores_open("scores.bin", NULL));
    ck_assert_ptr_null(ReferenceScores_window(NULL, "chr1", "G1", 100, 160, '+', 110, 5));
    ReferenceScores_close(NULL);

    long windows = 5, refused = 5;
    ReferenceScores_counts(NULL, &windows, &refused);
    ck_assert_int_eq(windows, 0);
    ck_assert_int_eq(refused, 0);
} END_TEST

#endif

Suite *reference_scores_suite() {
    Suite *s = suite_create("Reference scores");
    TCase *tc = tcase_create("reference_scores");

#ifdef CPLICEAI_HAS_REFERENCE_SCORES
    tcase_add_test(tc, check_reference_scores__is_supported_by_this_build);
    tcase_add_test(tc, check_reference_scores__a_window_inside_the_gene_holds_its_scores);
    tcase_add_test(tc, check_reference_scores__a_window_over_the_start_of_the_gene_has_no_sites_before_it);
    tcase_add_test(tc, check_reference_scores__a_window_over_the_end_of_the_gene_has_no_sites_after_it);
    tcase_add_test(tc, check_reference_scores__a_reverse_gene_is_in_ascending_order);
    tcase_add_test(tc, check_reference_scores__calls_for_different_genes_in_turn_give_the_scores_of_each);
    tcase_add_test(tc, check_reference_scores__a_gene_that_is_not_in_the_scores_gives_none);
    tcase_add_test(tc, check_reference_scores__scores_made_from_another_annotation_are_not_used_for_the_genes_that_differ);
    tcase_add_test(tc, check_reference_scores__counts_the_windows_and_the_genes_that_were_refused);
    tcase_add_test(tc, check_reference_scores__a_window_away_from_the_gene_gives_none);
    tcase_add_test(tc, check_reference_scores__without_scores_there_is_no_window);
    tcase_add_test(tc, check_reference_scores__scores_made_for_another_fasta_are_refused);
    tcase_add_test(tc, check_reference_scores__a_missing_file_or_one_that_is_not_scores_is_refused);
#else
    tcase_add_test(tc, check_reference_scores__a_build_without_them_says_so);
#endif

    suite_add_tcase(s, tc);
    return s;
}
