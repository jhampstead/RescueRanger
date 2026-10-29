#include "mocks.h"

#include <stdio.h>
#include <stdlib.h>

// Every record mock_bcf1_t made and mocks_free_all has not freed yet
static bcf1_t **live_mocks = NULL;
static size_t n_live_mocks = 0;
static size_t m_live_mocks = 0;

static void remember_mock(bcf1_t *v) {
    if (n_live_mocks == m_live_mocks) {
        m_live_mocks = m_live_mocks ? 2 * m_live_mocks : 16;
        live_mocks = realloc(live_mocks, m_live_mocks * sizeof(bcf1_t *));
        if (live_mocks == NULL) abort();
    }

    live_mocks[n_live_mocks++] = v;
}

void mocks_free_all(void) {
    for (size_t i = 0; i < n_live_mocks; i++) bcf_destroy(live_mocks[i]);

    free(live_mocks);
    live_mocks = NULL;
    n_live_mocks = m_live_mocks = 0;
}

size_t mocks_live_count(void) {
    return n_live_mocks;
}

TCase *mocked_tcase(const char *name) {
    TCase *tc = tcase_create(name);
    tcase_add_checked_fixture(tc, NULL, mocks_free_all);
    return tc;
}

bcf1_t * mock_bcf1_t(const uint64_t position, const char *ref, const char *alt) {
    bcf1_t *v = bcf_init();
    remember_mock(v);

    v->pos = position;
    v->rlen = strlen(ref);
    v->n_allele = 2;
    v->d = (bcf_dec_t) {0};

    // The allele pointers are stored, not copied, and nothing writes through them. The field is not const in htslib.
    v->d.allele = malloc(2 * sizeof(char *));
    v->d.allele[0] = (char *) ref;
    v->d.allele[1] = (char *) alt;

    return v;
}

// The htslib calls that build the header are taken to return 0 on success. A header that silently misses a line would make the
// tests fail in confusing ways, so a failure stops the test run here.
static void require_success(const int result, const char *call) {
    if (result != 0) {
        fprintf(stderr, "mock_hdr: %s failed (returned %d)\n", call, result);
        abort();
    }
}

bcf_hdr_t *mock_hdr(void) {
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    require_success(bcf_hdr_append(hdr, "##contig=<ID=chr1,length=1000000>"), "bcf_hdr_append (contig)");
    require_success(bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"), "bcf_hdr_append (INFO/AF)");
    require_success(bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"), "bcf_hdr_append (FORMAT/GT)");
    require_success(bcf_hdr_add_sample(hdr, "sample"), "bcf_hdr_add_sample");
    require_success(bcf_hdr_sync(hdr), "bcf_hdr_sync");

    return hdr;
}

bcf1_t *mock_genotyped_bcf1_t(bcf_hdr_t *hdr, const int64_t position, const char *alleles, const int gt0, const int gt1, const bool phased, const float af) {
    bcf1_t *v = bcf_init();

    v->rid = bcf_hdr_name2id(hdr, "chr1");
    v->pos = position;
    v->n_sample = bcf_hdr_nsamples(hdr);
    bcf_update_alleles_str(hdr, v, alleles);

    // The phase flag is on the second allele, like it is when the genotype is parsed from "a|b"
    int32_t genotype[2] = { bcf_gt_unphased(gt0), phased ? bcf_gt_phased(gt1) : bcf_gt_unphased(gt1) };
    bcf_update_genotypes(hdr, v, genotype, 2);

    if (af >= 0) {
        float value = af;
        bcf_update_info_float(hdr, v, "AF", &value, 1);
    }

    return v;
}

