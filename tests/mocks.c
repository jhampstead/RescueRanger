#include "mocks.h"

bcf1_t * mock_bcf1_t(const uint64_t position, char *ref, char *alt) {
    bcf1_t *v = bcf_init();

    v->pos = position;
    v->rlen = strlen(ref);
    v->n_allele = 2;
    v->d = (bcf_dec_t) {0};

    v->d.allele = malloc(2 * sizeof(char *));
    v->d.allele[0] = ref;
    v->d.allele[1] = alt;

    return v;
}

