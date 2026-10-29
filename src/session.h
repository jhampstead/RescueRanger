#ifndef SESSION_HEADER
#define SESSION_HEADER

#include <stdio.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <cpliceai/predict.h>

#include "bcftools/gff.h"
#include "record.h"
#include "reference_scores.h"

typedef struct {
    Model     *models;
    gff_t     *gff;
    faidx_t   *fai;
    /* The precomputed reference predictions, or NULL: then the reference predictions are made from the sequence. */
    ReferenceScores *reference;
    htsFile   *fp_in;
    htsFile   *fp_out;
    /* Read from fp_in, extended with our INFO tags, then written to fp_out. */
    bcf_hdr_t *hdr;
    HdrFormat *vep_fmt;
    FILE      *scores_tsv;
    int        prediction_radius;
} Session;

/*
 * Releases whatever the session holds. Every field is guarded, so this also
 * unwinds a partially-completed Session_open.
 */
void Session_close(Session *s);

#endif
