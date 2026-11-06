#ifndef SESSION_HEADER
#define SESSION_HEADER

#include <stdio.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <cpliceai/predict.h>
#include <rescue/rescue.h>

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
    /* The raw scores, or NULL: then they are not written. */
    FILE      *scores_tsv;
    int        prediction_radius;
    /* The name of the FORMAT field the genotypes are read from, such as GT (see Genotype_read). Owned by the session. */
    char      *gt_tag;
} Session;

/*
 * Loads the models, the GFF, the FASTA and the reference scores, reads the consequence annotation from hdr and adds the output
 * INFO tags to it. hdr is the header of the input, opened by the caller; on success the session points at it (s->hdr). The files
 * (fp_in, fp_out, scores_tsv) are opened by the caller, because the program and the bcftools plugin open them differently.
 * Returns EXIT_SUCCESS, or EXIT_FAILURE after logging why and closing what was opened; then hdr is not freed.
 */
int Session_open(const RescueOptions *opts, bcf_hdr_t *hdr, Session *s);

/*
 * Releases whatever the session holds. Every field is guarded, so this also
 * unwinds a partially-completed Session_open.
 */
void Session_close(Session *s);

#endif
