#include "session.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "logging/log.h"
#include "phase.h"
#include "tags.h"
#include "utils.h"

int Session_open(const RescueOptions *opts, bcf_hdr_t *hdr, Session *s) {
    *s = (Session){0};

    if (!opts->model_dir || !opts->gff_fname || !opts->fasta_fname) {
        log_error("The models, the GFF and the FASTA are required");
        goto fail;
    }

    if (opts->prediction_radius < 1 || opts->prediction_radius > MAX_PREDICTION_RADIUS) {
        log_error("--prediction-radius must be between 1 and %d, got %d", MAX_PREDICTION_RADIUS, opts->prediction_radius);
        goto fail;
    }
    s->prediction_radius = opts->prediction_radius;

    AnnotationKind kind;
    switch (opts->annotation) {
        case RR_ANNOTATION_VEP: kind = ANNOTATION_VEP; break;
        case RR_ANNOTATION_BCFTOOLS: kind = ANNOTATION_BCFTOOLS; break;
        default:
            log_error("--annotation must be vep or bcftools, got %d", (int) opts->annotation);
            goto fail;
    }

    /* The genotypes are checked before anything is loaded: without them no record can be put on a haplotype. */
    const char *gt_tag = Genotype_tag_name(opts->gt_tag ? opts->gt_tag : "FORMAT/GT");
    if (gt_tag == NULL) {
        log_error("--gt-tag must be a FORMAT field, got %s", opts->gt_tag);
        goto fail;
    }
    if (bcf_hdr_nsamples(hdr) < 1) {
        log_error("The VCF has no samples, so there are no genotypes to phase the variants with");
        goto fail;
    }
    const int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, gt_tag);
    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) {
        log_error("The VCF header has no FORMAT/%s for the genotypes", gt_tag);
        goto fail;
    }
    if (strcmp(gt_tag, "GT") != 0 && bcf_hdr_id2type(hdr, BCF_HL_FMT, gt_id) != BCF_HT_STR) {
        log_error("FORMAT/%s must be of Type=String, holding genotypes like those of GT (0|1)", gt_tag);
        goto fail;
    }
    if (!(s->gt_tag = strdup(gt_tag))) {
        log_error("Out of memory");
        goto fail;
    }

    const bool use_reference_scores = opts->reference_scores && opts->reference_scores[0] != '\0';

    log_info("Loading SpliceAI models from: %s", opts->model_dir);
    if (!(s->models = load_models(opts->model_dir))) {
        log_error("Could not load SpliceAI models from: %s", opts->model_dir);
        goto fail;
    }

    /* gff_init and gff_parse report their own failures through bcftools' error(), which exits. */
    log_info("Loading GFF from: %s", opts->gff_fname);
    s->gff = gff_init(opts->gff_fname);
    gff_parse(s->gff);

    log_info("Loading FASTA from: %s", opts->fasta_fname);
    if (!(s->fai = fai_load(opts->fasta_fname))) {
        log_error("Could not load FASTA index for: %s", opts->fasta_fname);
        goto fail;
    }

    /* If they are given but cannot be used (another genome than they were made from), the run stops: predicting the reference from
     * the sequence instead would be quietly doing something else than was asked. Leaving the option out does that. A gene that is not
     * in them like the annotation's is predicted from the sequence, with a warning. */
    if (use_reference_scores) {
        log_info("Loading reference scores from: %s", opts->reference_scores);
        if (!(s->reference = ReferenceScores_open(opts->reference_scores, s->fai))) goto fail;
    } else {
        log_info("No reference scores given: the reference predictions are made from the sequence");
    }

    if (!(s->vep_fmt = parse_csq_hdr_fmt(hdr, kind, opts->csq_tag))) {
        log_error("Could not read the consequence annotation (%s); annotate with VEP or bcftools csq first.",
                  kind == ANNOTATION_VEP ? "vep" : "bcftools");
        goto fail;
    }

    if (add_rescue_tags(hdr) != 0) goto fail;
    s->hdr = hdr; /* Only once nothing can fail any more: on failure the header stays the caller's */

    return EXIT_SUCCESS;

fail:
    Session_close(s);
    return EXIT_FAILURE;
}

void Session_close(Session *s) {
    if (s->scores_tsv) fclose(s->scores_tsv);
    if (s->fp_out)     hts_close(s->fp_out);
    if (s->fp_in)      hts_close(s->fp_in);
    if (s->vep_fmt)    destroy_hdr_fmt(s->vep_fmt);
    if (s->hdr)        bcf_hdr_destroy(s->hdr);
    if (s->reference)  ReferenceScores_close(s->reference); /* Before the index it reads the sequence with */
    if (s->fai)        fai_destroy(s->fai);
    if (s->gff)        gff_destroy(s->gff);
    if (s->models)     destroy_models(s->models);
    free(s->gt_tag);
    *s = (Session){0};
}
