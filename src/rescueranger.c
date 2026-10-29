#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include <klib/kvec.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <cpliceai/predict.h>

#include "bcftools/gff.h"

#include "logging/log.h"
#include "isoforms.h"
#include "phase.h"
#include "record.h"
#include "record_block.h"
#include "session.h"
#include "splice_prediction.h"
#include "tags.h"
#include "triggers.h"
#include "utils.h"

#define REQUIRED_ARGS \
    REQUIRED_STRING_ARG(variants, "variants", "VCF file containing variants to predict for using SpliceAI") \
    REQUIRED_STRING_ARG(model_dir, "model_dir", "Directory containing SpliceAI models") \
    REQUIRED_STRING_ARG(gff, "gff_file", "GFF file matching the provided VCF file") \
    REQUIRED_STRING_ARG(fasta, "fasta_file", "Reference fasta file matching the provided VCF file") \
    REQUIRED_STRING_ARG(variants_output, "output", "Annotated VCF") \
    REQUIRED_STRING_ARG(scores_output, "output", "TSV with raw SpliceAI scores")

#define OPTIONAL_ARGS \
    OPTIONAL_STRING_ARG(af_tag, "INFO/AF", "--af-tag", "tag", "Allele frequency tag") \
    OPTIONAL_STRING_ARG(gt_tag, "FORMAT/GT", "--gt-tag", "tag", "Phased genotype tag") \
    OPTIONAL_DOUBLE_ARG(af_threshold, 0.01, "--af-threshold", "threshold", "Allele frequency threshold", 3) \
    OPTIONAL_STRING_ARG(annotation, "vep", "--annotation", "vep|bcftools", "Program that annotated the consequences of the variants: Ensembl VEP or bcftools csq") \
    OPTIONAL_STRING_ARG(csq_tag, "", "--csq-tag", "tag", "INFO field with the consequences (default CSQ for vep, BCSQ for bcftools)") \
    OPTIONAL_UINT_ARG(prediction_radius, 50, "--prediction-radius", "radius", "Prediction window radius in bases (non-default values give scores not comparable to SpliceAI)") \
    OPTIONAL_STRING_ARG(reference_scores, "", "--reference-scores", "file", "Precomputed reference scores (cpliceai_reference). Optional: without them, and for genes that are not in them like the annotation's, the reference predictions are made from the sequence while running")

#define BOOLEAN_ARGS \
    BOOLEAN_ARG(help, "-h", "Show help") \
    BOOLEAN_ARG(exclude_unphased, "--exclude-unphased", "Leave out heterozygous variants with an unknown phase (by default they are on both haplotypes)")

#include "easyargs.h"

static void init_logging(void) {
    log_set_level(LOG_INFO);
    /* TensorFlow logs independently of log.c; 1 drops its INFO output. */
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);
}

static int Session_open(const args_t *args, Session *s) {
    *s = (Session){0};

    if (args->prediction_radius < 1 || args->prediction_radius > MAX_PREDICTION_RADIUS) {
        log_error("--prediction-radius must be between 1 and %d, got %u", MAX_PREDICTION_RADIUS, args->prediction_radius);
        goto fail;
    }
    s->prediction_radius = (int) args->prediction_radius;

    AnnotationKind kind;
    if (strcmp(args->annotation, "vep") == 0) kind = ANNOTATION_VEP;
    else if (strcmp(args->annotation, "bcftools") == 0) kind = ANNOTATION_BCFTOOLS;
    else {
        log_error("--annotation must be vep or bcftools, got %s", args->annotation);
        goto fail;
    }

    const bool use_reference_scores = args->reference_scores[0] != '\0';

    log_info("Loading SpliceAI models from: %s", args->model_dir);
    if (!(s->models = load_models(args->model_dir))) {
        log_error("Could not load SpliceAI models from: %s", args->model_dir);
        goto fail;
    }

    /* gff_init and gff_parse report their own failures through bcftools' error(), which exits. */
    log_info("Loading GFF from: %s", args->gff);
    s->gff = gff_init(args->gff);
    gff_parse(s->gff);

    log_info("Loading FASTA from: %s", args->fasta);
    if (!(s->fai = fai_load(args->fasta))) {
        log_error("Could not load FASTA index for: %s", args->fasta);
        goto fail;
    }

    /* If they are given but cannot be used (another genome than they were made from), the run stops: predicting the reference from
     * the sequence instead would be quietly doing something else than was asked. Leaving the option out does that. A gene that is not
     * in them like the annotation's is predicted from the sequence, with a warning. */
    if (use_reference_scores) {
        log_info("Loading reference scores from: %s", args->reference_scores);
        if (!(s->reference = ReferenceScores_open(args->reference_scores, s->fai))) goto fail;
    } else {
        log_info("No reference scores given: the reference predictions are made from the sequence");
    }

    log_info("Loading variants from: %s", args->variants);
    if (!(s->fp_in = hts_open(args->variants, "r"))) {
        log_error("Could not open variants file: %s", args->variants);
        goto fail;
    }
    if (!(s->hdr = bcf_hdr_read(s->fp_in))) {
        log_error("Could not read VCF header from: %s", args->variants);
        goto fail;
    }
    if (!(s->vep_fmt = parse_csq_hdr_fmt(s->hdr, kind, args->csq_tag))) {
        log_error("Could not read the consequence annotation of %s (%s); annotate with VEP or bcftools csq first.",
                  args->variants, args->annotation);
        goto fail;
    }

    log_info("Initialising annotated VCF as %s", args->variants_output);
    if (!(s->fp_out = hts_open(args->variants_output, "w"))) {
        log_error("Could not open annotated VCF for writing: %s", args->variants_output);
        goto fail;
    }
    if (add_rescueranger_tags(s->hdr) != 0) goto fail;
    if (bcf_hdr_write(s->fp_out, s->hdr) != 0) {
        log_error("Could not write VCF header to: %s", args->variants_output);
        goto fail;
    }

    log_info("Initialising raw scores TSV as %s", args->scores_output);
    if (!(s->scores_tsv = fopen(args->scores_output, "w"))) {
        log_error("Could not open raw scores TSV for writing: %s", args->scores_output);
        goto fail;
    }
    if (fprintf(s->scores_tsv,
                "CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t"
                "PHASE\tOFFSET\tSOURCE\tTYPE\tDONOR\tACCEPTOR\tNEITHER\n") < 0) {
        log_error("Could not write header to raw scores TSV: %s", args->scores_output);
        goto fail;
    }

    return EXIT_SUCCESS;

fail:
    Session_close(s);
    return EXIT_FAILURE;
}

// Writes the records to the output. Returns 0, or -1 after logging the first record that could not be written.
static int write_records(const Session *s, const Records *records) {
    for (size_t i = 0; i < kv_size(*records); i++) {
        bcf1_t *v = kv_A(*records, i)->v;
        if (bcf_write(s->fp_out, s->hdr, v) != 0) {
            log_error("Could not write the record on %s:%li to the output, stopping", bcf_hdr_id2name(s->hdr, v->rid), (long) v->pos + 1);
            return -1;
        }
    }

    return 0;
}

int main(int argc, char *argv[]) {
    init_logging();

    // Parse arguments
    args_t args = make_default_args();
    if (!parse_args(argc, argv, &args) || args.help) {
        print_help(argv[0]);
        return EXIT_FAILURE;
    }

    Session session;
    if (Session_open(&args, &session) != EXIT_SUCCESS) return EXIT_FAILURE;

    RecordBuffer buffer = { .v = bcf_init(), .unused = false, .failed = false };
    RecordBlock block = { 0 };

    int status;
    while((status = RecordBlock_next(&session, &buffer, &block)) == 0) {
        // If there are no coding variants in the entire block, there's nothing to do
        if (!block.is_coding) {
            const int written = write_records(&session, &block.records);
            RecordBlock_destroy(block);
            if (written != 0) { status = 2; break; }
            continue;
        }

        // Split on phase
        // Every phase gets its own copies of the transcripts, so that each has its own splice sites for isoforms and flags
        for (int phase = 0; phase < NUM_PHASES; phase++) {
            RecordBlock_begin_phase(&block);

            // Get records for phase. Variants that conflict with an earlier variant are removed first, so that they
            // do not get triggers and do not take part in a reading frame conflict.
            Records records = Record_select_for_phase(&block.records, session.hdr, phase, !args.exclude_unphased);
            filter_overlapping_variants(&records, session.hdr);

            Triggers_compute_for_phase(&records, session.hdr, args.af_tag, args.af_threshold, phase);

            init_splice_sites_for_phase(&records);

            predict_splice_changes_for_phase(&session, phase, &records);

            build_isoforms_for_phase(&session, phase, &records);

            // Clean up records
            kv_destroy(records);

            // The copies of this phase are destroyed together with their splice sites. The genes and the transcripts
            // they are copied from belong to the annotation of the block.
            RecordBlock_end_phase(&block);
        }

        write_tags(&session, &block.records);

        // Write out records
        const int written = write_records(&session, &block.records);
        RecordBlock_destroy(block);
        if (written != 0) { status = 2; break; }
    }

    bcf_destroy(buffer.v);
    Session_close(&session);

    return status == 1 ? EXIT_SUCCESS : EXIT_FAILURE; // 1 is the end of the input, anything else is an error
}
