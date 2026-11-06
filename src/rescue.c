// The stand-alone command line program. Everything it does goes through the public interface of the library, like the bcftools
// plugin does (bcftools +rescue).

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <rescue/rescue.h>

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
    rr_set_log_level(RR_LOG_INFO);
    /* TensorFlow logs independently of the library; 1 drops its INFO output. */
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);
}

int main(int argc, char *argv[]) {
    init_logging();

    // Parse arguments
    args_t args = make_default_args();
    if (!parse_args(argc, argv, &args) || args.help) {
        print_help(argv[0]);
        return EXIT_FAILURE;
    }

    RescueOptions opts;
    rr_opts_init(&opts);
    opts.model_dir = args.model_dir;
    opts.gff_fname = args.gff;
    opts.fasta_fname = args.fasta;
    opts.reference_scores = args.reference_scores;
    opts.csq_tag = args.csq_tag;
    opts.af_tag = args.af_tag;
    opts.gt_tag = args.gt_tag;
    opts.af_threshold = args.af_threshold;
    opts.prediction_radius = (int) args.prediction_radius;
    opts.exclude_unphased = args.exclude_unphased;
    if (rr_annotation_parse(args.annotation, &opts.annotation) != 0) {
        fprintf(stderr, "--annotation must be vep or bcftools, got %s\n", args.annotation);
        return EXIT_FAILURE;
    }

    int ret = EXIT_FAILURE;
    htsFile *fp_in = NULL, *fp_out = NULL;
    bcf_hdr_t *hdr = NULL;
    FILE *scores_tsv = NULL;
    Rescue *rr = NULL;

    if (!(fp_in = hts_open(args.variants, "r"))) {
        fprintf(stderr, "Could not open variants file: %s\n", args.variants);
        goto done;
    }
    if (!(hdr = bcf_hdr_read(fp_in))) {
        fprintf(stderr, "Could not read VCF header from: %s\n", args.variants);
        goto done;
    }
    if (!(rr = rr_init(&opts, hdr))) goto done;

    if (!(fp_out = hts_open(args.variants_output, "w"))) {
        fprintf(stderr, "Could not open annotated VCF for writing: %s\n", args.variants_output);
        goto done;
    }
    if (bcf_hdr_write(fp_out, hdr) != 0) {
        fprintf(stderr, "Could not write VCF header to: %s\n", args.variants_output);
        goto done;
    }

    if (!(scores_tsv = fopen(args.scores_output, "w"))) {
        fprintf(stderr, "Could not open raw scores TSV for writing: %s\n", args.scores_output);
        goto done;
    }
    if (rr_set_scores_output(rr, scores_tsv) != 0) goto done;

    if (rr_annotate(rr, fp_in, fp_out) == 0) ret = EXIT_SUCCESS;

done:
    rr_destroy(rr);
    if (scores_tsv && fclose(scores_tsv) != 0) ret = EXIT_FAILURE;
    if (fp_out && hts_close(fp_out) != 0) ret = EXIT_FAILURE;
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp_in) hts_close(fp_in);

    return ret;
}
