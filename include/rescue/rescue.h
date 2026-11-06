/*  rescue.h -- predict the rescue of protein-truncating variants.

    The public interface of librescue. The variants of a phased VCF that
    were annotated with the consequences of Ensembl VEP or bcftools csq are
    checked for splice changes with SpliceAI-style models (cpliceai), per
    haplotype. A protein-truncating variant (stop gained, frameshift, splice
    donor or acceptor) that is rescued by another variant of the haplotype, or
    that induces an exon skip or an intron retention, is tagged, together with
    whether the transcript is predicted to be degraded by nonsense-mediated
    decay.

    Typical use, with the error handling left out:

        RescueOptions opts;
        rr_opts_init(&opts);
        opts.model_dir = "models/"; opts.gff_fname = "annot.gff3.gz"; opts.fasta_fname = "ref.fa";

        htsFile *in = hts_open("in.vcf.gz", "r");
        bcf_hdr_t *hdr = bcf_hdr_read(in);
        Rescue *rr = rr_init(&opts, hdr);         // adds the output INFO tags to hdr

        htsFile *out = hts_open("-", "w");
        bcf_hdr_write(out, hdr);
        int ret = rr_annotate(rr, in, out);     // the records after the header

        rr_destroy(rr);
        hts_close(out); bcf_hdr_destroy(hdr); hts_close(in);

    Only the functions and types in this file are exported by the library. */

#ifndef RESCUE_H
#define RESCUE_H

#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
#  define RR_EXPORT
#else
#  define RR_EXPORT __attribute__((visibility("default")))
#endif

#include <rescue/version.h> // RR_VERSION, from the version in CMakeLists.txt

// The version of the library that is linked, which can differ from RR_VERSION of the header that was compiled against
RR_EXPORT const char *rr_version(void);

// Which program wrote the consequence annotation of the input
typedef enum RescueAnnotation {
    RR_ANNOTATION_VEP = 0,      // Ensembl VEP: INFO/CSQ
    RR_ANNOTATION_BCFTOOLS = 1, // bcftools csq: INFO/BCSQ
} RescueAnnotation;

// Parses "vep" or "bcftools". Returns 0, or -1 for anything else.
RR_EXPORT int rr_annotation_parse(const char *name, RescueAnnotation *annotation);

// The options of a run. Strings are copied by rr_init, so they only have to live until it returns.
typedef struct RescueOptions {
    const char *model_dir;          // Directory with the SpliceAI models. Required
    const char *gff_fname;          // GFF3 annotation matching the VCF, in the format of bcftools csq. Required
    const char *fasta_fname;        // Indexed reference sequence. Required
    const char *reference_scores;   // Precomputed reference scores (cpliceai_reference), or NULL to predict the reference from the sequence
    RescueAnnotation annotation;     // [RR_ANNOTATION_VEP]
    const char *csq_tag;            // INFO field with the consequences, NULL or "" for CSQ (VEP) or BCSQ (bcftools) [NULL]
    const char *af_tag;             // Allele frequency tag ["INFO/AF"]
    const char *gt_tag;             // Phased genotype tag, a FORMAT field: GT, or a String field written like GT ("0|1") ["FORMAT/GT"]
    double af_threshold;            // Allele frequency threshold, 0 to 1 [0.01]
    int prediction_radius;          // Prediction window radius in bases; other values than 50 give scores not comparable to SpliceAI [50]
    int exclude_unphased;           // Leave out heterozygous variants with an unknown phase; by default they are on both haplotypes [0]
} RescueOptions;

// Sets the defaults given above. The required fields are set to NULL.
RR_EXPORT void rr_opts_init(RescueOptions *opts);

typedef struct Rescue Rescue;

// Loads the models, the annotation, the reference and the reference scores, and reads the consequence annotation from the
// header of the input. The INFO tags that are written to the records are added to hdr, which is then the header to write the
// output with. hdr stays the caller's and must stay valid until rr_destroy. Returns NULL after logging why not.
//
// The GFF is read with the GFF reader of bcftools, which ends the program, as bcftools does, when the GFF cannot be parsed.
RR_EXPORT Rescue *rr_init(const RescueOptions *opts, bcf_hdr_t *hdr);

// Writes the raw SpliceAI scores to fp as TSV, starting with its header line. fp stays the caller's and must stay open until
// rr_destroy. Without it no raw scores are written. Returns 0, or -1 if the header line could not be written.
RR_EXPORT int rr_set_scores_output(Rescue *rr, FILE *fp);

// Reads the records from in, which must be positioned after the header given to rr_init, and writes every one of them to out
// with the tags added, in the input order. The header must have been written to out. Records with more than one ALT allele
// are passed through untouched; split them first, for example with bcftools norm -m-.
// Returns 0 at the end of the input, or -1 when the input could not be read, an annotation could not be built or a record could
// not be written (the records before it are written).
RR_EXPORT int rr_annotate(Rescue *rr, htsFile *in, htsFile *out);

// Frees what rr_init loaded. The header, the files and the scores output are not closed. NULL is allowed.
RR_EXPORT void rr_destroy(Rescue *rr);

// What the library logs to stderr. Messages below the level are left out [RR_LOG_TRACE: everything]
enum { RR_LOG_TRACE, RR_LOG_DEBUG, RR_LOG_INFO, RR_LOG_WARN, RR_LOG_ERROR, RR_LOG_FATAL };
RR_EXPORT void rr_set_log_level(int level);

#ifdef __cplusplus
}
#endif

#endif
