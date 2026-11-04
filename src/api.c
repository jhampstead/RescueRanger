// The public interface of the library (include/rescueranger/rescueranger.h): the options, opening a Session, and the order in
// which the stages run on the blocks of records.

#include <rescueranger/rescueranger.h>

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
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

struct RescueRanger {
    Session s; // fp_in, fp_out, hdr and scores_tsv are the caller's, see rr_destroy
    char *af_tag;
    float af_threshold;
    bool include_unphased;
};

const char *rr_version(void) {
    return RR_VERSION;
}

int rr_annotation_parse(const char *name, RescueRangerAnnotation *annotation) {
    if (name == NULL) return -1;
    if (strcmp(name, "vep") == 0) { *annotation = RR_ANNOTATION_VEP; return 0; }
    if (strcmp(name, "bcftools") == 0) { *annotation = RR_ANNOTATION_BCFTOOLS; return 0; }
    return -1;
}

void rr_opts_init(RescueRangerOptions *opts) {
    *opts = (RescueRangerOptions){
        .model_dir = NULL,
        .gff_fname = NULL,
        .fasta_fname = NULL,
        .reference_scores = NULL,
        .annotation = RR_ANNOTATION_VEP,
        .csq_tag = NULL,
        .af_tag = "INFO/AF",
        .gt_tag = "FORMAT/GT",
        .af_threshold = 0.01,
        .prediction_radius = 50,
        .exclude_unphased = 0,
    };
}

void rr_set_log_level(int level) {
    log_set_level(level);
}

void rr_destroy(RescueRanger *rr) {
    if (rr == NULL) return;

    // Session_close also closes the files and the header, which belong to the caller here
    rr->s.fp_in = NULL;
    rr->s.fp_out = NULL;
    rr->s.hdr = NULL;
    rr->s.scores_tsv = NULL;
    Session_close(&rr->s);

    free(rr->af_tag);
    free(rr);
}

RescueRanger *rr_init(const RescueRangerOptions *opts, bcf_hdr_t *hdr) {
    if (opts->af_threshold < 0 || opts->af_threshold > 1) {
        log_error("--af-threshold must be between 0 and 1, got %f", opts->af_threshold);
        return NULL;
    }

    RescueRanger *rr = calloc(1, sizeof(RescueRanger));
    if (rr == NULL) {
        log_error("Out of memory");
        return NULL;
    }
    rr->af_threshold = (float) opts->af_threshold;
    rr->include_unphased = !opts->exclude_unphased;
    if (!(rr->af_tag = strdup(opts->af_tag ? opts->af_tag : "INFO/AF"))) {
        log_error("Out of memory");
        goto fail;
    }

    if (Session_open(opts, hdr, &rr->s) != EXIT_SUCCESS) goto fail;

    return rr;

fail:
    rr_destroy(rr);
    return NULL;
}

int rr_set_scores_output(RescueRanger *rr, FILE *fp) {
    if (fprintf(fp,
                "CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t"
                "PHASE\tOFFSET\tSOURCE\tTYPE\tDONOR\tACCEPTOR\tNEITHER\n") < 0) {
        log_error("Could not write the header of the raw scores");
        return -1;
    }
    rr->s.scores_tsv = fp;
    return 0;
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

// Runs the stages on the blocks of records, one block at a time. A block is a run of records that share a gene and so
// a sequence to predict on, see record_block.h.
int rr_annotate(RescueRanger *rr, htsFile *in, htsFile *out) {
    Session *s = &rr->s;
    s->fp_in = in;
    s->fp_out = out;

    RecordBuffer buffer = { .v = bcf_init(), .unused = false, .failed = false };
    RecordBlock block = { 0 };

    int status;
    while ((status = RecordBlock_next(s, &buffer, &block)) == 0) {
        // If there are no coding variants in the entire block, there's nothing to do
        if (!block.is_coding) {
            const int written = write_records(s, &block.records);
            RecordBlock_destroy(block);
            if (written != 0) { status = 2; break; }
            continue;
        }

        // Split on phase. Every phase gets its own copies of the transcripts, so that each has its own splice sites for isoforms and flags
        for (int phase = 0; phase < NUM_PHASES; phase++) {
            RecordBlock_begin_phase(&block);

            // Get records for phase. Variants that conflict with an earlier variant are removed first, so that they
            // do not get triggers and do not take part in a reading frame conflict.
            Records records = Record_select_for_phase(&block.records, s->hdr, s->gt_tag, phase, rr->include_unphased);
            filter_overlapping_variants(&records, s->hdr);

            Triggers_compute_for_phase(&records, s->hdr, rr->af_tag, rr->af_threshold, phase);
            init_splice_sites_for_phase(&records);
            predict_splice_changes_for_phase(s, phase, &records);
            build_isoforms_for_phase(s, phase, &records);

            kv_destroy(records);

            // The copies of this phase are destroyed together with their splice sites. The genes and the transcripts
            // they are copied from belong to the annotation of the block.
            RecordBlock_end_phase(&block);
        }

        write_tags(s, &block.records);

        const int written = write_records(s, &block.records);
        RecordBlock_destroy(block);
        if (written != 0) { status = 2; break; }
    }
    bcf_destroy(buffer.v);

    s->fp_in = NULL;
    s->fp_out = NULL;

    return status == 1 ? 0 : -1; // 1 is the end of the input, anything else is an error
}
