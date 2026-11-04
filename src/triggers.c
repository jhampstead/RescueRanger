#include "triggers.h"

#include <stdlib.h>
#include <string.h>

#include "klib/kvec.h"
#include "utils.h"

int has_af_trigger(const bcf1_t *record, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold) {
    float *af_value = NULL;
    int ndst = 0;

    // Set based on INFO/AF. The buffer is allocated by htslib, also when there is no value.
    const int n_values = bcf_get_info_float(hdr, (bcf1_t *) record, af_tag, &af_value, &ndst);
    const int triggered = n_values > 0 && *af_value <= af_threshold;

    free(af_value);

    return triggered;
}

void set_vep_triggers(Records records, int transcript_id) {
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);

        for (int vi = 0; vi < r->vep->n; vi++) {
            TscriptCsqs *transcript_consequence = kv_A(*r->vep, vi);
            Transcript *transcript = kv_A(*r->tscripts, vi);

            if (transcript == NULL || transcript->id != transcript_id) continue;

            for (int ci = 0; ci < transcript_consequence->num_consequences; ci++) {
                const char *consequence = transcript_consequence->consequences[ci];
                if ((strcmp(consequence, "non_coding_transcript_variant") == 0) | (strcmp(consequence, "NMD_transcript_variant") == 0)) continue;
                if ((strcmp(consequence, SPLICE_DONOR_VARIANT) == 0) | (strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0)) transcript_consequence->triggers[SPLICE_TRIGGER] = CANONICAL_SPLICE_VALUE;
                if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) transcript_consequence->triggers[FRAMESHIFT_TRIGGER] = 1;
                if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) transcript_consequence->triggers[STOP_TRIGGER] = STOP_RESCUE_VALUE;
            }
        }
    }
}

void set_af_triggers(Records records, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold) {
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);
        int triggered = has_af_trigger(r->v, (bcf_hdr_t *) hdr, af_tag, af_threshold);

        if (!triggered) continue;

        for (int vi = 0; vi < r->vep->n; vi++) {
            TscriptCsqs *transcript_consequences = kv_A(*r->vep, vi);

            if (transcript_consequences->triggers[SPLICE_TRIGGER]) continue;

            transcript_consequences->triggers[SPLICE_TRIGGER] = AF_SPLICE_VALUE;
        }
    }
}

int Trigger_value(const TscriptCsqs *csq, const int trigger, const int phase) {
    if (csq->triggers[trigger]) return csq->triggers[trigger];
    if (trigger == STOP_TRIGGER && csq->reading_frame[phase]) return STOP_GAIN_VALUE;

    return 0;
}

// Marks the consequences of the record that belong to the transcript with a reading frame conflict in the phase
static void set_reading_frame_trigger(Record *r, const uint32_t transcript_id, const int phase) {
    for (int ti = 0; ti < r->tscripts->n; ti++) {
        const Transcript *record_transcript = kv_A(*r->tscripts, ti);
        if (record_transcript == NULL || record_transcript->id != transcript_id) continue;

        kv_A(*r->vep, ti)->reading_frame[phase] = 1;
    }
}

// TODO: Rename this to something other than MNV, as it really just checks for possible reading frame conflicts not checked by VEP
void set_mnv_triggers(Records records, const Transcript *transcript, const int phase) {
    size_t n_exons = kv_size(transcript->exons);

    // Figure out the positions in the coding transcript
    int64_t transcript_start = transcript->start;
    int64_t positions[records.n];
    int position_record_map[records.n];
    size_t n_positions = 0;
    for (int ri = 0; ri < records.n; ri++) {
        Record *r = kv_A(records, ri);
        int64_t record_pos = variant_position(&r->alleles); // A deletion counts at its first deleted base, an insertion at the base it is attached to

        bool is_coding = true;
        int skipped_intron_bases = 0;
        for (int ei = 0; ei < n_exons; ei++) {
            Range exon = *kv_A(transcript->exons, ei)->ref_bounds;

            // Either intron or UTR, neither is coding
            if (record_pos < exon.start) {
                is_coding = false;
                break;
            }

            if (record_pos >= exon.end) {
                // UTR if beyond last exon
                if (ei == n_exons-1) {
                    is_coding = false;
                    break;
                }

                skipped_intron_bases += kv_A(transcript->introns, ei).len;
                continue;
            }

            break; // Falls within current exon, we're done!
        }

        if (!is_coding) continue;

        position_record_map[n_positions] = ri;
        positions[n_positions++] = record_pos - transcript_start - skipped_intron_bases;
    }

    if (!n_positions) return; // No coding positions, nothing to be done!

    for (int pi = 0; pi < n_positions-1; pi++) { // n_positions is size_t, so without the if size 0 -> return this can overflow
        int64_t curr = positions[pi];
        int64_t next = positions[pi+1];

        // Cannot be inframe!
        if (next - curr > 2) continue;

        // Both records of the pair, for their consequences for this transcript
        set_reading_frame_trigger(kv_A(records, position_record_map[pi]), transcript->id, phase);
        set_reading_frame_trigger(kv_A(records, position_record_map[pi+1]), transcript->id, phase);
    }
}

void Triggers_compute_for_phase(const Records *phase_records, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold, const int phase) {
    // Start the phase from nothing: the conflicts of an earlier computation of the phase may not be there any more
    for (int i = 0; i < phase_records->n; i++) {
        Record *r = kv_A(*phase_records, i);
        if (r->tscripts == NULL) continue;

        for (int vi = 0; vi < r->vep->n; vi++) kv_A(*r->vep, vi)->reading_frame[phase] = 0;
    }

    // The first copy of each transcript, used to look at the exons of each transcript once. It does not own the transcripts.
    khash_t(TranscriptHash) *transcripts = kh_init_TranscriptHash();
    for (int i = 0; i < phase_records->n; i++) {
        Record *r = kv_A(*phase_records, i);
        if (r->tscripts == NULL) continue;

        for (int ti = 0; ti < r->tscripts->n; ti++) {
            Transcript *t = kv_A(*r->tscripts, ti);
            if (t == NULL) continue;

            int is_new;
            khint32_t hash_index = kh_put_TranscriptHash(transcripts, t->id, &is_new);
            if (!is_new) continue;
            kh_val(transcripts, hash_index) = t;
        }
    }

    // Triggers need to checked per transcript
    for (khiter_t k = kh_begin(transcripts); k != kh_end(transcripts); ++k) {
        if (!kh_exist(transcripts, k)) continue;

        Transcript *transcript = kh_val(transcripts, k);

        // Find all records for current transcript
        Records records;
        kv_init(records);
        for (int ri = 0; ri < phase_records->n; ri++) {
            Record *r = kv_A(*phase_records, ri);

            if (r->tscripts == NULL) continue;

            for (int ti = 0; ti < r->tscripts->n; ti++) {
                Transcript *record_transcript = kv_A(*r->tscripts, ti);
                if (record_transcript == NULL || record_transcript->id != transcript->id) continue;

                kv_push(Record *, records, r);
            }
        }

        // TODO: Separate rescue and induced ptv triggers to get multi purpose functions for better readability
        set_vep_triggers(records, transcript->id);
        set_af_triggers(records, hdr, af_tag, af_threshold);
        set_mnv_triggers(records, transcript, phase);

        kv_destroy(records);
    }

    kh_destroy_TranscriptHash(transcripts);
}
