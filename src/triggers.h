#ifndef TRIGGERS_HEADER
#define TRIGGERS_HEADER

#include <stdbool.h>
#include <htslib/vcf.h>

#include "record.h"

// Indices into TscriptCsqs.triggers
#define SPLICE_TRIGGER 0
#define FRAMESHIFT_TRIGGER 1
#define STOP_TRIGGER 2

#define CANONICAL_SPLICE_VALUE 1
#define AF_SPLICE_VALUE 2

#define STOP_RESCUE_VALUE 1
#define STOP_GAIN_VALUE 2

// VEP consequences that set a trigger
#define SPLICE_DONOR_VARIANT "splice_donor_variant"
#define SPLICE_ACCEPTOR_VARIANT "splice_acceptor_variant"
#define SPLICE_REGION_VARIANT "splice_region_variant"
#define FRAMESHIFT_VARIANT "frameshift_variant"
#define STOP_GAINED_VARIANT "stop_gained"

// Whether INFO/<af_tag> of the record is present and at most af_threshold
int has_af_trigger(const bcf1_t *record, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold);

// Sets the triggers of the VEP consequences that belong to the transcript, from the consequence names
void set_vep_triggers(Records records, int transcript_id);

// Sets the splice trigger of every consequence of a record that is rare enough, unless it already has one
void set_af_triggers(Records records, const bcf_hdr_t *hdr, const char *af_tag, const float af_threshold);

// The reading frame trigger: when two records are at most two coding bases apart, together they can change a codon that VEP
// checks per variant. Both records get reading_frame[phase] on their consequence for the transcript. The records must be sorted
// by position. Whether two records conflict depends on who is on the haplotype, so the trigger is kept per phase.
void set_mnv_triggers(Records records, const Transcript *transcript, int phase);

// The value of a trigger of a consequence for a phase: the value from the consequences and the allele frequency, and for
// STOP_TRIGGER, when there is none, STOP_GAIN_VALUE if the consequence has a reading frame conflict in the phase. A stop that
// VEP found (STOP_RESCUE_VALUE) comes first. The code that reads the triggers of a phase goes through this.
int Trigger_value(const TscriptCsqs *csq, int trigger, int phase);

// Computes the triggers for the records of one phase. phase_records are the records that are on the haplotype of the phase
// and have coding transcripts (see Record_select_for_phase), after the records that conflict with an earlier variant have
// been removed, so a removed record neither gets triggers nor takes part in a reading frame conflict. For every transcript
// this sets the VEP, allele frequency and reading frame triggers of the records that have it. The reading frame triggers of
// these records for this phase are cleared first, so computing a phase again starts from nothing; other phases are not touched.
// Every record must hold its copies of the transcripts for this phase.
void Triggers_compute_for_phase(const Records *phase_records, const bcf_hdr_t *hdr, const char *af_tag, float af_threshold, int phase);

#endif
