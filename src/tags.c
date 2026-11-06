#include "tags.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging/log.h"
#include "triggers.h"
#include "utils.h"

#define TAG_SIZE 10
int add_rescue_tags(bcf_hdr_t *hdr) {
    const char *tags[TAG_SIZE] = INFO_TAGS;
    const char *descriptions[TAG_SIZE] = INFO_TAG_DESCS;

    for (int i = 0; i < TAG_SIZE; i++) {
        if (bcf_hdr_append(hdr, descriptions[i]) != 0) {
            log_error("Failed to append description for tag %s to vcf header.", tags[i]);
            return -1;
        }
    }

    return 0;
}
#undef TAG_SIZE

kstring_t generate_rescue_type_tag(Isoforms *isoforms_struct) { // INFO/RESCUE_TYPE rescue_type|probability|affected_bases
    kstring_t tag = {0};

    for (int i = 0; i < isoforms_struct->n_rescue_types; i++) {
        if (tag.l > 0) kputc(':', &tag);

        kputw(isoforms_struct->rescue_type_probabilities->rescues[i], &tag);
        kputc('|', &tag);

        enum RescueType rescue_type = isoforms_struct->rescue_type_probabilities->rescue_types[i];
        switch(rescue_type) {
            case COMPLETE_RESCUE:
                kputs("complete_rescue", &tag);
                break;
            case POSSIBLE_REINITIATION_SITE:
                kputs("possible_reinitiation", &tag);
                break;
            case FRAME_RESTORING_INDEL:
                kputs("frame_restoring_indel", &tag);
                break;
            case PARTIAL_EXON_SKIPPING:
                kputs("partial_exon_skipping", &tag);
                break;
            case PARTIAL_INTRON_RETENTION:
                kputs("partial_intron_retention", &tag);
                break;
            case EXON_SKIPPING:
                kputs("exon_skipping", &tag);
                break;
            case INTRON_RETENTION:
                kputs("intron_retention", &tag);
                break;
            case NO_RESCUE:
                kputs("no_rescue", &tag);
                break;
        }

        kputc('|', &tag);
        char prob_s[7];
        sprintf(prob_s, "%.4f", isoforms_struct->rescue_type_probabilities->rescue_probabilities[i]);
        kputsn(prob_s, 6, &tag);
        kputc('|', &tag);
        kputw(isoforms_struct->rescue_type_probabilities->affected_bases[i], &tag);

    }

    return tag;
}

void get_prediction_flags(uint32_t flags, kstring_t *concatenated_flags) {
    const char *flag_names[] = {
        "unphased_variant_in_prediction",
        "skipped_variant_in_prediction",
        "variant_larger_than_splice_window",
        "rescued_by_splice_frameshift",
        "possible_reinitiation",
        "manufactured_splice_rescue",
        "haplotype_specific_donor",
        "haplotype_specific_acceptor",
        "cryptic_donor",
        "cryptic_acceptor",
        "weak_canonical_donor",
        "weak_canonical_acceptor",
        "ablated_canonical_donor",
        "ablated_canonical_acceptor",
        "ablated_cryptic_donor",
        "ablated_cryptic_acceptor",
        "upstream_ptv",
        "5utr_variant",
        "3utr_variant",
        "last_exon",
        "single_exon",
        "non_constitutive",
        "in_frame",
        "start_proximal",
        "long_exon",
        "distal"
    };

    size_t flag_count = sizeof(flag_names) / sizeof(flag_names[0]);
    int is_first = 1;

    for (size_t i = 0; i < flag_count; i++) {
        if (flags & (1 << i)) {
            if (!is_first) {
                kputs(":", concatenated_flags); // append delimiter if not the first flag
            }
            kputs(flag_names[i], concatenated_flags); // append the flag name
            is_first = 0;
        }
    }
}

static bool has_splice_trigger(const Record *r) {
    for (int j = 0; j < r->vep->n; j++) {
        if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) return true;
    }
    return false;
}

// Whether a consequence of the record has a trigger, in any phase
static bool has_any_trigger(const Record *r) {
    for (int j = 0; j < r->vep->n; j++) {
        TscriptCsqs *csq = kv_A(*r->vep, j);
        if (csq->triggers[SPLICE_TRIGGER] || csq->triggers[FRAMESHIFT_TRIGGER]) return true;

        for (int p = 0; p < NUM_PHASES; p++) {
            if (Trigger_value(csq, STOP_TRIGGER, p)) return true;
        }
    }
    return false;
}

void format_gene_scores(const Record *r, Score **scores, kstring_t *out) {
    for (int g = 0; g < r->genes->n; g++) {
        if (g > 0) kputc(',', out);

        for (int p = 0; p < NUM_PHASES; p++) {
            if (p > 0) kputc('&', out);

            Score *score = scores[g * NUM_PHASES + p];
            kputs(r->v->d.allele[1], out);
            if (score == NULL) {
                kputc('|', out); kputs(kv_A(*r->genes, g)->name, out);
                kputs("|.|.|.|.|.|.|.|.", out);
                continue;
            }
            kputc('|', out); kputs(score->gene, out);

            char ds_scores[24];
            snprintf(ds_scores, sizeof(ds_scores), "%.2f|%.2f|%.2f|%.2f", score->ag, score->al, score->dg, score->dl);
            kputc('|', out); kputs(ds_scores, out);

            kputc('|', out); kputw(score->ag_idx, out);
            kputc('|', out); kputw(score->al_idx, out);
            kputc('|', out); kputw(score->dg_idx, out);
            kputc('|', out); kputw(score->dl_idx, out);
        }
    }
}

static void write_score_tag(const Session *s, const Record *r, Score **scores, const char *tag) {
    if (!has_splice_trigger(r)) return;
    if (!variant_fits_window(&r->alleles, s->prediction_radius)) return;

    kstring_t value = {0};
    format_gene_scores(r, scores, &value);
    bcf_update_info_string(s->hdr, r->v, tag, value.s);
    free(value.s);
}

static void write_rescue_tags(const Session *s, const Record *r) {
    int triggered = 0;
    int has_ptv_effect = 0;
    for (int j = 0; j < r->vep->n; j++) {
        TscriptCsqs *csq = kv_A(*r->vep, j);
        if (csq->triggers[SPLICE_TRIGGER] == CANONICAL_SPLICE_VALUE || csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE) triggered = 1;
        if (r->isoforms->valid && r->isoforms->ptv_probability > 0) has_ptv_effect = 1;
    }

    if (!triggered && !has_ptv_effect) return;
    if (!variant_fits_window(&r->alleles, s->prediction_radius)) return;

    kstring_t rescue_tag = {0};
    kstring_t rescue_prob_tag = {0};
    kstring_t rescue_type_tag = {0};

    for (int ti = 0; ti < r->transcript_genes->n; ti++) {
        if (ti > 0) {
            kputc(',', &rescue_tag);
            kputc(',', &rescue_prob_tag);
            kputc(',', &rescue_type_tag);
        }

        if (kv_A(*r->transcript_genes, ti) == NULL) {
            kputc('.', &rescue_tag);
            kputc('.', &rescue_prob_tag);
            kputc('.', &rescue_type_tag);
            continue;
        }

        TscriptCsqs *csq = kv_A(*r->vep, ti);

        // No Rescue triggers at all
        if ((csq->triggers[SPLICE_TRIGGER] != CANONICAL_SPLICE_VALUE) & (!csq->triggers[FRAMESHIFT_TRIGGER]) & (csq->triggers[STOP_TRIGGER] != STOP_RESCUE_VALUE)) {
            kputc('.', &rescue_tag);
            kputc('.', &rescue_prob_tag);
            kputc('.', &rescue_type_tag);
            continue;
        }

        for (int p = 0; p < NUM_PHASES; p++) {
            if (p > 0) {
                kputc('&', &rescue_tag);
                kputc('&', &rescue_prob_tag);
                kputc('&', &rescue_type_tag);
            }

            Isoforms isoforms = r->isoforms[ti * NUM_PHASES + p];
            if (!isoforms.valid | (
                (((Trigger_value(csq, STOP_TRIGGER, p) == STOP_GAIN_VALUE) | (csq->triggers[SPLICE_TRIGGER] == AF_SPLICE_VALUE)) & (isoforms.rescue == COMPLETE_RESCUE)) &
                !((csq->triggers[SPLICE_TRIGGER] == CANONICAL_SPLICE_VALUE) | (csq->triggers[FRAMESHIFT_TRIGGER]) | (csq->triggers[STOP_TRIGGER] == STOP_RESCUE_VALUE)))
            ) {
                kputc('.', &rescue_tag);
                kputc('.', &rescue_prob_tag);
                kputs(".|.|.|.", &rescue_type_tag);
                continue;
            }

            kputw(isoforms.rescue, &rescue_tag);

            char prob_s[7];
            sprintf(prob_s, "%.4f", isoforms.rescue_probability);
            kputsn(prob_s, 6, &rescue_prob_tag);

            kstring_t transcript_phase_tag = generate_rescue_type_tag(&isoforms);
            kputsn(transcript_phase_tag.s, transcript_phase_tag.l, &rescue_type_tag);
            free(transcript_phase_tag.s);
        }
    }

    bcf_update_info_string(s->hdr, r->v, RESCUE_TAG, rescue_tag.s);
    bcf_update_info_string(s->hdr, r->v, RESCUE_PROB_TAG, rescue_prob_tag.s);
    bcf_update_info_string(s->hdr, r->v, RESCUE_TYPE_TAG, rescue_type_tag.s);
    free(rescue_tag.s);
    free(rescue_prob_tag.s);
    free(rescue_type_tag.s);
}

static void write_ptv_tags(const Session *s, const Record *r) {
    if (!has_any_trigger(r)) return;
    if (!variant_fits_window(&r->alleles, s->prediction_radius)) return;

    kstring_t ptv_tag = {0};
    kstring_t ptv_prob_tag = {0};
    kstring_t ptv_type_tag = {0};

    for (int ti = 0; ti < r->transcript_genes->n; ti++) {
        if (ti > 0) {
            kputc(',', &ptv_tag);
            kputc(',', &ptv_prob_tag);
            kputc(',', &ptv_type_tag);
        }

        for (int p = 0; p < NUM_PHASES; p++) {
            if (p > 0) {
                kputc('&', &ptv_tag);
                kputc('&', &ptv_prob_tag);
                kputc('&', &ptv_type_tag);
            }

            Isoforms isoforms = r->isoforms[ti * NUM_PHASES + p];

            if (!isoforms.valid) {
                kputc('.', &ptv_tag);
                kputc('.', &ptv_prob_tag);
                kputc('.', &ptv_type_tag);
                continue;
            }

            kputw(isoforms.ptv, &ptv_tag);

            char prob_s[7];
            sprintf(prob_s, "%.4f", isoforms.ptv_probability);
            kputsn(prob_s, 6, &ptv_prob_tag);

            switch(isoforms.ptv_type) {
                case STOP_GAINED:
                    kputs("stop_gained", &ptv_type_tag);
                    break;
                case FRAMESHIFT:
                    kputs("frameshift", &ptv_type_tag);
                    break;
                case LONG_INFRAME_INDEL:
                    kputs("long_inframe_indel", &ptv_type_tag);
                    break;
                case NO_PTV:
                    kputs("no_ptv", &ptv_type_tag);
                    break;
            }
        }
    }

    bcf_update_info_string(s->hdr, r->v, PTV_TAG, ptv_tag.s);
    bcf_update_info_string(s->hdr, r->v, PTV_PROB_TAG, ptv_prob_tag.s);
    bcf_update_info_string(s->hdr, r->v, PTV_TYPE_TAG, ptv_type_tag.s);
    free(ptv_tag.s);
    free(ptv_prob_tag.s);
    free(ptv_type_tag.s);
}

static void write_prediction_flags_tag(const Session *s, const Record *r) {
    if (r->transcript_genes == NULL) return;
    if (!has_any_trigger(r)) return;
    if (!variant_fits_window(&r->alleles, s->prediction_radius)) return;

    kstring_t tags = {0};

    for (int t = 0; t < r->transcript_genes->n; t++) {
        if (t > 0) kputc(',', &tags);

        for (int p = 0; p < NUM_PHASES; p++) {
            if (p > 0) kputc('&', &tags);

            uint32_t flags = r->flags[t * NUM_PHASES + p];
            if (!flags) {
                kputc('.', &tags);
                continue;
            }

            kstring_t flag_names = {0};
            get_prediction_flags(flags, &flag_names);
            kputs(flag_names.s, &tags);
            free(flag_names.s);
        }
    }

    bcf_update_info_string(s->hdr, r->v, PREDICTION_FLAGS_TAG, tags.s);
    free(tags.s);
}

void write_tags(const Session *s, const Records *records) {
    for (int i = 0; i < records->n; i++) {
        Record *r = kv_A(*records, i);

        // A record without coding state (no coding transcripts, or skipped like one with more than one ALT) takes no part in any
        // phase and is written unchanged. Its triggers can still be set from the VEP consequences.
        if (r->tscripts == NULL) continue;

        write_score_tag(s, r, r->native_scores, SPLICEAI_TAG);
        write_score_tag(s, r, r->diff_scores, CS_TAG);
        write_score_tag(s, r, r->haplotype_scores, SPLICEAI_HAP_TAG);
        write_rescue_tags(s, r);
        write_ptv_tags(s, r);
        write_prediction_flags_tag(s, r);
    }
}
