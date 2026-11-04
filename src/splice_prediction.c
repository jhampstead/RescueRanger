#include "splice_prediction.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cpliceai/range.h>
#include <cpliceai/utils.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>

#include "bcftools/gff.h"
#include "haplotype.h"
#include "logging/log.h"
#include "phase.h"
#include "scores.h"
#include "triggers.h"
#include "utils.h"

void predict_for_seq(Model *models, const char *seq, int seq_len, int strand, float *predictions[], int *num_predictions) {
#ifdef CPLICEAI_HAS_PREDICT_PADDED_SEQUENCE
    // The library encodes the sequence, predicts and, on the negative strand, turns the predictions around. Its strand is
    // a character, and an unknown strand is read as the positive one, as it was here
    if (predict_padded_sequence(models, seq, seq_len, strand == STRAND_REV ? '-' : '+', predictions, num_predictions) != EXIT_SUCCESS) {
        log_info("Predict failed.");
        *predictions = NULL;
        *num_predictions = 0;
    }
#else
    // The older cpliceai has no predict_padded_sequence, and its one_hot_encode allocates the encoding
    float *ohe;
    int enc_len = one_hot_encode(seq, seq_len, &ohe);

    if (strand == STRAND_REV) reverse_encoding(ohe, enc_len);

    if (predict(models, enc_len, 1, &ohe, num_predictions, predictions)) {
        log_info("Predict failed.");
        free(ohe);
        *predictions = NULL;
        *num_predictions = 0;
        return;
    }
    free(ohe);

    if (strand == STRAND_REV) reverse_prediction(*predictions, *num_predictions, 3);
#endif
}

void filter_to_most_likely_prediction(size_t num_probabilities, float *probabilities) {
    for (int i = 0; i < num_probabilities; i += NUM_SCORES) {
        float max = 0;
        int max_idx = 0;

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if (probabilities[pred_idx] <= max) continue;

            max = probabilities[pred_idx];
            max_idx = pred_idx;
        }

        for (int pred_idx = i; pred_idx < i+NUM_SCORES; pred_idx++) {
            if (pred_idx == max_idx) continue;

            probabilities[pred_idx] = 0;
        }
    }
}

int transform_to_posterior_probabilities(const size_t num_probabilities, const float *probabilties, float *posterior_probabilities) {
    float sum = 0;
    int num_values = 0;

    for (int i = 0; i < num_probabilities; i++) {
        if (probabilties[i] == 0) continue;

        num_values++;
        sum += probabilties[i];
    }
    memset(posterior_probabilities, 0.0, num_probabilities * sizeof(float));

    for (int i = 0; i < num_probabilities; i++) {
        if (probabilties[i] == 0) continue;

        posterior_probabilities[i] = probabilties[i] / sum;
    }

    return num_values;
}

Range sequence_in_window(uint32_t position, Range prediction_boundaries, Records *rs, const int prediction_radius) {
    int64_t gene_start = prediction_boundaries.start;
    int64_t gene_end = prediction_boundaries.end;
    int64_t vposition = position;

    if (rs != NULL) {
        int start_pred_offset = haplotype_offset(rs, range_before(vposition, prediction_radius));
        int end_pred_offset = haplotype_offset(rs, (Range) { vposition, vposition + prediction_radius + 1 }); // The rest of the window, from the site on
        gene_start -= start_pred_offset;
        gene_end -= (haplotype_offset(rs, prediction_boundaries) - end_pred_offset);
        vposition -= haplotype_offset(rs, (Range) { prediction_boundaries.start, position });

    }
    int32_t start_offset = variant_position_in_sequence(prediction_radius) + (gene_start - vposition);
    int32_t end_offset = variant_position_in_sequence(prediction_radius) - (gene_end - (vposition + 1)); // gene_end is exclusive

    return (Range) { start_offset, end_offset };
}

float * predict_on_sequence(Model *models, const char *seq, const int width, uint32_t position, Range prediction_boundaries, int strand, Records *rs, const int prediction_radius) {
    const Range offsets = sequence_in_window(position, prediction_boundaries, rs, prediction_radius);

    const char *shifted_seq = seq - offsets.start;
    Range padding_boundary = { offsets.start > 0 ? offsets.start : 0, offsets.end > 0 ? offsets.end : 0 };
    char *padded_seq = pad_sequence(shifted_seq, padding_boundary, width);
    float *predictions; int num_predictions;
    predict_for_seq(models, padded_seq, width, strand, &predictions, &num_predictions);
    free(padded_seq);

    return predictions;
}

// The pairs of the reference and the haplotype prediction of one kind of site (ACCEPTOR_POS or DONOR_POS) over the window
static void get_prediction_pairs_for_site(int site, float *ref_predictions, size_t ref_len, float *hap_alt_predictions, size_t hap_alt_len, Records *rs, Range pred_range, size_t *map_size, PairPrediction **map) {
    float *ref_site_predictions = calloc(ref_len, sizeof(float));
    for (int i = 0; i < ref_len; i++) ref_site_predictions[i] = ref_predictions[NUM_SCORES * i + site];

    float *ref_site_prediction_posteriors = calloc(ref_len, sizeof(float));
    transform_to_posterior_probabilities(ref_len, ref_site_predictions, ref_site_prediction_posteriors);

    float *halt_site_predictions = calloc(hap_alt_len, sizeof(float));
    for (int i = 0; i < hap_alt_len; i++) halt_site_predictions[i] = hap_alt_predictions[NUM_SCORES * i + site];

    float *halt_site_prediction_posteriors = calloc(hap_alt_len, sizeof(float));
    transform_to_posterior_probabilities(hap_alt_len, halt_site_predictions, halt_site_prediction_posteriors);

    PairPrediction *ref_alt_map = NULL;
    size_t ref_alt_map_size = 0;
    parse_prediction_pairs(ref_site_predictions, ref_site_prediction_posteriors, halt_site_predictions, halt_site_prediction_posteriors, rs, pred_range, &ref_alt_map, &ref_alt_map_size);

    free(ref_site_predictions);
    free(ref_site_prediction_posteriors);
    free(halt_site_predictions);
    free(halt_site_prediction_posteriors);

    *map = ref_alt_map;
    *map_size = ref_alt_map_size;
}

static void kv_push_new_probability(SpliceProbabilities *ps, SpliceProbability p) {
    for (int i = 0; i < ps->n; i++) if (p.position == kv_A(*ps, i).position) return;

    kv_push(SpliceProbability, *ps, p);
}

static void get_splice_probabilities_from_map(Record *r, enum SpliceType type, int record_phase, int *gene_mask, Range pred_range, size_t map_size, PairPrediction *map, SpliceProbabilities *probabilities) {
    for (int t = 0; t < r->tscripts->n; t++) {
        if (gene_mask[t] == 0) continue;

        Gene *transcript_gene = kv_A(*r->transcript_genes, t);
        Transcript *transcript = kv_A(*r->tscripts, t);

        int *canonical_sites = calloc(map_size, sizeof(int));

        for (int e = 0; e < kv_size(transcript->exons); e++) {
            Exon *exon = kv_A(transcript->exons, e);

            int exon_boundary_position = type == DONOR ? Exon_donor_base(exon, transcript->strand) : Exon_acceptor_base(exon, transcript->strand);

            if ((type == DONOR) & ((transcript->strand == STRAND_REV && e == 0) | (transcript->strand == STRAND_FWD && e >= kv_size(transcript->exons)-1))) continue;
            if ((type == ACCEPTOR) & ((transcript->strand == STRAND_REV && e >= kv_size(transcript->exons)-1) | (transcript->strand == STRAND_FWD && e == 0))) continue;
            if (!((exon_boundary_position >= pred_range.start) & (exon_boundary_position < pred_range.end))) continue;

            for (int i = 0; i < map_size; i++) {
                PairPrediction pair = map[i];
                if (pair.ref.missing) continue;
                if (pair.ref.position != exon_boundary_position) continue;

                canonical_sites[i] = 1;

                if (pair.ref.prior < 0.01) {
                    type == DONOR ? set_record_flag(r, t, record_phase, WEAK_CANONICAL_DONOR) : set_record_flag(r, t, record_phase, WEAK_CANONICAL_ACCEPTOR);

                    if (pair.alt.missing | (pair.alt.prior < 0.01)) continue;

                    SpliceProbability p = { r, type, pair.ref.prior, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, exon_boundary_position };
                    kv_push_new_probability(probabilities, p);
                } else {
                    // If splice site has weakened
                    if (pair.ref.prior - pair.alt.prior > 0.1) {
                        type == DONOR ? set_record_flag(r, t, record_phase, ABLATED_CANONICAL_DONOR) : set_record_flag(r, t, record_phase, ABLATED_CANONICAL_ACCEPTOR);
                    }

                    SpliceProbability p = { r, type, pair.ref.prior, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, exon_boundary_position};
                    kv_push_new_probability(probabilities, p);
                }
            }
        }

        for (int i = 0; i < map_size; i++) {
            if (canonical_sites[i]) continue; // Canonical sites have already been checked, only handle cryptics here

            PairPrediction pair = map[i];

            if (pair.ref.prior == 0) {
                if (pair.alt.prior == 0) continue;

                type == DONOR ? set_record_flag(r, t, record_phase, CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, CRYPTIC_ACCEPTOR);
            } else {
                if (pair.alt.prior > 0) {
                    type == DONOR ? set_record_flag(r, t, record_phase, CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, CRYPTIC_ACCEPTOR);
                } else {
                    type == DONOR ? set_record_flag(r, t, record_phase, ABLATED_CRYPTIC_DONOR) : set_record_flag(r, t, record_phase, ABLATED_CRYPTIC_ACCEPTOR);
                }
            }

            SpliceProbability p = { r, type, pair.ref.prior, pair.alt.prior, pair.alt.posterior, pair.alt.missing ? pair.ref.position : pair.alt.position, 0 };
            kv_push_new_probability(probabilities, p);
        }

        free(canonical_sites);
    }
}

int check_proximity(uint32_t predicted_position, uint32_t canonical_site_position) {
    return abs((int32_t)predicted_position - (int32_t)canonical_site_position) <= PROXIMITY_THRESHOLD;
}

void kv_push_new_record(Records *rs, Record *r) {
    for (int i = 0; i < rs->n; i++) if (r == kv_A(*rs, i)) return;

    kv_push(Record *, *rs, r);
}

void store_splice_sites(SpliceProbability prediction, int *predict_for_pair, int *predict_for_competitor, SpliceSites *splice_sites) {
    if (prediction.type == DONOR && prediction.posterior_probability == 0 && prediction.position == splice_sites->canonical_donors[0].position) {
        kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
        splice_sites->canonical_donors[0].predicted = 1;
        splice_sites->canonical_donors[0].ref_prior_probability = prediction.ref_prior_probability;
        splice_sites->canonical_donors[0].prior_probability = prediction.prior_probability;
        splice_sites->canonical_donors[0].posterior_probability = prediction.posterior_probability;
        splice_sites->canonical_donors[0].lost = 1;
        *predict_for_pair = 1;
        *predict_for_competitor = 1;
    } else if (prediction.type == ACCEPTOR && prediction.posterior_probability == 0 && prediction.position == splice_sites->canonical_acceptors[0].position) {
        kv_push_new_record((Records *) splice_sites->canonical_acceptors[0].records_ptr, prediction.record);
        splice_sites->canonical_acceptors[0].predicted = 1;
        splice_sites->canonical_acceptors[0].ref_prior_probability = prediction.ref_prior_probability;
        splice_sites->canonical_acceptors[0].prior_probability = prediction.prior_probability;
        splice_sites->canonical_acceptors[0].posterior_probability = prediction.posterior_probability;
        splice_sites->canonical_acceptors[0].lost = 1;
        *predict_for_pair = 1;
        *predict_for_competitor = 1;
    } else if (prediction.type == DONOR & prediction.posterior_probability > 0) {
        if (splice_sites->canonical_donors[0].gff_canonical && prediction.posterior_probability == 1.0 && prediction.position == splice_sites->canonical_donors[0].position) { // If canonical exists and is the only splice site within the window
            kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
            splice_sites->canonical_donors[0].predicted = 1;
            splice_sites->canonical_donors[0].ref_prior_probability = prediction.ref_prior_probability;
            splice_sites->canonical_donors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_donors[0].posterior_probability = prediction.posterior_probability;
            *predict_for_pair = 1;
            *predict_for_competitor = 1;
        } else if (splice_sites->canonical_donors[0].gff_canonical && check_proximity(prediction.position, splice_sites->canonical_donors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical donor site
            int canonical_donor_exists = 0;
            int canonical_donor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_donors; i++) {
                if (splice_sites->canonical_donors[i].position == prediction.position) {
                    canonical_donor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->canonical_donors[i].records_ptr, prediction.record);
                    if (splice_sites->canonical_donors[i].predicted == 1) {
                        canonical_donor_predicted = 1;
                    }
                }
            }

            if (canonical_donor_exists & canonical_donor_predicted) return;
            if (canonical_donor_exists & !canonical_donor_predicted) {
                *predict_for_pair = 1;
                for (int i = 0; i < splice_sites->n_canonical_donors; i++) {
                    if (splice_sites->canonical_donors[i].position == prediction.position) {
                        kv_push_new_record((Records *) splice_sites->canonical_donors[0].records_ptr, prediction.record);
                        splice_sites->canonical_donors[i].predicted = 1;
                        splice_sites->canonical_donors[i].ref_prior_probability = prediction.ref_prior_probability;
                        splice_sites->canonical_donors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_donors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_donors >= splice_sites->m_canonical_donors) {
                    resize_splice_site(&splice_sites->canonical_donors, &splice_sites->m_canonical_donors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_donors[splice_sites->n_canonical_donors++] = site;
            }
        // Store cryptic donor site
        } else {
            // Check if prediction position already exists as a cryptic donor site
            int cryptic_donor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_donors; i++) {
                if (splice_sites->cryptic_donors[i].position == prediction.position) {
                    cryptic_donor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->cryptic_donors[i].records_ptr, prediction.record);
                }
            }

            if (!cryptic_donor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_donors >= splice_sites->m_cryptic_donors) {
                    resize_splice_site(&splice_sites->cryptic_donors, &splice_sites->m_cryptic_donors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);

                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_donors[splice_sites->n_cryptic_donors++] = site;
            }
        }
    } else if (prediction.type == ACCEPTOR & prediction.posterior_probability > 0) {
        // Store canonical acceptor site
        if (splice_sites->canonical_acceptors[0].gff_canonical && prediction.posterior_probability == 1.0 && prediction.position == splice_sites->canonical_acceptors[0].position) { // If canonical exists and is the only splice site within the window
            kv_push_new_record((Records *) splice_sites->canonical_acceptors[0].records_ptr, prediction.record);
            splice_sites->canonical_acceptors[0].predicted = 1;
            splice_sites->canonical_acceptors[0].ref_prior_probability = prediction.ref_prior_probability;
            splice_sites->canonical_acceptors[0].prior_probability = prediction.prior_probability;
            splice_sites->canonical_acceptors[0].posterior_probability = prediction.posterior_probability;
            *predict_for_pair = 1;
            *predict_for_competitor = 1;
        } else if (splice_sites->canonical_acceptors[0].gff_canonical & check_proximity(prediction.position, splice_sites->canonical_acceptors[0].position)) { // If canonical exists and site is compensatory for canonical
            // Check if prediction position already exists as a canonical acceptor site
            int canonical_acceptor_exists = 0;
            int canonical_acceptor_predicted = 0;
            for (int i = 0; i < splice_sites->n_canonical_acceptors; i++) {
                if (splice_sites->canonical_acceptors[i].position == prediction.position) {
                    canonical_acceptor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->canonical_acceptors[i].records_ptr, prediction.record);
                    if (splice_sites->canonical_acceptors[i].predicted == 1) {
                        canonical_acceptor_predicted = 1;
                    }
                }
            }

            if (canonical_acceptor_exists & canonical_acceptor_predicted) return;
            if (canonical_acceptor_exists & !canonical_acceptor_predicted) {
                *predict_for_pair = 1;
                for (int i = 0; i < splice_sites->n_canonical_acceptors; i++) {
                    if (splice_sites->canonical_acceptors[i].position == prediction.position) {
                        kv_push_new_record((Records *) splice_sites->canonical_acceptors[i].records_ptr, prediction.record);
                        splice_sites->canonical_acceptors[i].predicted = 1;
                        splice_sites->canonical_acceptors[i].ref_prior_probability = prediction.ref_prior_probability;
                        splice_sites->canonical_acceptors[i].prior_probability = prediction.prior_probability;
                        splice_sites->canonical_acceptors[i].posterior_probability = prediction.posterior_probability;
                    }
                }
            } else {
                *predict_for_pair = 1;
                if (splice_sites->n_canonical_acceptors >= splice_sites->m_canonical_acceptors) {
                    resize_splice_site(&splice_sites->canonical_acceptors, &splice_sites->m_canonical_acceptors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->canonical_acceptors[splice_sites->n_canonical_acceptors++] = site;
            }
        // Store cryptic acceptor site
        } else {
            // Check if prediction position already exists as a cryptic acceptor site
            int cryptic_acceptor_exists = 0;
            for (int i = 0; i < splice_sites->n_cryptic_acceptors; i++) {
                if (splice_sites->cryptic_acceptors[i].position == prediction.position) {
                    cryptic_acceptor_exists = 1;

                    kv_push_new_record((Records *) splice_sites->cryptic_acceptors[i].records_ptr, prediction.record);
                }
            }

            if (!cryptic_acceptor_exists) {
                *predict_for_competitor = 1;
                if (splice_sites->n_cryptic_acceptors >= splice_sites->m_cryptic_acceptors) {
                    resize_splice_site(&splice_sites->cryptic_acceptors, &splice_sites->m_cryptic_acceptors);
                }
                Records *splice_records = malloc(sizeof(Records));
                kv_init(*splice_records);
                kv_push(Record *, *splice_records, prediction.record);
                SpliceSite site = { splice_records, prediction.position, prediction.canonical_position, 1, 0, 0, prediction.ref_prior_probability, prediction.prior_probability, prediction.posterior_probability };
                splice_sites->cryptic_acceptors[splice_sites->n_cryptic_acceptors++] = site;
            }
        }
    }
}

// Writes the positions of the window where one of the predictions has a probability above the rescue threshold to the raw
// scores file
static void write_predictions(FILE *scores_tsv, const char *chr, const int64_t pos, const char *ref, const char *alt, const char *gene_name, const char *transcript_name, const int phase, const char *type, const char *pred_tags[], const int n_predictions, const float *predictions[], const int prediction_radius) {
    if (scores_tsv == NULL) return; // The raw scores are optional

    for (int i = 0; i < pred_width(prediction_radius); i++) {
        const int offset = i - prediction_radius;

        bool should_print = false;
        for (int j = 0; j < n_predictions; j++) {
            const float acceptor = predictions[j][i * NUM_SCORES + ACCEPTOR];
            const float donor = predictions[j][i * NUM_SCORES + DONOR];

            if ((acceptor > SPLICEAI_RESCUE_THRESHOLD) | (donor > SPLICEAI_RESCUE_THRESHOLD)) {
                should_print = true;
            }

        }
        for (int j = 0; j < n_predictions; j++) {
            if (!should_print) continue;

            const float neither = predictions[j][i * NUM_SCORES];
            const float acceptor = predictions[j][i * NUM_SCORES + ACCEPTOR];
            const float donor = predictions[j][i * NUM_SCORES + DONOR];

            fprintf(scores_tsv, "%s\t%li\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
                    chr, pos + 1, ref, alt, gene_name, transcript_name, phase, offset, pred_tags[j], type, donor, acceptor, neither);
        }
    }
}

void init_exon_splice_sites(Exon *e, const size_t ei, const Transcript *t) {
    SpliceSites *sites = malloc(sizeof(SpliceSites));
    init_splice_site(&sites->canonical_acceptors, &sites->n_canonical_acceptors, &sites->m_canonical_acceptors);
    init_splice_site(&sites->canonical_donors, &sites->n_canonical_donors, &sites->m_canonical_donors);
    init_splice_site(&sites->cryptic_acceptors, &sites->n_cryptic_acceptors, &sites->m_cryptic_acceptors);
    init_splice_site(&sites->cryptic_donors, &sites->n_cryptic_donors, &sites->m_cryptic_donors);

    if (t->strand == STRAND_FWD || t->strand == STRAND_REV) {
        const bool reverse = t->strand == STRAND_REV;
        const bool has_acceptor = reverse ? ei < kv_size(t->exons)-1 : ei > 0;
        const bool has_donor = reverse ? ei > 0 : ei < kv_size(t->exons)-1;

        if (has_acceptor) {
            Records *splice_records = malloc(sizeof(Records));
            kv_init(*splice_records);
            const int64_t acceptor = Exon_acceptor_base(e, t->strand);
            sites->canonical_acceptors[sites->n_canonical_acceptors++] = (SpliceSite) { splice_records, acceptor, acceptor, 0, 1, 0, 1.0, 1.0, 1.0 };
        }
        if (has_donor) {
            Records *splice_records = malloc(sizeof(Records));
            kv_init(*splice_records);
            const int64_t donor = Exon_donor_base(e, t->strand);
            sites->canonical_donors[sites->n_canonical_donors++] = (SpliceSite) { splice_records, donor, donor, 0, 1, 0, 1.0, 1.0, 1.0 };
        }
    }

    e->splice_sites = sites;
}

void init_splice_sites_for_phase(const Records *records) {
    for (int i = 0; i < records->n; i++) {
        Record *r = kv_A(*records, i);

        for (int ti = 0; ti < r->tscripts->n; ti++) {
            Transcript *t = kv_A(*r->tscripts, ti);
            Gene *gene = kv_A(*r->transcript_genes, ti);
            if (t == NULL || gene == NULL) continue;
            if (!record_overlaps_gene(r, gene)) continue;

            for (int ei = 0; ei < kv_size(t->exons); ei++) {
                init_exon_splice_sites(kv_A(t->exons, ei), ei, t);
            }
        }
    }
}

// Predicts the splice sites of a sequence that is padded to the width of the model. The caller frees the predictions. (Not
// predict_padded_sequence: cpliceai has a function of that name, which is what predict_for_seq calls.)
static float *predict_padded_window(Model *models, char *seq, const Range padding, const int width, const int strand, int *num_predictions) {
    char *padded = pad_sequence(seq, padding, width);

    float *predictions = NULL;
    predict_for_seq(models, padded, strlen(padded), strand, &predictions, num_predictions);
    free(padded);

    return predictions;
}

float *precomputed_reference_predictions(const Session *s, const Gene *gene, const char *chrom, const int64_t position) {
    return ReferenceScores_window(s->reference, chrom, gene->name, gene->coords.start, gene->coords.end, gene->strand == STRAND_REV ? '-' : '+', position, s->prediction_radius);
}

// How much the haplotype changes the score of a splice site compared to the variant on its own
static float score_difference(const float *hap_alt_predictions, const float *alt_predictions, const int index, const int type, const int prediction_radius) {
    return hap_alt_predictions[(index + prediction_radius) * NUM_SCORES + type] - alt_predictions[(index + prediction_radius) * NUM_SCORES + type];
}

// Sets a flag for every transcript of the record that belongs to the gene
static void flag_transcripts_of_gene(Record *r, const Gene *gene, const int phase, const uint32_t flag) {
    for (int ti = 0; ti < r->transcript_genes->n; ti++) {
        Gene *transcript_gene = kv_A(*r->transcript_genes, ti);

        if (transcript_gene == NULL) continue;
        if (transcript_gene->gff_id != gene->gff_id) continue;

        set_record_flag(r, ti, phase, flag);
    }
}

bool locate_prediction_on_exon(const Transcript *t, const size_t ei, SpliceProbability *prediction) {
    Range positions;
    if (!Exon_splice_window(t, ei, prediction->type == DONOR, &positions)) return false;

    Range splice_position = (Range) { prediction->position, prediction->position+1 };
    Range intersection;
    if (!intersect_ranges(splice_position, positions, &intersection) || !range_len(intersection)) return false;

    Exon *e = kv_A(t->exons, ei);
    if (prediction->type == DONOR) {
        prediction->canonical_position = e->splice_sites->canonical_donors[0].position;
    } else { // ACCEPTOR
        prediction->canonical_position = e->splice_sites->canonical_acceptors[0].position;
    }

    return true;
}

// Stores every predicted splice site on every exon whose splice window it lies in. This can ask for the pair or the
// competitor of a site to be predicted too, which is reported through the two flags.
static void store_probabilities_on_exons(const Transcript *t, const SpliceProbabilities *probabilities, int *predict_for_pair, int *predict_for_competitor) {
    for (int pi = 0; pi < probabilities->n; pi++) {
        SpliceProbability prediction = probabilities->a[pi];

        for (int ei = 0; ei < kv_size(t->exons); ei++) {
            if (!locate_prediction_on_exon(t, ei, &prediction)) continue;

            store_splice_sites(prediction, predict_for_pair, predict_for_competitor, kv_A(t->exons, ei)->splice_sites);
        }
    }
}

void collect_splice_probabilities(Record *r, const int phase, int *current_gene, const Range pred_range, const int prediction_radius, float *ref_predictions, float *hap_alt_predictions, const int hap_offset, Records *hap_in, SpliceProbabilities *probabilities) {
    // Get acceptor pairs
    PairPrediction *ref_alt_acceptor_map = NULL;
    size_t ref_alt_acceptor_map_size = 0;
    get_prediction_pairs_for_site(ACCEPTOR_POS, ref_predictions, pred_width(prediction_radius), hap_alt_predictions, pred_width(prediction_radius)-hap_offset, hap_in, pred_range, &ref_alt_acceptor_map_size, &ref_alt_acceptor_map);
    get_splice_probabilities_from_map(r, ACCEPTOR, phase, current_gene, pred_range, ref_alt_acceptor_map_size, ref_alt_acceptor_map, probabilities);
    free(ref_alt_acceptor_map);


    // Get donor pairs
    PairPrediction *ref_alt_donor_map = NULL;
    size_t ref_alt_donor_map_size = 0;
    get_prediction_pairs_for_site(DONOR_POS, ref_predictions, pred_width(prediction_radius), hap_alt_predictions, pred_width(prediction_radius)-hap_offset, hap_in, pred_range, &ref_alt_donor_map_size, &ref_alt_donor_map);
    get_splice_probabilities_from_map(r, DONOR, phase, current_gene, pred_range, ref_alt_donor_map_size, ref_alt_donor_map, probabilities);
    free(ref_alt_donor_map);
}

// The kinds of splice sites that are predicted next to a predicted one
enum NeighbourKind {
    PAIR_SITE,       // The site it is spliced to: the acceptor after a donor, and the donor after an acceptor
    COMPETITOR_SITE  // The canonical site of the same kind on the next exon, which competes with it
};

// What predicting a neighbouring site needs from the record and the gene it is predicted for
typedef struct NeighbourContext {
    const Session *session;
    int phase;
    Record *r;
    Gene *gene;
    const char *chr;
    Records *just_alt;  // Only the record
    Records *hap_in;    // The record and the variants on its haplotype
    kstring_t *alt;     // The gene sequence with only the record applied
    kstring_t *hap_alt; // The gene sequence with the haplotype applied
} NeighbourContext;

// Predicts the pair or the competitor of a predicted splice site, which is the canonical site on the next exon, and
// stores the splice sites that this finds. Does nothing for a site that has been predicted before.
// Returns false if there is no neighbouring site to predict, or if it has been predicted already, which is where the
// caller moves on to the next exon.
static bool predict_neighbouring_site(const NeighbourContext *ctx, Transcript *t, const Gene *tgene, int *current_gene, const SpliceProbability *prediction, const size_t ei, const int direction, const enum NeighbourKind kind, int *predict_for_pair, int *predict_for_competitor) {
    const Session *s = ctx->session;
    const bool donor = prediction->type == DONOR;

    // A pair is on the next exon in the direction of the transcript after a donor, and on the previous one before an
    // acceptor. The competitor of a site is on the exon on the other side of that.
    const int step = (kind == PAIR_SITE) == donor ? direction : -direction;
    const int neighbour = (int) ei + step;
    if (neighbour < 0 || neighbour >= kv_size(t->exons)) return false;

    // A pair is of the other kind than the predicted site, and a competitor of the same kind
    const bool site_is_donor = (kind == COMPETITOR_SITE) == donor;
    SpliceSites *neighbour_sites = kv_A(t->exons, neighbour)->splice_sites;
    SpliceSite site = site_is_donor ? neighbour_sites->canonical_donors[0] : neighbour_sites->canonical_acceptors[0];

    kv_push_new_record((Records *) site.records_ptr, prediction->record);

    if (site.predicted) return false;

    Range pred_range = range_around(site.position, s->prediction_radius);
    float *ref_predictions = precomputed_reference_predictions(s, tgene, ctx->chr, site.position);
    if (ref_predictions == NULL) ref_predictions = predict_on_sequence(s->models, ctx->gene->seq.s, model_width(s->prediction_radius), site.position, tgene->coords, tgene->strand, NULL, s->prediction_radius);

    int pred_offset = haplotype_offset(ctx->just_alt, pred_range);
    float *alt_predictions = predict_on_sequence(s->models, ctx->alt->s, model_width(s->prediction_radius)-pred_offset, site.position, tgene->coords, tgene->strand, ctx->just_alt, s->prediction_radius);
    float *corrected_alt_predictions = match_predictions_to_reference(alt_predictions, pred_range, ctx->just_alt, pred_array_len(s->prediction_radius));

    int hap_pred_offset = haplotype_offset(ctx->hap_in, pred_range);
    float *hap_alt_predictions = predict_on_sequence(s->models, ctx->hap_alt->s, model_width(s->prediction_radius)-hap_pred_offset, site.position, tgene->coords, tgene->strand, ctx->hap_in, s->prediction_radius);
    float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, pred_range, ctx->hap_in, pred_array_len(s->prediction_radius));

    const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
    const float *predictions[3] = { ref_predictions, corrected_alt_predictions, corrected_hap_alt_predictions };
    write_predictions(s->scores_tsv, ctx->chr, site.position, ctx->r->v->d.allele[0], ctx->r->v->d.allele[1], ctx->gene->name, gff_id2string(s->gff, 0, t->id), ctx->phase+1, kind == PAIR_SITE ? "pair" : "competitor", prediction_names, 3, predictions, s->prediction_radius);

    filter_to_most_likely_prediction(pred_array_len(s->prediction_radius), ref_predictions);
    filter_to_most_likely_prediction((pred_width(s->prediction_radius) - hap_pred_offset) * NUM_SCORES, hap_alt_predictions);

    SpliceProbabilities probabilities;
    kv_init(probabilities);
    collect_splice_probabilities(ctx->r, ctx->phase, current_gene, pred_range, s->prediction_radius, ref_predictions, hap_alt_predictions, hap_pred_offset, ctx->hap_in, &probabilities);

    store_probabilities_on_exons(t, &probabilities, predict_for_pair, predict_for_competitor);

    kv_destroy(probabilities);
    free(ref_predictions);
    free(alt_predictions);
    free(corrected_alt_predictions);
    free(hap_alt_predictions);
    free(corrected_hap_alt_predictions);

    return true;
}

// For a record whose haplotype changes a splice site by enough: finds the splice sites of the window (the reference
// predictions and the haplotype predictions are filtered to the most likely prediction of every position first), places them on
// the exons of the transcripts of the gene, and predicts the pair and the competitor of every site that asks for it. The
// predictions are not freed here, the caller owns them.
static void predict_pairs_and_competitors(const Session *s, const int phase, Record *r, Gene *gene, const char *chr, const int64_t centre, float *ref_predictions, float *hap_alt_predictions, const int num_ref_predictions, const int num_hap_alt_predictions, const int pred_offset, Records *hap_in, Records *just_alt, kstring_t *alt, kstring_t *hap_alt) {
    // Remove predictions that are most likely to be neutral
    filter_to_most_likely_prediction(num_ref_predictions, ref_predictions);
    filter_to_most_likely_prediction(num_hap_alt_predictions, hap_alt_predictions);


    int in_gene[r->tscripts->n]; // Per transcript, whether it belongs to the gene
    for (int ti = 0; ti < r->tscripts->n; ti++) {
        Transcript *t = kv_A(*r->tscripts, ti);
        Gene *tgene = kv_A(*r->transcript_genes, ti);

        in_gene[ti] = (t != NULL && tgene->gff_id == gene->gff_id);
    }

    Range pred_range = range_around(centre, s->prediction_radius);
    SpliceProbabilities splice_probabilities;
    kv_init(splice_probabilities);
    collect_splice_probabilities(r, phase, in_gene, pred_range, s->prediction_radius, ref_predictions, hap_alt_predictions, pred_offset, hap_in, &splice_probabilities);

    NeighbourContext ctx = { s, phase, r, gene, chr, just_alt, hap_in, alt, hap_alt };

    for (int ti = 0; ti < r->tscripts->n; ti++) {
        Transcript *t = kv_A(*r->tscripts, ti);
        Gene *tgene = kv_A(*r->transcript_genes, ti);

        if (t == NULL || tgene->gff_id != gene->gff_id) continue;

        int only_this_transcript[r->tscripts->n]; // Per transcript, whether it is the current one
        for (int k = 0; k < r->tscripts->n; k++) only_this_transcript[k] = 0;
        only_this_transcript[ti] = 1;

        for (int si = 0; si < splice_probabilities.n; si++) {
            SpliceProbability prediction = splice_probabilities.a[si];
            for (int ei = 0; ei < kv_size(t->exons); ei++) {
                int direction = 0;
                if (t->strand == STRAND_FWD) direction = 1;
                if (t->strand == STRAND_REV) direction = -1;
                if (t->strand == STRAND_UNK) continue;

                if (!locate_prediction_on_exon(t, ei, &prediction)) continue;

                int predict_for_pair = 0;
                int predict_for_competitor = 0;
                store_splice_sites(prediction, &predict_for_pair, &predict_for_competitor, kv_A(t->exons, ei)->splice_sites);

                // Without a pair to predict, the competitor is skipped too
                if (predict_for_pair && !predict_neighbouring_site(&ctx, t, tgene, only_this_transcript, &prediction, ei, direction, PAIR_SITE, &predict_for_pair, &predict_for_competitor)) continue;

                // Storing the sites of a pair can ask for a competitor to be predicted too
                if (predict_for_competitor) predict_neighbouring_site(&ctx, t, tgene, only_this_transcript, &prediction, ei, direction, COMPETITOR_SITE, &predict_for_pair, &predict_for_competitor);
            }
        }
    }

    kv_destroy(splice_probabilities);
}

void predict_splice_changes_for_phase(const Session *s, const int phase, const Records *records) {
    for (int i = 0; i < records->n; i++) {
        Record *r = kv_A(*records, i);

        int splice_triggered = 0;
        for (int j = 0; j < r->vep->n; j++) {
            if (kv_A(*r->vep, j)->triggers[SPLICE_TRIGGER]) splice_triggered = 1;
        }
        if (!splice_triggered) continue;

        log_info("%s:%li (splice)", bcf_hdr_id2name(s->hdr, r->v->rid), r->v->pos + 1);

        if (!variant_fits_window(&r->alleles, s->prediction_radius)) {
            r->skipped = 1;
            log_warn("Skipping on %s:%li, variant is too long (%d alt bases, %d ref bases, --prediction-radius is %d)", bcf_hdr_id2name(s->hdr, r->v->rid), r->v->pos + 1, r->alleles.alt_len, r->alleles.ref_len, s->prediction_radius);
            continue;
        }

        for (int g = 0; g < r->genes->n; g++) {
            Gene *gene = kv_A(*r->genes, g);
            if (gene == NULL) continue;

            for (int j = 0; j < records->n; j++) {
                Record *other = kv_A(*records, j);
                if (has_matching_phase(other->v, s->hdr, s->gt_tag, phase, false)) {
                    continue;
                }

                other->unphased = 1;
            }

            if (!record_overlaps_gene(r, gene)) continue;

            // Find records for the current allele
            Records hap_in; // Variants for current allele including current variant
            Records hap_ex; // Variant for current allele excluding current variant
            find_haplotype_records(records, r, gene, true, &hap_in, &hap_ex);

            Records just_alt;
            kv_init(just_alt);
            kv_push(Record *, just_alt, r);

            kstring_t alt;
            get_haplotype_sequence(gene->seq, &just_alt, gene->coords, &alt);

            kstring_t hap_alt;
            get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_alt);


            // The prediction windows are centred on the variant's first base: for an indel the base it is attached to (the VCF
            // position of a normalised indel)
            const int64_t centre = variant_footprint(&r->alleles).start;

            // REF
            int32_t start = variant_position_in_sequence(s->prediction_radius) - (centre - gene->coords.start);
            int32_t end = variant_position_in_sequence(s->prediction_radius) - (gene->coords.end - (centre+1)); // Open ended, so + 1

            char *ref_seq = gene->seq.s - start; // Pointer to existing memory, no free necessary
            if (start < 0) start = 0;
            if (end < 0) end = 0;

            // The haplotype with the unphased records, which hap_in leaves out
            Records hap_with_unphased;
            find_haplotype_records(records, r, gene, false, &hap_with_unphased, NULL);
            flag_phasing_in_haplotype(r, &hap_with_unphased, range_around(centre, s->prediction_radius), s->prediction_radius, phase);
            kv_destroy(hap_with_unphased);

            // The reference predictions come from the precomputed reference scores when there are some and they have this gene
            const char *chr = bcf_seqname(s->hdr, r->v);
            int num_ref_predictions = pred_array_len(s->prediction_radius);
            float *ref_predictions = precomputed_reference_predictions(s, gene, chr, centre);
            if (ref_predictions == NULL) {
                ref_predictions = predict_padded_window(s->models, ref_seq, (Range) { start, end }, model_width(s->prediction_radius), gene->strand, &num_ref_predictions);
            }

            // ALT, padded with the start and end of the REF window: the width and the gene end shrink by the same net offset, so
            // the padding on the right is the same as for the REF
            int alt_pred_offset = allele_offset(&r->alleles, range_around(centre, s->prediction_radius));

            int32_t alt_start = variant_position_in_sequence(s->prediction_radius) - (centre - gene->coords.start);

            char *alt_seq = alt.s - alt_start;

            int num_alt_predictions = 0;
            float *alt_predictions = predict_padded_window(s->models, alt_seq, (Range) { start, end }, model_width(s->prediction_radius) - alt_pred_offset, gene->strand, &num_alt_predictions);

            float *corrected_alt_predictions = match_predictions_to_reference(alt_predictions, range_around(centre, s->prediction_radius), &just_alt, num_ref_predictions);

            // HAP ALT
            int pred_offset = prediction_offset(&r->alleles, &hap_in, s->prediction_radius);

            int64_t gene_start = gene->coords.start;
            int64_t gene_end = gene->coords.end - haplotype_offset(&hap_in, gene->coords);
            int64_t vposition = centre - haplotype_offset(&hap_ex, (Range) {gene->coords.start, centre});

            int64_t start_offset = haplotype_offset(&hap_in, range_before(centre, s->prediction_radius));
            int32_t hap_alt_start = variant_position_in_sequence(s->prediction_radius) + (gene_start - vposition) - start_offset;
            int32_t hap_alt_end = variant_position_in_sequence(s->prediction_radius) - (gene_end - (vposition + 1)) - (pred_offset - start_offset); // gene_end is exclusive

            char *hap_alt_seq = hap_alt.s - hap_alt_start;
            if (hap_alt_start < 0) hap_alt_start = 0;
            if (hap_alt_end < 0) hap_alt_end = 0;

            int num_hap_alt_predictions = 0;
            float *hap_alt_predictions = predict_padded_window(s->models, hap_alt_seq, (Range) {hap_alt_start, hap_alt_end}, model_width(s->prediction_radius) - pred_offset, gene->strand, &num_hap_alt_predictions);
            float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, range_around(centre, s->prediction_radius), &hap_in, num_ref_predictions);


            // Calculate scores between ref and alt
            Score score = calculate_delta_scores(r->v->d.allele[1], gene->name, ref_predictions, corrected_alt_predictions, num_ref_predictions, s->prediction_radius);
            r->native_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
            *r->native_scores[g * NUM_PHASES + phase] = score;


            // Calculate scores between ref and hap_alt
            Score hap_score = calculate_delta_scores(r->v->d.allele[1], gene->name, ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, s->prediction_radius);
            r->haplotype_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
            *r->haplotype_scores[g * NUM_PHASES + phase] = hap_score;

            // Score diff;
            r->diff_scores[g * NUM_PHASES + phase] = malloc(sizeof(Score));
            Score diff = {0};
            diff.alt = r->v->d.allele[1];
            diff.gene = gene->name;
            diff.ag_idx = hap_score.ag_idx;
            diff.ag = score_difference(corrected_hap_alt_predictions, corrected_alt_predictions, diff.ag_idx, ACCEPTOR, s->prediction_radius);

            diff.al_idx = hap_score.al_idx;
            diff.al = score_difference(corrected_hap_alt_predictions, corrected_alt_predictions, diff.al_idx, ACCEPTOR, s->prediction_radius);

            diff.dg_idx = hap_score.dg_idx;
            diff.dg = score_difference(corrected_hap_alt_predictions, corrected_alt_predictions, diff.dg_idx, DONOR, s->prediction_radius);

            diff.dl_idx = hap_score.dl_idx;
            diff.dl = score_difference(corrected_hap_alt_predictions, corrected_alt_predictions, diff.dl_idx, DONOR, s->prediction_radius);


            *r->diff_scores[g * NUM_PHASES + phase] = diff;

            if (
                diff.ag > SPLICEAI_RESCUE_THRESHOLD |
                diff.ag < -SPLICEAI_RESCUE_THRESHOLD |
                diff.al > SPLICEAI_RESCUE_THRESHOLD |
                diff.al < -SPLICEAI_RESCUE_THRESHOLD
            ) {
                flag_transcripts_of_gene(r, gene, phase, HAPLOTYPE_SPECIFIC_ACCEPTOR);
            }

            if (
                diff.dg > SPLICEAI_RESCUE_THRESHOLD |
                diff.dg < -SPLICEAI_RESCUE_THRESHOLD |
                diff.dl > SPLICEAI_RESCUE_THRESHOLD |
                diff.dl < -SPLICEAI_RESCUE_THRESHOLD
            ) {
                flag_transcripts_of_gene(r, gene, phase, HAPLOTYPE_SPECIFIC_DONOR);
            }

            const char *prediction_names[3] = { "ref", "alt", "hap_alt" };
            const float *predictions[3] = { ref_predictions, corrected_alt_predictions, corrected_hap_alt_predictions };
            write_predictions(s->scores_tsv, chr, r->v->pos, r->v->d.allele[0], r->v->d.allele[1], gene->name, ".", phase+1, "canonical", prediction_names, 3, predictions, s->prediction_radius);

            // If there's less than 0.1 change between ref and hap alt, don't predict for competitors and pairs
            const bool small_change = hap_score.ag < SPLICEAI_RESCUE_THRESHOLD && hap_score.al < SPLICEAI_RESCUE_THRESHOLD && hap_score.dg < SPLICEAI_RESCUE_THRESHOLD && hap_score.dl < SPLICEAI_RESCUE_THRESHOLD;
            if (!small_change) {
                predict_pairs_and_competitors(s, phase, r, gene, chr, centre, ref_predictions, hap_alt_predictions, num_ref_predictions, num_hap_alt_predictions, pred_offset, &hap_in, &just_alt, &alt, &hap_alt);
            }

            free(ref_predictions);
            free(alt_predictions);
            free(corrected_alt_predictions);
            free(hap_alt_predictions);
            free(corrected_hap_alt_predictions);

            free(alt.s);
            free(hap_alt.s);

            kv_destroy(just_alt);
            kv_destroy(hap_in);
            kv_destroy(hap_ex);
        }
    }
}
