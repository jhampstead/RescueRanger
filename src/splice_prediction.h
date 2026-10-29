#ifndef SPLICE_PREDICTION_HEADER
#define SPLICE_PREDICTION_HEADER

#include <stdbool.h>
#include <stdint.h>
#include <klib/kvec.h>
#include <cpliceai/predict.h>

#include "record.h"
#include "session.h"

// A predicted splice site this many bases or fewer from the canonical site of an exon is a compensatory site for it
#define PROXIMITY_THRESHOLD 15

enum SpliceType {
    ACCEPTOR = 1,
    DONOR = 2
};

// A splice site that a prediction found, with the record whose prediction it is. canonical_position is the position of the
// canonical site of the exon it belongs to, or 0 while that is not known.
typedef struct SpliceProbability {
    Record *record;
    enum SpliceType type;
    float ref_prior_probability;
    float prior_probability;
    float posterior_probability;
    uint32_t position;
    uint32_t canonical_position;
} SpliceProbability;

typedef kvec_t(SpliceProbability) SpliceProbabilities;

// Predicts the splice sites of an encoded sequence of a strand with the models. The caller frees the predictions.
void predict_for_seq(Model *models, const char *seq, int seq_len, int strand, float *predictions[], int *num_predictions);

// The reference predictions of the 2 * prediction_radius + 1 bases around the 0-based position, from the precomputed reference
// scores of the session, or NULL when the session has none or they do not have the gene (then the caller predicts the sequence).
// The same layout as the predictions of predict_on_sequence: NUM_SCORES floats per base, which the caller frees.
float *precomputed_reference_predictions(const Session *s, const Gene *gene, const char *chrom, const int64_t position);

// Predicts the window of a model width around the position in a sequence. rs, if given, are the records the sequence carries:
// they move the window. The caller frees the predictions.
float *predict_on_sequence(Model *models, const char *seq, const int width, uint32_t position, Range prediction_boundaries, int strand, Records *rs, const int prediction_radius);

// Of the NUM_SCORES values of every position, keeps the largest and sets the others to 0. num_probabilities is the length of
// the array.
void filter_to_most_likely_prediction(size_t num_probabilities, float *probabilities);

// Scales the values that are not 0 so that they sum to 1, and returns how many there are. posterior_probabilities gets the
// scaled values, and 0 where the probability is 0.
int transform_to_posterior_probabilities(const size_t num_probabilities, const float *probabilties, float *posterior_probabilities);

// Whether a predicted position is within PROXIMITY_THRESHOLD of a canonical site
int check_proximity(uint32_t predicted_position, uint32_t canonical_site_position);

// Adds a record to a vector of records unless it is in it
void kv_push_new_record(Records *rs, Record *r);

// Stores a predicted splice site on the sites of an exon. A lost canonical site is marked lost and predicted; a canonical site
// that is the only one in the window is marked predicted; a site next to the canonical one becomes a canonical site of its own,
// and any other site becomes a cryptic one; a site that is stored already only gets the record added. predict_for_pair is set when
// the pair of the site is to be predicted, and predict_for_competitor when its competitor is.
void store_splice_sites(SpliceProbability prediction, int *predict_for_pair, int *predict_for_competitor, SpliceSites *splice_sites);

// Whether a predicted splice site lies in the splice window of an exon (see Exon_splice_window). If it does, its canonical
// position is set to the position of the canonical site of that exon.
bool locate_prediction_on_exon(const Transcript *t, const size_t ei, SpliceProbability *prediction);

// Gets the probabilities of the acceptor sites and the donor sites of the transcripts that current_gene marks (one flag per
// transcript of the record) in the window pred_range, by comparing the predictions of the reference with those of the haplotype
// hap_in, and appends them to probabilities (a position once). Flags the record for the sites it finds: weak, ablated and
// cryptic ones. The predictions are NUM_SCORES values per position of the window. hap_offset is how many bases shorter the
// haplotype is than the reference in the window.
void collect_splice_probabilities(Record *r, const int phase, int *current_gene, const Range pred_range, const int prediction_radius, float *ref_predictions, float *hap_alt_predictions, const int hap_offset, Records *hap_in, SpliceProbabilities *probabilities);

// Gives an exon the canonical splice sites it has in the GFF. ei is the index of the exon in the transcript. Whether a
// deletion removes a site is not decided here, it follows from the predictions.
void init_exon_splice_sites(Exon *e, const size_t ei, const Transcript *t);

// Stage 1 of a phase: sets the canonical splice sites of the exons of every transcript of every record. Each transcript copy is
// initialised exactly once, because every record entry has its own copy for the phase.
void init_splice_sites_for_phase(const Records *records);

// Stage 2 of a phase: predicts the effect of every record that triggers a splice rescue check on the splice sites of its genes. This
// compares the reference, the record on its own and the record on its haplotype, stores the scores and flags, and stores the cryptic
// and competing splice sites that are found on the exons of the transcripts, together with the pairs and competitors of those.
void predict_splice_changes_for_phase(const Session *s, const int phase, const Records *records);

#endif
