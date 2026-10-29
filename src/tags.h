#ifndef TAGS_HEADER
#define TAGS_HEADER

#include <stdint.h>
#include <htslib/kstring.h>
#include <htslib/vcf.h>
#include <cpliceai/utils.h>

#include "record.h"
#include "session.h"

// The INFO tags that are written to the records, with their header lines
#undef SPLICEAI_HAP_TAG
#undef SPLICEAI_HAP_DESC
#define SPLICEAI_HAP_TAG "SpliceAI_Haplotype"
#define SPLICEAI_HAP_DESC "##INFO=<ID=SpliceAI_Haplotype,Number=.,Type=String,Description=\"CpliceAIv1.0 haplotype variant annotation from reference to complete haplotype. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"

#define CS_TAG "CS"
#define CS_DESC "##INFO=<ID=CS,Number=.,Type=String,Description=\"Change in splicing between haplotype + variant and haplotype at positions DP_AG, DP_AL, DP_DG, DP_DL. Format: CS_AG|CS_AL|CS_DG|CS_DL|CS_SUM\">"

#define RESCUE_TAG "RESCUE"
#define RESCUE_DESC "##INFO=<ID=RESCUE,Number=.,Type=String,Description=\"Boolean rescue prediction\">"

#define RESCUE_PROB_TAG "RESCUE_PROB"
#define RESCUE_PROB_DESC "##INFO=<ID=RESCUE_PROB,Number=.,Type=String,Description=\"Probability of rescue\">"

#define RESCUE_TYPE_TAG "RESCUE_TYPE"
#define RESCUE_TYPE_DESC "##INFO=<ID=RESCUE_TYPE,Number=.,Type=String,Description=\"Information on type of rescue. Format: RESCUE_TYPE|EXON_NUMBER|PROBABILITY|AFFECTED_BASES\">"

#define PTV_TAG "PTV"
#define PTV_DESC "##INFO=<ID=PTV,Number=.,Type=String,Description=\"Newly predicted PTV consequence\">"

#define PTV_PROB_TAG "PTV_PROB"
#define PTV_PROB_DESC "##INFO=<ID=PTV_PROB,Number=.,Type=String,Description=\"Probability of induced PTV\">"

#define PTV_TYPE_TAG "PTV_TYPE"
#define PTV_TYPE_DESC "##INFO=<ID=PTV_TYPE,Number=.,Type=String,Description=\"Type of PTV\">"

#define PREDICTION_FLAGS_TAG "PREDICTION_FLAGS"
#define PREDICTION_FLAGS_DESC "##INFO=<ID=PREDICTION_FLAGS,Number=.,Type=String,Description=\"Additional flags for the context in which the prediction happened.\">"

// The tags in the order of their header lines, and the header lines
#define INFO_TAGS { SPLICEAI_TAG, SPLICEAI_HAP_TAG, CS_TAG, RESCUE_TAG, RESCUE_PROB_TAG, RESCUE_TYPE_TAG, PTV_TAG, PTV_PROB_TAG, PTV_TYPE_TAG, PREDICTION_FLAGS_TAG }
#define INFO_TAG_DESCS { SPLICEAI_DESC, SPLICEAI_HAP_DESC, CS_DESC, RESCUE_DESC, RESCUE_PROB_DESC, RESCUE_TYPE_DESC, PTV_DESC, PTV_PROB_DESC, PTV_TYPE_DESC, PREDICTION_FLAGS_DESC }

// Appends the header lines of the output tags to the header. Returns 0, or -1 if a line could not be appended.
int add_rescueranger_tags(bcf_hdr_t *hdr);

// The value of INFO/RESCUE_TYPE for one transcript and phase: for every rescue type of the isoforms
// rescue|type|probability|affected_bases, joined with ':'. The caller frees the string.
kstring_t generate_rescue_type_tag(Isoforms *isoforms_struct);

// Appends the names of the set bits of an enum PredictionFlags value to concatenated_flags, in bit order, joined with ':'
void get_prediction_flags(uint32_t flags, kstring_t *concatenated_flags);

// One SpliceAI-style value per gene and phase: genes separated by ',', phases by '&'. `scores` is one of the record's
// native, diff or haplotype score arrays (genes * NUM_PHASES slots). A slot that was not scored gives dots.
void format_gene_scores(const Record *r, Score **scores, kstring_t *out);

// Adds the output INFO tags to every record, from the isoforms, scores and flags the earlier stages left on it. The records are
// written by the caller.
void write_tags(const Session *s, const Records *records);

#endif
