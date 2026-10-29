#ifndef REFERENCE_SCORES_HEADER
#define REFERENCE_SCORES_HEADER

#include <stdint.h>
#include <htslib/faidx.h>

// The precomputed reference predictions of cpliceai (the file cpliceai_reference writes), as an alternative to predicting the
// reference sequence of every window with the model. They are optional: wherever they are not there, or do not cover the gene, the
// caller predicts the sequence as it always has, so the program runs without them.
//
// There is no gene regions file to give. The scores are looked up by the gene of the annotation (its contig, name, coordinates and
// strand), and a gene that is not in the scores exactly like that is not read from them. That check is what keeps scores made
// from another annotation from being used, and it covers every window that is read, which a fingerprint of a gene set would not.
//
// This header and its source are the only place that includes the gene and reference headers of cpliceai. Those define a Gene of
// their own, which record.h's Gene would clash with, so nothing here is in terms of either.
typedef struct ReferenceScores ReferenceScores;

// Whether this build has them: the cpliceai it is built against has to offer them. A build without them refuses the option.
int ReferenceScores_supported(void);

// Opens the reference scores and checks them against the FASTA (by contig names and lengths) they were made from: scores for
// another genome would give wrong predictions with no sign of it. Returns NULL, after logging why, if they cannot be used. The
// FASTA index must outlive the result.
ReferenceScores *ReferenceScores_open(const char *reference_path, faidx_t *fai);

// Logs how many windows were read and how many genes were not in the scores, then frees everything
void ReferenceScores_close(ReferenceScores *scores);

// The reference predictions of the 2 * radius + 1 bases around the 0-based position pos: NUM_SCORES floats per base, in the order
// cpliceai predicts them (none of the two, acceptor, donor), ascending whatever the strand, which the caller frees. Positions
// outside the gene have no splice site (none of the two is 1).
//
// The gene is the annotation's: its name, 0-based half-open coordinates and strand ('+' or '-') must be those of a gene in the
// reference scores, and if they are not the result is NULL, with a warning the first time for that gene (it says what the scores
// have under that name, if anything, so that a difference in how the coordinates are made shows at once). The caller then predicts
// the sequence. The scores of one gene are kept after the call, so calls for the same gene one after the other do not read it again.
float *ReferenceScores_window(ReferenceScores *scores, const char *chrom, const char *gene_name, int64_t gene_start, int64_t gene_end,
                              char strand, int64_t pos, int radius);

// How many windows were given, and how many different genes were refused, so far
void ReferenceScores_counts(const ReferenceScores *scores, long *windows, long *refused_genes);

#endif
