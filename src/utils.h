#ifndef RESCUERANGER_UTILS_H
#define RESCUERANGER_UTILS_H

#include <limits.h>
#include <stdint.h>
#include <htslib/vcf.h>

#include <klib/kstring.h>
#include <klib/kvec.h>

#include "record.h"

// Bases of sequence context the model needs on each side of the prediction window.
#ifndef SPLICEAI_CONTEXT
#define SPLICEAI_CONTEXT 5000
#endif

// Largest prediction radius for which pred_array_len (the first derived size to overflow) fits in an int.
#define MAX_PREDICTION_RADIUS (((INT_MAX / NUM_SCORES) - 1) / 2)

// Half-open range [pos - offset, pos + offset + 1). Not clamped, so it can start below 0.
static inline Range range_around(const int64_t pos, const int offset) {
    return (Range) { pos - offset, pos + offset + 1 };
}

// Half-open range [pos - offset, pos), the left half of range_around without pos itself.
static inline Range range_before(const int64_t pos, const int offset) {
    return (Range) { pos - offset, pos };
}

// Number of positions in a prediction window of radius `prediction_radius`.
static inline int pred_width(const int prediction_radius) {
    return 2 * prediction_radius + 1;
}

// Length of the model input sequence for a prediction window of radius `prediction_radius`.
static inline int model_width(const int prediction_radius) {
    return pred_width(prediction_radius) + 2 * SPLICEAI_CONTEXT;
}

// Position within the model input sequence that corresponds to the variant's genomic position.
static inline int variant_position_in_sequence(const int prediction_radius) {
    return SPLICEAI_CONTEXT + prediction_radius;
}

// Number of floats in a prediction array (NUM_SCORES per position).
static inline int pred_array_len(const int prediction_radius) {
    return pred_width(prediction_radius) * NUM_SCORES;
}

Range variant_to_range(const bcf1_t *v, const int allele_index);

// True when neither the ALT bases nor the REF bases of the anchor-free view are more than the prediction radius. The shared
// anchor base of an indel does not count, so an indel of exactly the radius is kept.
bool variant_fits_window(const Alleles *a, const int prediction_radius);

// The bases a variant touches: its REF bases after the shared prefix and, for an indel (a shared prefix), the base before
// them that it is attached to. An insertion at p is the one base [p - 1, p). For a normalised record this is its VCF REF range.
Range variant_footprint(const Alleles *a);

// Where a variant counts when it is classified against exons, coding sequence and UTRs: its first changed base (the first
// base after the shared prefix), so a deletion counts at its first deleted base. An insertion changes no REF base and counts
// at the base it is attached to. Window centres and prediction matching do not use this: they stay on the attached base,
// the start of variant_footprint.
int64_t variant_position(const Alleles *a);

// Whether the variant touches a base of the range (its footprint shares a base with it)
bool variant_overlaps(const Alleles *a, const Range range);

// Whether a variant belongs to a region such as a gene: some REF base it changes lies in the region. The base an indel is
// attached to is not changed, so a deletion anchored on the last base of the region that only deletes bases after it is not in
// it (variant_overlaps says it is). An insertion lies between two bases and is in the region when both are, so an insertion at
// either edge (before the first base, after the last) is not.
// Use this for membership of a gene and of a gene's haplotype; windows and prediction matching use variant_overlaps.
bool variant_in_region(const Alleles *a, const Range region);

// The one primitive of the haplotype coordinates: how far one variant moves the coordinate b further down the haplotype, the
// bases it deletes before b minus the bases it inserts before b. The REF and ALT bases are paired one to one from the left and
// move nothing; the REF bases left over are deleted and count when they lie before b; the ALT bases left over are inserted at
// the point after the last REF base (between the bases point - 1 and point) and count when the point is before b, or at b when
// insertions_at_b is set (which side of a border an insertion exactly on it belongs to depends on what is being mapped). Nothing
// is clipped to a gene. allele_offset, haplotype_offset, haplotype_shift and region_to_haplotype are built on it.
int64_t variant_shift(const Alleles *a, const int64_t b, const bool insertions_at_b);

// How many bases shorter the haplotype is than the reference within bounds because of the variant, by its pieces: the REF and
// ALT bases are paired one to one from the left and move nothing, the REF bases left over are deleted and count when they lie in
// bounds, and the ALT bases left over are inserted after the last paired base and count when the base they are attached to (for
// a pure insertion the base before it) is in bounds. So a window edge through an MNV or a complex indel only counts the part of
// it that is in bounds, and a substitution is always 0. An empty or inverted range holds nothing.
int allele_offset(const Alleles *a, const Range bounds);

int haplotype_offset(const Records *rs, const Range bounds);

// How many bases further down the haplotype moves the coordinate b of a gene: the bases deleted before b minus the bases
// inserted before b, over the records that are in the gene (variant_in_region). The pieces of a variant are its paired REF and
// ALT bases (substitutions, which move nothing), the REF bases left over (a deletion, counted when before b) and the ALT bases
// left over (an insertion at the point after the last REF base, counted when before b). An insertion exactly at b counts only
// with insertions_at_b, because which side of a region border it belongs to depends on the region. Unlike allele_offset a
// variant that straddles b is split: only the part before b counts.
int64_t haplotype_shift(const Records *rs, const Range gene, const int64_t b, const bool insertions_at_b);

// Maps a region of an isoform (a UTR or an intron, in reference coordinates) to haplotype coordinates. An insertion exactly on
// the border of a region goes to the side that does not hold the feature that bounds the region, so it never splits it: it is
// in a UTR (start or stop codon stays whole) and outside an intron (the splice site is in the intron, the insertion is in the
// exon).
Range region_to_haplotype(const Records *rs, const Range gene, const Range region, const bool is_utr);

// A stretch of reference bases (not a region with a feature at its border, like region_to_haplotype takes) in haplotype
// coordinates. Deleted bases drop out of it and inserted bases inside it are part of it. An insertion at a border belongs to
// the side that holds the base it is attached to, which is the base before it in the direction of the transcript: on the forward
// strand (and an unknown one) the lower neighbour, so at the start border it is outside the range and at the end border inside;
// on the reverse strand the higher neighbour, so at the start border it is inside and at the end border outside. An empty range
// stays empty.
Range reference_range_to_haplotype(const Records *rs, const Range gene, const Range range, const uint32_t strand);

// The haplotype offset in the prediction window of a variant, which is centred on its first base (for an indel the base it
// is attached to, the VCF position of a normalised indel)
int prediction_offset(const Alleles *a, const Records *rs, const int prediction_radius);

// Raises each of the NUM_SCORES values of best to the highest value of its kind among the length floats of insert_start, which
// are consecutive entries of NUM_SCORES scores
void get_best_insertion_score(const float *insert_start, const int length, float *best);

float *match_predictions_to_reference(const float *predictions, const Range pred_range, const Records *rs, const size_t num_reference_predictions);

#endif
