#ifndef HAPLOTYPE_HEADER
#define HAPLOTYPE_HEADER

#include <stdbool.h>
#include "klib/kstring.h"
#include "cpliceai/range.h"
#include "record.h"
#include "utils.h"

void get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap);

// Whether the variant of the record overlaps the gene: a base it changes lies in it (see variant_in_region)
bool record_overlaps_gene(const Record *r, const Gene *gene);

// The records of the phase that overlap the gene, which together make up the haplotype of the gene, in the order of the
// phase. hap_in holds them including r, and hap_ex, if given, holds them without r. With skip_other_unphased, records that are
// unphased do not belong to the haplotype of another record, so they only count for themselves.
// The caller frees hap_in, and hap_ex if given, with kv_destroy.
void find_haplotype_records(const Records *phase_records, const Record *r, const Gene *gene, const bool skip_other_unphased, Records *hap_in, Records *hap_ex);

// Appends the sequence of an isoform to out: the haplotype sequence of the gene without the non-coding regions that are
// spliced out. hap_seq starts at gene_start, and hap_ncrs holds the regions in transcript order, in haplotype coordinates.
// Regions that are not used (retained introns) stay in the sequence. Nothing is appended for an unknown strand.
// The sequence is always returned in ascending reference order, whatever the strand.
void splice_haplotype_sequence(const char *hap_seq, const NonCodingRegions *hap_ncrs, const int64_t gene_start, const uint32_t strand, kstring_t *out);

// The stretches of the haplotype sequence that splice_haplotype_sequence keeps, as offsets from gene_start, ascending, without
// the empty ones. The same walk over the regions: a used region that starts before the end of the previous one is skipped, and
// nothing is kept after the last used region. Empty for an unknown strand. The caller frees the result with kv_destroy.
Ranges isoform_segments(const NonCodingRegions *hap_ncrs, const int64_t gene_start, const uint32_t strand);

// The length of the isoform sequence the segments make
int64_t isoform_length(const Ranges *segments);

// A range in haplotype coordinates as a range of indices into the isoform sequence, as it is read in the direction of the
// transcript (on the reverse strand the sequence of splice_haplotype_sequence is reversed). The bases of the range that are not
// in the isoform drop out; a range with none is empty, at the place where its bases would have been. An empty or inverted
// range is empty.
Range hap_range_to_isoform(const Ranges *segments, const int64_t gene_start, const uint32_t strand, const Range hap_range);

// The NMD escape regions of the reference transcript as ranges of indices into the isoform sequence: every region goes through
// reference_range_to_haplotype (the haplotype records of the gene are hap_in; an insertion at a border goes with the base it is
// attached to in the direction of the transcript, which depends on the strand) and hap_range_to_isoform, and a region with no
// base in the isoform is left out. Every kind of region keeps its own list. out is initialised here and freed by the caller with
// NmdRegions_destroy.
void nmd_escape_ranges_in_isoform(const NmdRegions *regions, const Records *hap_in, const Range gene, const Ranges *segments, const uint32_t strand, NmdRegions *out);

// Whether position is in one of the half-open ranges
bool position_in_ranges(const Ranges *ranges, const int64_t position);

typedef kvec_t(NonCodingRegions) NonCodingRegionSets;

// Initialises out with every combination of one region from each non-empty group, keeping the regions of a combination
// in group order. The choice in the first group varies fastest. Empty groups are skipped, so no non-empty group gives no
// combinations. Regions are copied by value: their records vectors stay shared with the groups and are not freed by
// NonCodingRegionSets_destroy.
void enumerate_region_combinations(const NonCodingRegions *groups, const size_t n_groups, NonCodingRegionSets *out);
void NonCodingRegionSets_destroy(NonCodingRegionSets *sets);

// The offset of the first base of the codon that holds offset
int64_t codon_start(const int64_t offset);

// From the start of the codon that holds offset up to offset + length
Range codon_window(const int64_t offset, const int64_t length);

// codon_window for an offset that may lie just outside the coding sequence of sequence_length bases (a deletion of its first or
// last bases places the position at -1 or at the length): the offset is first moved into [0, sequence_length - 1] and the
// window ends no later than the sequence. An empty window for an empty sequence.
Range codon_window_within(const int64_t offset, const int64_t length, const int64_t sequence_length);

// Whether one of the codons starting at start, start + 3, ... below start + length is a stop codon. Every codon is read in
// full, also the last one, but only a codon that lies entirely inside the sequence (0 to sequence_length) is read: a scan that
// runs past the end, or starts outside the sequence, stops there.
bool has_stop_codon(const char *seq, const int64_t sequence_length, const int64_t start, const int64_t length);

// The same scan as has_stop_codon, giving the index of the first base of the first stop codon, or -1 if there is none. With
// ignore_final_codon a stop codon that ends exactly at the end of the sequence is not counted: for the sequence of an isoform
// that is the original stop codon, which is not a new stop.
int64_t first_stop_codon(const char *seq, const int64_t sequence_length, const int64_t start, const int64_t length, const bool ignore_final_codon);

// Where the stop codon of a PTV is, as the index of its first base in the isoform sequence (read in the direction of the
// transcript), or -1 if there is none. The original stop codon, the last three bases of the sequence, never counts. The
// position is where the change begins: its first changed base, the first base after a deletion, or the base an insertion is
// attached to. A position outside the sequence (a deletion of the last bases puts it at the length) is moved to the nearest base.
//
// ptv_stop_in_codons looks at the codons the change touches: changed_length bases from the position (for an insertion the
// attached base and the inserted bases; at least one base is always looked at, so a deletion looks at the codon of the junction).
// The window is measured from the position, not from the start of its codon, so a change that begins late in a codon and runs on
// into the next one is not cut short.
int64_t ptv_stop_in_codons(const char *seq, const int64_t length, const int64_t position, const int64_t changed_length);

// ptv_stop_downstream looks from the codon of the position to the end of the sequence, which is where a frameshift reads on. The
// codons before it are those of the reference, which has no stop before its own.
int64_t ptv_stop_downstream(const char *seq, const int64_t length, const int64_t position);

// The escape region a stop is in: NMD_START_PROXIMAL, NMD_LONG_EXON or NMD_DISTAL, tested in that order (the start-proximal
// stretch can overlap a long exon or the distal region), and 0 for a stop in none of them or for no stop (stop_position below 0)
uint32_t nmd_flag_for_stop(const int64_t stop_position, const NmdRegions *escape_ranges);

// Reverse complements seq in place. Only A, C, G and T are complemented; other characters only move.
void reverse_complement(char *seq, const size_t length);

// The prior probability of an isoform: over its used regions, the product of the prior of every donor and acceptor, or of
// 1 - prior for a site that is lost (posterior 0)
float isoform_prior(const NonCodingRegions *ncrs);


#endif
