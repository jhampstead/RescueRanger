#ifndef ISOFORM_POSITION_HEADER
#define ISOFORM_POSITION_HEADER

#include <stdint.h>
#include <cpliceai/range.h>

#include "record.h"

// Where variants and positions of a gene lie in an isoform: its sequence without the non-coding regions that are spliced
// out. ncrs are the regions of the isoform in transcript order, in reference coordinates; regions that are not used
// (retained introns) are skipped. gene is the half-open gene range.

// Places a variant in the isoform, by what it does to the coding sequence (variant_coding_effect: the REF and ALT bases are
// paired from the left, the coding part of the variant counts and the non-coding part does not, on both strands). Returns 0
// when some part of the variant is coding and 1 otherwise. For a coding variant out_offset is how many bases shorter it makes
// the coding sequence (coding bases deleted minus bases inserted into it), and out_position is the index into the coding
// sequence (in the reading direction of the strand) of the base it is placed on: its first changed coding base in genomic order
// (for an insertion, the coding base it is attached to), see isoform_coding_index. A deleted base takes the index of the first
// remaining base at or after it in genomic order. A deletion that removes the last coding bases in genomic order has no such
// base and is placed just outside the sequence: at its length on the forward strand and at -1 on the reverse strand. The caller
// has to keep an index it reads from leaving the sequence. rs are the records of the haplotype of the gene, including the
// variant. For a non-coding variant out_position is its reference position (variant_position) and out_offset is 0.
int get_range_in_isoform(const Alleles *a, const NonCodingRegions *ncrs, const Records *rs, const int strand, const Range gene, int64_t *out_position, int *out_offset);

// Whether the base b is coding in the isoform: in no used non-coding region
bool isoform_base_is_coding(const NonCodingRegions *ncrs, const int64_t b);

// What a variant does to the coding sequence of an isoform. The REF and ALT bases are paired one to one from the left in
// genomic coordinates, on either strand. The paired bases are substitutions; the REF bases left over are a deletion at the
// right-hand end of the block, and the ALT bases left over an insertion after the last paired base. Only the coding part of
// the variant counts: a base is coding when it is in the gene and in no used non-coding region (bases outside the gene,
// UTRs and spliced introns are not), and an insertion point is coding strictly inside the CDS and at an intron border, and
// not at a UTR border or at the edge of the gene.
typedef struct CodingEffect {
    bool coding;     // some paired base, deleted base or the insertion point is coding
    int substituted; // coding bases that are replaced one for one
    int deleted;     // coding bases the variant deletes
    int inserted;    // bases the variant inserts in the coding sequence
} CodingEffect;

// ncrs are the non-coding regions of the isoform in reference coordinates, gene the half-open gene range. The change of the
// coding sequence length is deleted - inserted.
CodingEffect variant_coding_effect(const Alleles *a, const NonCodingRegions *ncrs, const Range gene);

// The index of the reference base x in the coding sequence of the isoform of a haplotype: the number of coding bases that come
// before x in the reading direction of the strand (forward: the coding bases below x, reverse: the ones above it). Only the
// coding part of the variants counts: the coding bases they delete and the bases they insert at coding insertion points, as
// variant_coding_effect splits them, so a variant inside a spliced-out intron or UTR changes nothing. For a base the haplotype
// deletes it is the index of the base that follows it in reading direction. rs are the records of the haplotype of the gene
// (only the ones in the gene count), ncrs the non-coding regions of the isoform in reference coordinates.
int64_t isoform_coding_index(const Records *rs, const NonCodingRegions *ncrs, const Range gene, const int strand, const int64_t x);

// Whether a variant counts as being in the coding sequence of an exon, for the per-exon flags. A variant with REF bases counts
// at its first changed base. An insertion lies between two bases: inside the CDS it is in; at an edge of the CDS it is in only
// where the CDS ends at the end of the exon (an intron border, the splice site is in the intron) and not where the exon goes on
// as a UTR (the start or stop codon stays whole). False for an exon without CDS.
bool variant_in_cds_exon(const Alleles *a, const Exon *exon);

// The UTR (UTR5 or UTR3) a variant is in, or 0. A variant with REF bases is in a UTR when a base it changes lies in it (the
// attached base of an indel does not count, and a variant that only touches the UTR is not in it). An insertion is in the UTR
// when its point, between two bases, is in it or at either edge. ncrs are in reference coordinates and so is the variant: the
// other variants of the haplotype do not matter (the check used to move the variant into haplotype coordinates first, which put
// a variant after a long deletion in the wrong place).
enum UTRType is_utr_variant(const Alleles *a, const NonCodingRegions *ncrs);

enum ShiftType {
    INDEL = 1,
    SPLICE = 2
};

// A change of the length of an isoform, placed on it: what an indel of a record of the haplotype or a changed splice site does
typedef struct IsoformFrameshift {
    int coding;
    int partial;
    enum ShiftType shift_type;
    int64_t seq_pos; // Signed: a partial indel can sit before the first coding base
    int shift;
    int modified_bases;
    Record *r;
} IsoformVariant;

typedef kvec_t(IsoformVariant) IsoformVariants;

// Copies the shifts into unique, in order, leaving out every shift that has the shift_type, seq_pos and modified_bases of one
// copied before it. The record and the shift itself are not compared, so of two records that changed the same site the first one
// is kept. unique is initialised here and freed by the caller with kv_destroy.
void unique_isoform_shifts(const IsoformVariants *shifts, IsoformVariants *unique);

#endif
