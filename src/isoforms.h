#ifndef ISOFORMS_HEADER
#define ISOFORMS_HEADER

#include <stdbool.h>
#include <htslib/vcf.h>

#include "haplotype.h"
#include "isoform_position.h"
#include "record.h"
#include "session.h"

// A frameshift within this many bases of another one can be restored by it, and an in-frame shift of this many bases or more is
// too long to count as a rescue. Both are in bases of the coding sequence.
#define MAX_FRAMESHIFT_RADIUS 80
#define MAX_SPLICE_RADIUS 80

// What all isoforms of one transcript of a record share
typedef struct IsoformContext {
    const bcf_hdr_t *hdr;
    int phase;
    Record *r; // The record the isoforms are built for
    int ti; // The index of the transcript in the record
    const Transcript *t;
    const Gene *gene;
    const TscriptCsqs *csq; // The consequences of the record for the transcript
    const Records *hap_in; // The records that make up the haplotype of the gene, including r
    const kstring_t *hap_seq; // The haplotype sequence of the gene
} IsoformContext;

// The rescue probability, the rescue types with their probabilities and the rescue call of an Isoforms struct, from its isoforms
// (whose posterior probabilities are set). Isoforms with the same rescue type and the same affected bases are one entry. The
// struct needs m_rescue_types >= 1 and n_rescue_types 0; Isoforms_destroy frees what this allocates.
void calculate_rescue_probabilities(Isoforms *isoforms_struct);

// The probability of a PTV, the PTV call (more than half of the probability) and its type, from the isoforms with a PTV. The type
// is stop gained or frameshift, whichever has the larger probability. When the two are equal it is a long in-frame indel if that
// has more probability than they do, and none otherwise.
void calculate_ptv_probabilities(Isoforms *isoforms_struct);

// The rescue type of a shift that keeps the frame: a negative shift adds bases (intron retention), a positive one removes them
// (exon skipping), partially when the shift does not come from a competitor site
enum RescueType shift_rescue_type(const IsoformVariant *shift);

// Builds the non-coding regions of a transcript on the haplotype, as groups in transcription order: the 5' UTR, one group per
// intron between two coding exons with every way it can be spliced, and the 3' UTR. An isoform takes one region from every
// group. A transcript on an unknown strand gets no groups. Returns false if the transcript has no coding exon. The record is the
// one the isoforms are built for; it is put in the UTR regions. A lost splice site retains the intron: a region that is not used.
// The regions point at the splice sites of the exons, so they are only valid while the exons have them.
// The caller frees groups with free_noncoding_region_groups, also when false is returned.
bool build_noncoding_regions(const Transcript *t, const Gene *gene, Record *r, NonCodingRegionSets *groups);

// Frees the groups of build_noncoding_regions, with the record vectors of the UTRs
void free_noncoding_region_groups(NonCodingRegionSets *groups);

// Appends the splice shifts of the regions of an isoform. A region whose donor or acceptor is lost is an intron retention at
// the canonical donor base, for the records that changed the lost site or sites. Otherwise the donor and the acceptor each
// give the shift from their canonical position, for the records that changed them.
void collect_splice_shifts(const NonCodingRegions *ncrs, const Records *hap_in, const Transcript *t, const Gene *gene, IsoformVariants *out);

// The shifts of an isoform, sorted in transcription order: the indel of every record of the haplotype and the splice shifts of
// the regions, keeping only those that change the isoform (a coding indel, or a splice shift other than 0). Also flags the
// records that lie in a UTR of the isoform.
void collect_isoform_shifts(const IsoformContext *ctx, const NonCodingRegions *ncrs, IsoformVariants *shifts);

// The single NMD flag of a transcript and phase: the nmd_flag of the isoform with the highest posterior probability among those
// with nmd_ptv set (the first of them on a tie), or 0 if none has it
uint32_t nmd_flag_of_isoforms(const Isoforms *isoforms);

// Builds one isoform of the transcript from one non-coding region per group: splices the haplotype sequence, places the shifts
// of the haplotype on it, and decides whether the record is rescued on this isoform and whether the isoform has a new PTV.
// Sets the context flags of the records on the way. The NMD flag is kept on the isoform (nmd_ptv and nmd_flag), not set on the
// record. The posterior probability is left 0.
Isoform evaluate_isoform(const IsoformContext *ctx, const NonCodingRegions *ncrs);

// Stage 3 of a phase: for every record and every transcript the record has a splice, frameshift or stop trigger for, builds
// every isoform the haplotype can give and stores them in the record, with their posterior, rescue and PTV probabilities.
void build_isoforms_for_phase(const Session *s, const int phase, const Records *records);

#endif
