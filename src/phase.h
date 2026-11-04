#ifndef PHASE_HEADER
#define PHASE_HEADER

#include <stdbool.h>
#include <htslib/vcf.h>

#include "record.h"

// The first two alleles of the genotype of the first sample, and whether they are phased
typedef struct Genotype {
    int allele[2];
    bool phased;
} Genotype;

// The name of the FORMAT field of a genotype tag: "FORMAT/GT" and "GT" are both GT. Returns NULL for a tag of another kind,
// such as "INFO/GT".
const char *Genotype_tag_name(const char *gt_tag);

// Reads the genotype of the first sample from the FORMAT field gt_tag (a name, see Genotype_tag_name). GT is read as htslib
// stores genotypes; another field is read as a string in the notation of GT ("0|1", "1/0"), of which the first two alleles are
// taken, phased when they are separated by '|'. Returns false when the record has no usable genotype: the field is not on the
// record, it has fewer than two alleles, or one of the two is missing ('.', as in "./.").
bool Genotype_read(const bcf1_t *v, const bcf_hdr_t *hdr, const char *gt_tag, Genotype *gt);

// Whether a record is present on the haplotype of the given phase (0 or 1) of the first sample, by the genotype in the FORMAT
// field gt_tag (see Genotype_read). A record without a usable genotype is on neither haplotype.
// Homozygous genotypes are present on both haplotypes, even homozygous reference ones. An unphased heterozygous genotype
// is present on both when include_unphased is true, and on neither when it is false. Otherwise it is whether the
// haplotype carries a non-reference allele.
int has_matching_phase(const bcf1_t *v, const bcf_hdr_t *hdr, const char *gt_tag, const int phase, const bool include_unphased);

// The records of a block that are on the haplotype of the phase, in the order of the block. Records without coding
// transcripts (tscripts is NULL) are left out, and so are records without a usable genotype in gt_tag, with a warning in
// phase 0. The caller frees the result with kv_destroy, the records stay owned by the block.
Records Record_select_for_phase(const Records *block_records, const bcf_hdr_t *hdr, const char *gt_tag, const int phase, const bool include_unphased);

// Removes, in place and keeping the order, the records that conflict with the variant in front of them: the record before
// it, or an earlier one that reaches further (also when that one was removed itself). By REF and ALT bases after the shared
// prefix, a variant is an insertion (no REF bases), a deletion (no ALT bases) or a substitution (SNV, MNV, complex). It
// conflicts when it is a substitution changing a base the variant in front changed, a deletion starting before the variant
// in front ends, or an insertion strictly inside a deletion or at the same place as another insertion. Logs every removal.
void filter_overlapping_variants(Records *vs, const bcf_hdr_t *hdr);

// Flags the focal record of a prediction, on every transcript, for the phase: SKIPPED_VARIANT_IN_PREDICTION when a record of its
// haplotype is too long for the prediction radius (not variant_fits_window; r->skipped is not read, 8b only sets it when it
// reaches the record), UNPHASED_VARIANT_IN_PREDICTION when one is unphased. Only records in the prediction window count
// (see variant_in_region: a base the variant changes lies in it, an insertion has both its neighbours in it), the focal record
// too. window is a genomic range. The haplotype must include the unphased records.
// A focal record without transcript entries or flag slots is left alone.
void flag_phasing_in_haplotype(Record *focal, const Records *haplotype, const Range window, const int prediction_radius, const int phase);

#endif
