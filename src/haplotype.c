#include "haplotype.h"
#include <cpliceai/range.h>

void get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap) {
    const int sequence_offset = haplotype_offset(rs, gene);

    // Prepare kstring
    kstring_t haplotype_sequence = {0};
    // haplotype_sequence.l = reference_sequence.l - sequence_offset;
    // haplotype_sequence.m = haplotype_sequence.l + 1; 
    // haplotype_sequence.s = malloc(haplotype_sequence.m);

    int ref_idx = 0;
    int hap_idx = 0;

    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);

        // Get reference position, intersected to gene boundaries to manage start overlap
        Range ref = variant_to_range(r->v, 0);
        Range intersected_ref = {0}; 
        if (intersect_ranges(ref, gene, &intersected_ref) == 0) continue;
        shift_range(&ref, -gene.start);
        shift_range(&intersected_ref, -gene.start);
        int ref_segment = intersected_ref.start - ref_idx;

        // Get alt position, intersected to gene boundaries
        Range alt = variant_to_range(r->v, 1);
        Range intersected_alt = {0}; 
        int alt_segment = 0;
        if (intersect_ranges(alt, gene, &intersected_alt) && range_len(intersected_alt) > 0) {
            shift_range(&alt, -gene.start);
            shift_range(&intersected_alt, -gene.start);

            alt_segment = range_len(alt);
        }

        char *ref_seq_pos = reference_sequence.s + ref_idx;
        int skipped_alt_bases = (intersected_alt.start - alt.start);
        char *alt_seq_pos = r->v->d.allele[1] + skipped_alt_bases;
        if (alt_segment > 0) alt_segment -= skipped_alt_bases;

        // Offset start and segment lengths for the leading reference base when dealing with indels
        // int is_indel = ((ref.start == intersected_ref.start & range_len(ref) > 1) | (alt.start == intersected_alt.start & range_len(alt) > 1));
        int skipped_ref_base = (range_len(ref) != range_len(alt)) & (ref.start == intersected_ref.start);

        // Copy over reference sequence without any variation
        // if (ref_segment + is_indel > 0) strncpy(haplotype_sequence.s + hap_idx, ref_seq_pos, ref_segment + is_indel);
        if (ref_segment + skipped_ref_base > 0) kputsn(ref_seq_pos, ref_segment + skipped_ref_base, &haplotype_sequence);

        ref_idx += ref_segment + skipped_ref_base;
        hap_idx += ref_segment + skipped_ref_base;

        // Copy over variant
        // if (alt_segment - is_indel > 0) strncpy(haplotype_sequence.s + hap_idx, alt_seq_pos + is_indel, alt_segment - is_indel);
        if (alt_segment - skipped_ref_base > 0) kputsn(alt_seq_pos + skipped_ref_base, alt_segment - skipped_ref_base, &haplotype_sequence);

        ref_idx += range_len(intersected_ref) - skipped_ref_base;
        hap_idx += alt_segment - skipped_ref_base;
    }

    // Copy over remaining sequence from reference and complete sequence string
    int remainder = reference_sequence.l - ref_idx;
    if (remainder > 0) {
        haplotype_sequence.l = hap_idx;
        kputsn(reference_sequence.s + ref_idx, remainder, &haplotype_sequence);
    }
    // if (remainder > 0) strncpy(haplotype_sequence.s + hap_idx, reference_sequence.s + ref_idx, remainder);
    // haplotype_sequence.s[haplotype_sequence.l] = '\0';
 
    *out_hap = haplotype_sequence;
}
