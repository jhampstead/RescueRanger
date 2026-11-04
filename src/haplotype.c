#include "haplotype.h"
#include <cpliceai/range.h>

Ranges isoform_segments(const NonCodingRegions *hap_ncrs, const int64_t gene_start, const uint32_t strand) {
    Ranges segments;
    kv_init(segments);
    if (strand != STRAND_FWD && strand != STRAND_REV) return segments;

    const size_t n_regions = kv_size(*hap_ncrs);
    int64_t next = 0; // First base of the haplotype sequence that has not been handled yet
    for (size_t j = 0; j < n_regions; j++) {
        // The regions are walked in ascending order, which is reverse transcript order on the reverse strand
        const NonCodingRegion ncr = kv_A(*hap_ncrs, strand == STRAND_FWD ? j : n_regions - 1 - j);
        if (!ncr.is_used) continue; // A retained intron stays in the sequence

        const int64_t region_start = ncr.coords.start - gene_start;
        if (region_start < next) continue;

        if (region_start > next) kv_push(Range, segments, ((Range) { next, region_start }));
        next = ncr.coords.end - gene_start;
    }

    return segments;
}

void splice_haplotype_sequence(const char *hap_seq, const NonCodingRegions *hap_ncrs, const int64_t gene_start, const uint32_t strand, kstring_t *out) {
    Ranges segments = isoform_segments(hap_ncrs, gene_start, strand);

    for (size_t i = 0; i < kv_size(segments); i++) {
        const Range segment = kv_A(segments, i);
        kputsn(hap_seq + segment.start, segment.end - segment.start, out);
    }

    kv_destroy(segments);
}

int64_t isoform_length(const Ranges *segments) {
    int64_t length = 0;
    for (size_t i = 0; i < kv_size(*segments); i++) length += range_len(kv_A(*segments, i));

    return length;
}

// The number of bases of the isoform before the offset (from the gene start) of the haplotype sequence
static int64_t isoform_bases_before(const Ranges *segments, const int64_t offset) {
    int64_t kept = 0;
    for (size_t i = 0; i < kv_size(*segments); i++) {
        const Range segment = kv_A(*segments, i);
        if (offset <= segment.start) break;

        kept += (offset < segment.end ? offset : segment.end) - segment.start;
    }

    return kept;
}

Range hap_range_to_isoform(const Ranges *segments, const int64_t gene_start, const uint32_t strand, const Range hap_range) {
    const int64_t start = isoform_bases_before(segments, hap_range.start - gene_start);
    const int64_t end = hap_range.end <= hap_range.start ? start : isoform_bases_before(segments, hap_range.end - gene_start);

    if (strand == STRAND_REV) {
        const int64_t length = isoform_length(segments);
        return (Range) { length - end, length - start };
    }

    return (Range) { start, end };
}

void nmd_escape_ranges_in_isoform(const NmdRegions *regions, const Records *hap_in, const Range gene, const Ranges *segments, const uint32_t strand, NmdRegions *out) {
    kv_init(out->start_proximal);
    kv_init(out->long_exon);
    kv_init(out->distal);

    const Ranges *lists[3] = { &regions->start_proximal, &regions->long_exon, &regions->distal };
    Ranges *outs[3] = { &out->start_proximal, &out->long_exon, &out->distal };
    for (int l = 0; l < 3; l++) {
        for (size_t i = 0; i < kv_size(*lists[l]); i++) {
            const Range hap_range = reference_range_to_haplotype(hap_in, gene, kv_A(*lists[l], i), strand);
            const Range isoform_range = hap_range_to_isoform(segments, gene.start, strand, hap_range);
            if (isoform_range.end > isoform_range.start) kv_push(Range, *outs[l], isoform_range);
        }
    }
}

bool position_in_ranges(const Ranges *ranges, const int64_t position) {
    for (size_t i = 0; i < kv_size(*ranges); i++) {
        if (position >= kv_A(*ranges, i).start && position < kv_A(*ranges, i).end) return true;
    }

    return false;
}

void enumerate_region_combinations(const NonCodingRegions *groups, const size_t n_groups, NonCodingRegionSets *out) {
    kv_init(*out);

    for (size_t g = 0; g < n_groups; g++) {
        const NonCodingRegions group = groups[g];
        if (kv_size(group) == 0) continue;

        if (kv_size(*out) == 0) { // First non-empty group
            for (size_t j = 0; j < kv_size(group); j++) {
                NonCodingRegions set; kv_init(set);
                kv_push(NonCodingRegion, set, kv_A(group, j));
                kv_push(NonCodingRegions, *out, set);
            }
            continue;
        }

        // Copy the combinations so far once for every further region of the group, then extend copy j with region j
        const size_t prior = kv_size(*out);
        for (size_t j = 1; j < kv_size(group); j++) {
            for (size_t k = 0; k < prior; k++) {
                NonCodingRegions copy; kv_init(copy);
                kv_copy(NonCodingRegion, copy, kv_A(*out, k));
                kv_push(NonCodingRegions, *out, copy);
            }
        }

        for (size_t j = 0; j < kv_size(group); j++) {
            for (size_t k = 0; k < prior; k++) {
                kv_push(NonCodingRegion, kv_A(*out, j * prior + k), kv_A(group, j));
            }
        }
    }
}

void NonCodingRegionSets_destroy(NonCodingRegionSets *sets) {
    for (size_t i = 0; i < kv_size(*sets); i++) kv_destroy(kv_A(*sets, i));
    kv_destroy(*sets);
}

int64_t codon_start(const int64_t offset) {
    return (offset / 3) * 3;
}

Range codon_window(const int64_t offset, const int64_t length) {
    return (Range) { codon_start(offset), offset + length };
}

Range codon_window_within(const int64_t offset, const int64_t length, const int64_t sequence_length) {
    if (sequence_length <= 0) return (Range) { 0, 0 };

    int64_t inside = offset;
    if (inside >= sequence_length) inside = sequence_length - 1;
    if (inside < 0) inside = 0;

    Range window = codon_window(inside, length);
    if (window.end > sequence_length) window.end = sequence_length;

    return window;
}

int64_t first_stop_codon(const char *seq, const int64_t sequence_length, const int64_t start, const int64_t length, const bool ignore_final_codon) {
    static const char stop_codons[3][3] = { {'T', 'A', 'G'}, {'T', 'G', 'A'}, {'T', 'A', 'A'} };

    for (int64_t i = 0; i < length; i += 3) {
        if (start + i < 0 || start + i + 3 > sequence_length) continue;
        if (ignore_final_codon && start + i + 3 == sequence_length) continue;

        const char *codon = seq + start + i;
        for (int j = 0; j < 3; j++) {
            if (codon[0] == stop_codons[j][0] && codon[1] == stop_codons[j][1] && codon[2] == stop_codons[j][2]) return start + i;
        }
    }

    return -1;
}

bool has_stop_codon(const char *seq, const int64_t sequence_length, const int64_t start, const int64_t length) {
    return first_stop_codon(seq, sequence_length, start, length, false) >= 0;
}

// The position moved to the nearest base of a sequence of length bases (length above 0)
static int64_t clamp_to_sequence(const int64_t position, const int64_t length) {
    if (position >= length) return length - 1;
    if (position < 0) return 0;

    return position;
}

int64_t ptv_stop_in_codons(const char *seq, const int64_t length, const int64_t position, const int64_t changed_length) {
    if (length <= 0) return -1;

    const int64_t inside = clamp_to_sequence(position, length);
    const int64_t codon = codon_start(inside);
    const int64_t changed = changed_length > 0 ? changed_length : 1;

    return first_stop_codon(seq, length, codon, (inside - codon) + changed, true);
}

int64_t ptv_stop_downstream(const char *seq, const int64_t length, const int64_t position) {
    if (length <= 0) return -1;

    const int64_t codon = codon_start(clamp_to_sequence(position, length));

    return first_stop_codon(seq, length, codon, length - codon, true);
}

uint32_t nmd_flag_for_stop(const int64_t stop_position, const NmdRegions *escape_ranges) {
    if (stop_position < 0) return 0;

    if (position_in_ranges(&escape_ranges->start_proximal, stop_position)) return NMD_START_PROXIMAL;
    if (position_in_ranges(&escape_ranges->long_exon, stop_position)) return NMD_LONG_EXON;
    if (position_in_ranges(&escape_ranges->distal, stop_position)) return NMD_DISTAL;

    return 0;
}

static char complement(const char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return base;
    }
}

void reverse_complement(char *seq, const size_t length) {
    if (length == 0) return;

    for (size_t first = 0, last = length - 1; first < last; first++, last--) {
        const char tmp = complement(seq[first]);
        seq[first] = complement(seq[last]);
        seq[last] = tmp;
    }
    if (length % 2 == 1) seq[length / 2] = complement(seq[length / 2]);
}

float isoform_prior(const NonCodingRegions *ncrs) {
    float prior = 1.0;

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_used) continue; // A retained intron is not a splicing event

        if (ncr.donor != NULL) {
            float prob = ncr.donor->prior_probability;
            if (ncr.donor->posterior_probability == 0) prob = 1 - prob;
            prior *= prob;
        }
        if (ncr.acceptor != NULL) {
            float prob = ncr.acceptor->prior_probability;
            if (ncr.acceptor->posterior_probability == 0) prob = 1 - prob;
            prior *= prob;
        }
    }

    return prior;
}

bool record_overlaps_gene(const Record *r, const Gene *gene) {
    return variant_in_region(&r->alleles, gene->coords);
}

void find_haplotype_records(const Records *phase_records, const Record *r, const Gene *gene, const bool skip_other_unphased, Records *hap_in, Records *hap_ex) {
    kv_init(*hap_in);
    if (hap_ex != NULL) kv_init(*hap_ex);

    for (int j = 0; j < phase_records->n; j++) {
        Record *hap_r = kv_A(*phase_records, j);
        if (hap_r->genes == NULL) continue;
        if (skip_other_unphased && hap_r->unphased && hap_r != r) continue;

        for (int g = 0; g < hap_r->genes->n; g++) {
            Gene *hap_gene = kv_A(*hap_r->genes, g);
            if (hap_gene == NULL) continue;
            if (gene->gff_id != hap_gene->gff_id) continue;
            if (!variant_in_region(&hap_r->alleles, gene->coords)) continue;

            kv_push(Record *, *hap_in, hap_r);

            if (hap_ex != NULL && hap_r != r) kv_push(Record *, *hap_ex, hap_r);
        }
    }
}

void get_haplotype_sequence(const kstring_t reference_sequence, const Records *rs, const Range gene, kstring_t *out_hap) {
    kstring_t haplotype_sequence = {0};
    int64_t ref_idx = 0; // First base of the gene sequence that has not been handled yet

    for (int i = 0; i < rs->n; i++) {
        const Alleles *a = &kv_A(*rs, i)->alleles;
        if (!variant_in_region(a, gene)) continue;

        // The reference up to the change, which for an indel includes the base it is attached to
        const int64_t change_start = (a->pos > gene.start ? a->pos : gene.start) - gene.start;
        if (change_start > ref_idx) {
            kputsn(reference_sequence.s + ref_idx, change_start - ref_idx, &haplotype_sequence);
            ref_idx = change_start;
        }

        // The ALT bases from the gene start on. An ALT that runs past the gene end is kept whole.
        const int64_t first_alt = a->pos < gene.start ? gene.start - a->pos : 0;
        if (a->alt_len > first_alt) kputsn(a->alt + first_alt, a->alt_len - first_alt, &haplotype_sequence);

        // The REF bases the variant removes from the gene
        Range removed;
        if (a->ref_len > 0 && intersect_ranges((Range) { a->pos, a->pos + a->ref_len }, gene, &removed) != 0) {
            ref_idx += range_len(removed);
        }
    }

    const int64_t remainder = (int64_t) reference_sequence.l - ref_idx;
    if (remainder > 0) kputsn(reference_sequence.s + ref_idx, remainder, &haplotype_sequence);

    *out_hap = haplotype_sequence;
}
