#include "isoform_position.h"

#include <stdbool.h>
#include <stdlib.h>

#include "logging/log.h"
#include "utils.h"

static int64_t overlap_len(const Range a, const Range b) {
    const int64_t start = a.start > b.start ? a.start : b.start;
    const int64_t end = a.end < b.end ? a.end : b.end;

    return end > start ? end - start : 0;
}

// The coding bases of [start, end): the bases in the gene that are in no used non-coding region
static int64_t coding_bases(const int64_t start, const int64_t end, const NonCodingRegions *ncrs, const Range gene) {
    const Range bases = { start, end };
    int64_t n = overlap_len(bases, gene);

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_used) continue; // A retained intron is part of the isoform's sequence

        n -= overlap_len(bases, ncr.coords);
    }

    return n;
}

// Whether an insertion between the bases point - 1 and point is in the coding sequence
static bool coding_insertion_point(const int64_t point, const NonCodingRegions *ncrs, const Range gene) {
    if (point <= gene.start || point >= gene.end) return false;

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_used) continue;

        // Next to a UTR the start or stop codon stays whole, so the insertion is in the UTR. At an intron border the splice
        // site is in the intron, so the insertion is in the exon.
        if (ncr.is_utr && point >= ncr.coords.start && point <= ncr.coords.end) return false;
        if (!ncr.is_utr && point > ncr.coords.start && point < ncr.coords.end) return false;
    }

    return true;
}

CodingEffect variant_coding_effect(const Alleles *a, const NonCodingRegions *ncrs, const Range gene) {
    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;

    CodingEffect effect = {0};
    effect.substituted = coding_bases(a->pos, a->pos + paired, ncrs, gene);
    if (a->ref_len > paired) effect.deleted = coding_bases(a->pos + paired, a->pos + a->ref_len, ncrs, gene);
    if (a->alt_len > paired && coding_insertion_point(a->pos + paired, ncrs, gene)) effect.inserted = a->alt_len - paired;
    effect.coding = effect.substituted > 0 || effect.deleted > 0 || effect.inserted > 0;

    return effect;
}

int64_t isoform_coding_index(const Records *rs, const NonCodingRegions *ncrs, const Range gene, const int strand, const int64_t x) {
    const bool forward = strand == STRAND_FWD;

    // The coding bases of the reference before x in reading direction
    int64_t index = forward ? coding_bases(gene.start, x, ncrs, gene) : coding_bases(x + 1, gene.end, ncrs, gene);

    for (size_t i = 0; i < kv_size(*rs); i++) {
        const Alleles *a = &kv_A(*rs, i)->alleles;
        if (!variant_in_region(a, gene)) continue;

        // The REF bases left over after pairing are deleted: the coding ones before x leave the sequence
        const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
        const int64_t deleted_start = a->pos + paired;
        const int64_t deleted_end = a->pos + a->ref_len;
        if (deleted_end > deleted_start) {
            index -= forward ? coding_bases(deleted_start, deleted_end < x ? deleted_end : x, ncrs, gene)
                             : coding_bases(deleted_start > x + 1 ? deleted_start : x + 1, deleted_end, ncrs, gene);
        }

        // The ALT bases left over are inserted at the point after the last REF base, and are coding at a coding point
        if (a->alt_len > a->ref_len) {
            const int64_t point = a->pos + a->ref_len;
            const bool before = forward ? point <= x : point > x;
            if (before && coding_insertion_point(point, ncrs, gene)) index += a->alt_len - a->ref_len;
        }
    }

    return index;
}

// The used non-coding region that contains the base b, or NULL
static const NonCodingRegion *region_containing(const int64_t b, const NonCodingRegions *ncrs) {
    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion *ncr = &kv_A(*ncrs, i);
        if (ncr->is_used && b >= ncr->coords.start && b < ncr->coords.end) return ncr;
    }

    return NULL;
}

// The first coding base at or after b, or -1 when there is none in the gene
static int64_t next_coding_base(int64_t b, const NonCodingRegions *ncrs, const Range gene) {
    if (b < gene.start) b = gene.start;

    while (b < gene.end) {
        const NonCodingRegion *ncr = region_containing(b, ncrs);
        if (ncr == NULL) return b;

        b = ncr->coords.end;
    }

    return -1;
}

// The last coding base at or before b, or -1 when there is none in the gene
static int64_t previous_coding_base(int64_t b, const NonCodingRegions *ncrs, const Range gene) {
    if (b >= gene.end) b = gene.end - 1;

    while (b >= gene.start) {
        const NonCodingRegion *ncr = region_containing(b, ncrs);
        if (ncr == NULL) return b;

        b = ncr->coords.start - 1;
    }

    return -1;
}

// The coding base a variant with a coding effect is placed on: its first changed coding base in genomic order. When no REF base
// it changes is coding, only inserted bases are, and it sits on the nearest coding base to the left of the insertion point.
static int64_t first_changed_coding_base(const Alleles *a, const NonCodingRegions *ncrs, const Range gene) {
    if (a->ref_len > 0) {
        const int64_t base = next_coding_base(a->pos, ncrs, gene);
        if (base >= 0 && base < a->pos + a->ref_len) return base;
    }

    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
    const int64_t left = previous_coding_base(a->pos + paired - 1, ncrs, gene);

    return left >= 0 ? left : a->pos;
}

int get_range_in_isoform(const Alleles *a, const NonCodingRegions *ncrs, const Records *rs, const int strand, const Range gene, int64_t *out_position, int *out_offset) {
    const CodingEffect effect = variant_coding_effect(a, ncrs, gene);
    if (!effect.coding) {
        *out_position = variant_position(a);
        *out_offset = 0;
        return 1;
    }

    *out_offset = effect.deleted - effect.inserted;

    // The index of the base the variant is placed on. A substituted base, or the base an insertion is attached to, has its own
    // index. A deleted base has the index of the base that follows it in reading direction; this convention places the variant
    // on the first remaining base at or after it in genomic order, which on the reverse strand is the base before it in reading
    // direction, so one less.
    const int64_t first = first_changed_coding_base(a, ncrs, gene);
    int64_t index = isoform_coding_index(rs, ncrs, gene, strand, first);

    const int paired = a->ref_len < a->alt_len ? a->ref_len : a->alt_len;
    const bool deleted = first >= a->pos + paired && first < a->pos + a->ref_len;
    if (strand != STRAND_FWD && deleted) index--;

    *out_position = index;

    return 0;
}

bool isoform_base_is_coding(const NonCodingRegions *ncrs, const int64_t b) {
    return region_containing(b, ncrs) == NULL;
}

bool variant_in_cds_exon(const Alleles *a, const Exon *exon) {
    if (exon->cds == NULL) return false;
    const Range cds = *exon->cds;

    if (a->ref_len > 0) {
        const int64_t first = variant_position(a);
        return first >= cds.start && first < cds.end;
    }

    // An insertion lies between the bases pos - 1 and pos
    const int64_t point = a->pos;
    if (point > cds.start && point < cds.end) return true;
    if (point == cds.start) return exon->ref_bounds != NULL && exon->ref_bounds->start == cds.start;
    if (point == cds.end) return exon->ref_bounds != NULL && exon->ref_bounds->end == cds.end;

    return false;
}

enum UTRType is_utr_variant(const Alleles *a, const NonCodingRegions *ncrs) {
    // The UTRs are reference ranges, so the variant is compared with them at its reference position
    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        if (!ncr.is_utr) continue;

        if (a->ref_len > 0) {
            // A base the variant changes lies in the UTR. ranges_overlap is closed on both ends, so it would also match a
            // variant that only touches the UTR.
            if (a->pos < ncr.coords.end && a->pos + a->ref_len > ncr.coords.start) return ncr.is_utr;
        } else {
            // An insertion lies between two bases, and in the UTR when that point is in it or at either edge: at the CDS
            // side it does not split the start or stop codon (as in get_range_in_isoform)
            if (a->pos >= ncr.coords.start && a->pos <= ncr.coords.end) return ncr.is_utr;
        }
    }

    return 0;
}

void unique_isoform_shifts(const IsoformVariants *shifts, IsoformVariants *unique) {
    kv_init(*unique);

    for (size_t i = 0; i < kv_size(*shifts); i++) {
        const IsoformVariant shift = kv_A(*shifts, i);

        bool already_present = false;
        for (size_t j = 0; j < kv_size(*unique); j++) {
            const IsoformVariant kept = kv_A(*unique, j);
            if (shift.shift_type == kept.shift_type && shift.seq_pos == kept.seq_pos && shift.modified_bases == kept.modified_bases) {
                already_present = true;
            }
        }

        if (!already_present) kv_push(IsoformVariant, *unique, shift);
    }
}
