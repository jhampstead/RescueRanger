#include "phase.h"

#include <stdlib.h>
#include <cpliceai/range.h>

#include "logging/log.h"
#include "utils.h"

int has_matching_phase(const bcf1_t *v, const bcf_hdr_t *hdr, const int phase, const bool include_unphased) {
    int *gt = NULL; int ngts = 0;
    bcf_get_genotypes(hdr, (bcf1_t *) v, &gt, &ngts);

    // Homozygous
    if (bcf_gt_allele(gt[0]) == bcf_gt_allele(gt[1])) {
        free(gt);
        return 1;
    }

    // Unphased
    if (!bcf_gt_is_phased(gt[0]) & !bcf_gt_is_phased(gt[1])) {
        free(gt);
        return include_unphased;
    }

    int is_present = bcf_gt_allele(gt[phase]);
    free(gt);
    return is_present;
}

Records Record_select_for_phase(const Records *block_records, const bcf_hdr_t *hdr, const int phase, const bool include_unphased) {
    Records selected;
    kv_init(selected);

    for (int ri = 0; ri < block_records->n; ri++) {
        Record *r = kv_A(*block_records, ri);

        if (r->tscripts == NULL) continue;
        if (!has_matching_phase(r->v, hdr, phase, include_unphased)) continue;

        kv_push(Record *, selected, r);
    }

    return selected;
}

// What a variant does to the reference, by its REF and ALT bases after the shared prefix
enum ChangeKind { SUBSTITUTION, INSERTION, DELETION };

void filter_overlapping_variants(Records *vs, const bcf_hdr_t *hdr) {
    Range last = {0};
    enum ChangeKind last_kind = SUBSTITUTION;

    Range current = {0};
    enum ChangeKind current_kind = SUBSTITUTION;

    size_t n_kept = 0;
    for (size_t i = 0; i < kv_size(*vs); i++) {
        Record *r = kv_A(*vs, i);
        const Alleles *a = &r->alleles;

        // The previous variant becomes the one in front if it reaches further, also when it was removed
        if (current.end > last.end) {
            last = current;
            last_kind = current_kind;
        }

        current = (Range) { a->pos, a->pos + a->ref_len }; // The REF bases it changes, empty for an insertion
        current_kind = a->ref_len == 0 ? INSERTION : a->alt_len == 0 ? DELETION : SUBSTITUTION;

        bool conflicts;
        if (current_kind == INSERTION) {
            conflicts = (last_kind == DELETION && last.start < current.start && current.end < last.end)
                || (last_kind == INSERTION && last.start == current.start && last.end == current.end);
        } else if (current_kind == DELETION) {
            conflicts = current.start < last.end;
        } else {
            Range intersection;
            conflicts = intersect_ranges(current, last, &intersection) && range_len(intersection) > 0;
        }

        if (conflicts) {
            log_warn("Skipping on %s:%li, overlaps previous variant", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos + 1);
            continue;
        }

        kv_A(*vs, n_kept++) = r;
    }

    vs->n = n_kept;
}

void flag_phasing_in_haplotype(Record *focal, const Records *haplotype, const Range window, const int prediction_radius, const int phase) {
    if (focal->transcript_genes == NULL || focal->flags == NULL) return;

    bool has_skipped = false, has_unphased = false;
    for (size_t i = 0; i < kv_size(*haplotype); i++) {
        const Record *r = kv_A(*haplotype, i);
        if (!variant_in_region(&r->alleles, window)) continue;

        if (!variant_fits_window(&r->alleles, prediction_radius)) has_skipped = true;
        if (r->unphased) has_unphased = true;
    }

    for (size_t t = 0; t < kv_size(*focal->transcript_genes); t++) {
        if (has_skipped) set_record_flag(focal, (int) t, phase, SKIPPED_VARIANT_IN_PREDICTION);
        if (has_unphased) set_record_flag(focal, (int) t, phase, UNPHASED_VARIANT_IN_PREDICTION);
    }
}
