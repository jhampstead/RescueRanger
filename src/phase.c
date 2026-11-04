#include "phase.h"

#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <cpliceai/range.h>

#include "logging/log.h"
#include "utils.h"

const char *Genotype_tag_name(const char *gt_tag) {
    if (strncmp(gt_tag, "FORMAT/", 7) == 0) return gt_tag + 7;
    if (strchr(gt_tag, '/') != NULL) return NULL;
    return gt_tag;
}

// GT, as htslib stores it: allele numbers with a phase bit, the phase on the second allele
static bool read_gt(const bcf1_t *v, const bcf_hdr_t *hdr, Genotype *out) {
    int32_t *gt = NULL; int ngts = 0;
    const int n = bcf_get_genotypes(hdr, (bcf1_t *) v, &gt, &ngts);
    const int nsamples = bcf_hdr_nsamples(hdr);

    bool usable = n > 0 && nsamples > 0 && n / nsamples >= 2
        && gt[1] != bcf_int32_vector_end
        && !bcf_gt_is_missing(gt[0]) && !bcf_gt_is_missing(gt[1]);
    if (usable) {
        out->allele[0] = bcf_gt_allele(gt[0]);
        out->allele[1] = bcf_gt_allele(gt[1]);
        out->phased = bcf_gt_is_phased(gt[0]) || bcf_gt_is_phased(gt[1]);
    }

    free(gt);
    return usable;
}

// One allele number at *s, which is moved past it. A missing allele ('.') is not one.
static bool parse_allele(const char **s, int *allele) {
    if (!isdigit((unsigned char) **s)) return false;

    char *end;
    const long value = strtol(*s, &end, 10);
    if (value > INT_MAX) return false;

    *allele = (int) value;
    *s = end;
    return true;
}

// A genotype written as a string in the notation of GT: the first two alleles, the separator between them saying the phase
static bool parse_genotype_string(const char *s, Genotype *out) {
    if (!parse_allele(&s, &out->allele[0])) return false;
    if (*s != '|' && *s != '/') return false;
    out->phased = *s == '|';
    s++;
    if (!parse_allele(&s, &out->allele[1])) return false;
    return *s == '\0' || *s == '|' || *s == '/';
}

bool Genotype_read(const bcf1_t *v, const bcf_hdr_t *hdr, const char *gt_tag, Genotype *gt) {
    if (strcmp(gt_tag, "GT") == 0) return read_gt(v, hdr, gt);

    char **values = NULL; int nvalues = 0;
    const int n = bcf_get_format_string(hdr, (bcf1_t *) v, gt_tag, &values, &nvalues);
    const bool usable = n > 0 && bcf_hdr_nsamples(hdr) > 0 && parse_genotype_string(values[0], gt);

    if (values) {
        free(values[0]);
        free(values);
    }
    return usable;
}

int has_matching_phase(const bcf1_t *v, const bcf_hdr_t *hdr, const char *gt_tag, const int phase, const bool include_unphased) {
    Genotype gt;
    if (!Genotype_read(v, hdr, gt_tag, &gt)) return 0;

    // Homozygous
    if (gt.allele[0] == gt.allele[1]) return 1;

    // Unphased
    if (!gt.phased) return include_unphased;

    return gt.allele[phase];
}

Records Record_select_for_phase(const Records *block_records, const bcf_hdr_t *hdr, const char *gt_tag, const int phase, const bool include_unphased) {
    Records selected;
    kv_init(selected);

    for (int ri = 0; ri < block_records->n; ri++) {
        Record *r = kv_A(*block_records, ri);

        if (r->tscripts == NULL) continue;

        // Once per record, not in both phases
        Genotype gt;
        if (phase == 0 && !Genotype_read(r->v, hdr, gt_tag, &gt)) {
            log_warn("Skipping on %s:%li, no usable genotype in FORMAT/%s", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos + 1, gt_tag);
        }

        if (!has_matching_phase(r->v, hdr, gt_tag, phase, include_unphased)) continue;

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
