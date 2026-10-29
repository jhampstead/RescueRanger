#include "isoforms.h"

#include <stdlib.h>

#include <klib/kvec.h>
#include <klib/ksort.h>

#include "logging/log.h"
#include "scores.h"
#include "triggers.h"
#include "utils.h"

// seq_pos is the index into the coding sequence in reading direction, so the shifts of both strands are sorted ascending
#define frameshift_lt(a, b) ((a).seq_pos < (b).seq_pos)
KSORT_INIT(frameshift_sort, IsoformVariant, frameshift_lt)

void calculate_rescue_probabilities(Isoforms *isoforms_struct) {
    isoforms_struct->rescue_probability = 0.0f;
    RescueProbabilities *rescue_probs = malloc(sizeof(RescueProbabilities));

    rescue_probs->n_isoforms = malloc(isoforms_struct->m_rescue_types * sizeof(int32_t));
    rescue_probs->rescues = malloc(isoforms_struct->m_rescue_types * sizeof(uint8_t));
    rescue_probs->rescue_types = malloc(isoforms_struct->m_rescue_types * sizeof(enum RescueType));
    rescue_probs->rescue_probabilities = malloc(isoforms_struct->m_rescue_types * sizeof(float));
    rescue_probs->affected_bases = malloc(isoforms_struct->m_rescue_types * sizeof(int32_t));

    for (int i = 0; i < isoforms_struct->m_rescue_types; i++) {
        rescue_probs->n_isoforms[i] = 0;
        rescue_probs->rescues[i] = 0;
        rescue_probs->rescue_probabilities[i] = 0.0f;
        rescue_probs->rescue_types[i] = NO_RESCUE; // initialise as NO_RESCUE
        rescue_probs->affected_bases[i] = 0;
    }

    for (size_t i = 0; i < kv_size(isoforms_struct->isoforms); i++) {
        const Isoform *isoform = &kv_A(isoforms_struct->isoforms, i);

        // If rescue is predicted on the isoform, add isoform probability to the rescue probability [INFO/RESCUE_PROB]
        if (isoform->rescue == 1) {
            isoforms_struct->rescue_probability += isoform->posterior_probability;
        }

        // Get probabilities by rescue type [INFO/RESCUE_TYPE]. If the rescue type is already recorded, add the isoform and sum the probability
        int found = 0;
        for (int j = 0; j < isoforms_struct->n_rescue_types; j++) {
            if ((rescue_probs->rescue_types[j] != isoform->rescue_type) || (rescue_probs->affected_bases[j] != isoform->affected_bases)) continue;

            rescue_probs->n_isoforms[j] += 1;
            rescue_probs->rescue_probabilities[j] += isoform->posterior_probability;
            if (isoform->affected_bases > rescue_probs->affected_bases[j]) rescue_probs->affected_bases[j] = isoform->affected_bases;
            found = 1;
            break;
        }

        // If rescue type isn't found, add a new rescue type, dynamically reallocating memory if necessary
        if (found) continue;

        if (isoforms_struct->n_rescue_types >= isoforms_struct->m_rescue_types) {
            isoforms_struct->m_rescue_types *= 2; // double the size
            rescue_probs->n_isoforms = realloc(rescue_probs->n_isoforms, isoforms_struct->m_rescue_types * sizeof(int32_t));
            rescue_probs->rescues = realloc(rescue_probs->rescues, isoforms_struct->m_rescue_types * sizeof(uint8_t));
            rescue_probs->rescue_probabilities = realloc(rescue_probs->rescue_probabilities, isoforms_struct->m_rescue_types * sizeof(float));
            rescue_probs->rescue_types = realloc(rescue_probs->rescue_types, isoforms_struct->m_rescue_types * sizeof(enum RescueType));
            rescue_probs->affected_bases = realloc(rescue_probs->affected_bases, isoforms_struct->m_rescue_types * sizeof(int32_t));

            for (int k = isoforms_struct->n_rescue_types; k < isoforms_struct->m_rescue_types; k++) {
                rescue_probs->n_isoforms[k] = 0;
                rescue_probs->rescues[k] = 0;
                rescue_probs->rescue_probabilities[k] = 0.0f;
                rescue_probs->rescue_types[k] = NO_RESCUE; // initialise with NULL
                rescue_probs->affected_bases[k] = 0;
            }
        }

        rescue_probs->n_isoforms[isoforms_struct->n_rescue_types] = 1;
        rescue_probs->rescues[isoforms_struct->n_rescue_types] = isoform->rescue;
        rescue_probs->rescue_types[isoforms_struct->n_rescue_types] = isoform->rescue_type;
        rescue_probs->rescue_probabilities[isoforms_struct->n_rescue_types] = isoform->posterior_probability;
        rescue_probs->affected_bases[isoforms_struct->n_rescue_types] = isoform->affected_bases;
        isoforms_struct->n_rescue_types++;
    }

    isoforms_struct->rescue = (isoforms_struct->rescue_probability > 0.5) ? 1 : 0; // set rescue value on the basis of rescue_probability
    isoforms_struct->rescue_type_probabilities = rescue_probs; // assign probabilities per rescue type back to struct

}

void calculate_ptv_probabilities(Isoforms *isoforms_struct) {
    float ptv_probability = 0.0;
    float stop_gained_probability = 0.0;
    float frameshift_variant_probability = 0.0;
    float long_inframe_indel_probability = 0.0;

    for (size_t i = 0; i < kv_size(isoforms_struct->isoforms); i++) {
        const Isoform *isoform = &kv_A(isoforms_struct->isoforms, i);
        if (isoform->ptv == 0) continue;

        ptv_probability += isoform->posterior_probability;

        if (isoform->ptv_type == STOP_GAINED) {
            stop_gained_probability += isoform->posterior_probability;
        } else if (isoform->ptv_type == FRAMESHIFT) {
            frameshift_variant_probability += isoform->posterior_probability;
        } else if (isoform->ptv_type == LONG_INFRAME_INDEL) {
            long_inframe_indel_probability += isoform->posterior_probability;
        }
    }

    isoforms_struct->ptv_probability = ptv_probability;
    isoforms_struct->ptv = (isoforms_struct->ptv_probability > 0.5) ? 1 : 0;

    if (stop_gained_probability > frameshift_variant_probability) {
        isoforms_struct->ptv_type = STOP_GAINED;
    } else if (frameshift_variant_probability > stop_gained_probability) {
        isoforms_struct->ptv_type = FRAMESHIFT;
    } else if (long_inframe_indel_probability > stop_gained_probability && long_inframe_indel_probability > stop_gained_probability) {
        isoforms_struct->ptv_type = LONG_INFRAME_INDEL;
    } else {
        isoforms_struct->ptv_type = NO_PTV; // tie or no PTV
    }
}

// Splice sites further than this from the canonical site of the intron are not considered as alternatives
#define MAX_ALTERNATIVE_SITE_DISTANCE 500

typedef kvec_t(SpliceSite *) SpliceSitePtrs;

// Appends the sites that are still used (posterior above 0). With limit_distance, only the sites within
// MAX_ALTERNATIVE_SITE_DISTANCE of the canonical site are appended.
static void collect_sites(SpliceSite *sites, const size_t n, const SpliceSite *canonical, const bool limit_distance, SpliceSitePtrs *out) {
    for (size_t i = 0; i < n; i++) {
        SpliceSite *site = sites + i;
        if (site->posterior_probability == 0) continue;
        if (limit_distance && labs((int64_t) site->position - (int64_t) canonical->position) > MAX_ALTERNATIVE_SITE_DISTANCE) continue;

        kv_push(SpliceSite *, *out, site);
    }
}

// Appends the used canonical sites and then the used cryptic sites of an exon, donors or acceptors. Cryptic sites are always
// limited to the distance from the canonical site, canonical sites only with limit_canonical.
static void collect_exon_sites(const SpliceSites *sites, const bool donors, const SpliceSite *canonical, const bool limit_canonical, SpliceSitePtrs *out) {
    if (donors) {
        collect_sites(sites->canonical_donors, sites->n_canonical_donors, canonical, limit_canonical, out);
        collect_sites(sites->cryptic_donors, sites->n_cryptic_donors, canonical, true, out);
    } else {
        collect_sites(sites->canonical_acceptors, sites->n_canonical_acceptors, canonical, limit_canonical, out);
        collect_sites(sites->cryptic_acceptors, sites->n_cryptic_acceptors, canonical, true, out);
    }
}

// A UTR, which is always spliced out. The record is the one the isoforms are built for.
static NonCodingRegion utr_region(const enum UTRType type, const uint32_t strand, const Range gene, const Range boundary_cds, Record *r) {
    const Range utr = Utr_range(type, strand, gene, boundary_cds);

    NonCodingRegion ncr = { 1, utr, utr, NULL, NULL, type, 0 };
    kv_init(ncr.records);
    kv_push(Record *, ncr.records, r);

    return ncr;
}

// An intron that is retained. Its range runs between the canonical sites of the exons the donor and the acceptor belong to.
static NonCodingRegion retained_region(const uint32_t strand, SpliceSite *donor, SpliceSite *acceptor) {
    const Range intron = Intron_between_sites(strand, donor->canonical_position, acceptor->canonical_position);

    return (NonCodingRegion) { 0, intron, intron, donor, acceptor, 0, 0 };
}

// Appends every way the intron after exon ei (in transcription order, step is +1 forward and -1 reverse) can be spliced on
// the haplotype. If both canonical sites are lost the intron is retained. Otherwise the used donors of the exon and the used
// acceptors of the next exon are paired, and when a canonical site is ablated, the nearby sites of the coding exon beyond it
// are added. If one side has no used site, the intron is retained once for every used site of the other side.
// The regions point at the splice sites of the exons, so they are only valid while the exons have them.
static void add_intron_regions(const Transcript *t, const int ei, const int step, const int first_cds, const int last_cds, NonCodingRegions *group) {
    const int partner = ei + step;
    const int alternative_donor_exon = ei - step;
    const int alternative_acceptor_exon = partner + step;

    SpliceSite *canonical_donor = kv_A(t->exons, ei)->splice_sites->canonical_donors;
    const bool lost_donor = canonical_donor->predicted && canonical_donor->posterior_probability == 0;
    const bool ablated_donor = canonical_donor->predicted && (canonical_donor->ref_prior_probability - canonical_donor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;

    SpliceSite *canonical_acceptor = kv_A(t->exons, partner)->splice_sites->canonical_acceptors;
    const bool lost_acceptor = canonical_acceptor->predicted && canonical_acceptor->posterior_probability == 0;
    const bool ablated_acceptor = canonical_acceptor->predicted && (canonical_acceptor->ref_prior_probability - canonical_acceptor->prior_probability) >= SPLICEAI_RESCUE_THRESHOLD;

    if (lost_donor && lost_acceptor) {
        kv_push(NonCodingRegion, *group, retained_region(t->strand, canonical_donor, canonical_acceptor));
        return;
    }

    SpliceSitePtrs donors; kv_init(donors);
    collect_exon_sites(kv_A(t->exons, ei)->splice_sites, true, canonical_donor, false, &donors);
    if (ablated_donor && alternative_donor_exon >= first_cds && alternative_donor_exon <= last_cds) {
        collect_exon_sites(kv_A(t->exons, alternative_donor_exon)->splice_sites, true, canonical_donor, true, &donors);
    }

    SpliceSitePtrs acceptors; kv_init(acceptors);
    collect_exon_sites(kv_A(t->exons, partner)->splice_sites, false, canonical_acceptor, false, &acceptors);
    if (ablated_acceptor && alternative_acceptor_exon >= first_cds && alternative_acceptor_exon <= last_cds) {
        collect_exon_sites(kv_A(t->exons, alternative_acceptor_exon)->splice_sites, false, canonical_acceptor, true, &acceptors);
    }

    if (kv_size(donors) == 0 || kv_size(acceptors) == 0) {
        for (size_t i = 0; i < kv_size(donors); i++) {
            kv_push(NonCodingRegion, *group, retained_region(t->strand, kv_A(donors, i), canonical_acceptor));
        }
        for (size_t i = 0; i < kv_size(acceptors); i++) {
            kv_push(NonCodingRegion, *group, retained_region(t->strand, canonical_donor, kv_A(acceptors, i)));
        }
    } else {
        const Range canonical_intron = Intron_between_sites(t->strand, canonical_donor->canonical_position, canonical_acceptor->canonical_position);

        for (size_t i = 0; i < kv_size(donors); i++) {
            SpliceSite *donor = kv_A(donors, i);
            for (size_t j = 0; j < kv_size(acceptors); j++) {
                SpliceSite *acceptor = kv_A(acceptors, j);

                const Range intron = Intron_between_sites(t->strand, donor->position, acceptor->position);
                if (intron.start > intron.end) continue; // A donor beyond the acceptor cannot form an intron

                const int uses_competitor = canonical_donor->canonical_position != donor->canonical_position || canonical_acceptor->canonical_position != acceptor->canonical_position;
                kv_push(NonCodingRegion, *group, ((NonCodingRegion) { 1, intron, canonical_intron, donor, acceptor, 0, uses_competitor }));
            }
        }
    }

    kv_destroy(donors);
    kv_destroy(acceptors);
}

bool build_noncoding_regions(const Transcript *t, const Gene *gene, Record *r, NonCodingRegionSets *groups) {
    kv_init(*groups);

    int first_cds = -1, last_cds = -1;
    for (int ei = 0; ei < kv_size(t->exons); ei++) {
        if (kv_A(t->exons, ei)->cds == NULL) continue;

        if (first_cds < 0) first_cds = ei;
        last_cds = ei;
    }
    if (first_cds < 0) return false;
    if (t->strand != STRAND_FWD && t->strand != STRAND_REV) return true;

    // The first and last coding exon in transcription order
    const int step = t->strand == STRAND_FWD ? 1 : -1;
    const int transcript_first = step > 0 ? first_cds : last_cds;
    const int transcript_last = step > 0 ? last_cds : first_cds;

    for (int ei = transcript_first; ; ei += step) {
        Exon *e = kv_A(t->exons, ei);

        if (ei == transcript_first) {
            NonCodingRegions utr5; kv_init(utr5);
            kv_push(NonCodingRegion, utr5, utr_region(UTR5, t->strand, gene->coords, *e->cds, r));
            kv_push(NonCodingRegions, *groups, utr5);
        }

        if (ei == transcript_last) {
            NonCodingRegions utr3; kv_init(utr3);
            kv_push(NonCodingRegion, utr3, utr_region(UTR3, t->strand, gene->coords, *e->cds, r));
            kv_push(NonCodingRegions, *groups, utr3);
            break;
        }

        NonCodingRegions introns; kv_init(introns);
        add_intron_regions(t, ei, step, first_cds, last_cds, &introns);
        kv_push(NonCodingRegions, *groups, introns);
    }

    return true;
}

void free_noncoding_region_groups(NonCodingRegionSets *groups) {
    for (size_t g = 0; g < kv_size(*groups); g++) {
        NonCodingRegions *group = &kv_A(*groups, g);
        for (size_t i = 0; i < kv_size(*group); i++) kv_destroy(kv_A(*group, i).records);
    }

    NonCodingRegionSets_destroy(groups);
}

// The index of a gene base in the coding sequence of the isoform, in reading direction, with the coding part of the indels of
// the haplotype applied. Returns whether the base is coding.
static int isoform_position_of(const int64_t position, const NonCodingRegions *ncrs, const Records *hap_in, const Transcript *t, const Gene *gene, int64_t *isoform_position) {
    *isoform_position = isoform_coding_index(hap_in, ncrs, gene->coords, t->strand, position);

    return isoform_base_is_coding(ncrs, position);
}

// Appends one splice shift for every record that changed the splice site
static void add_splice_shifts(const SpliceSite *site, const int coding, const int partial, const int64_t isoform_position, const int shift, const int modified_bases, IsoformVariants *out) {
    const Records *records = (const Records *) site->records_ptr;

    for (size_t j = 0; j < kv_size(*records); j++) {
        IsoformVariant variant = { coding, partial, SPLICE, isoform_position, shift, modified_bases, kv_A(*records, j) };
        kv_push(IsoformVariant, *out, variant);
    }
}

void collect_splice_shifts(const NonCodingRegions *ncrs, const Records *hap_in, const Transcript *t, const Gene *gene, IsoformVariants *out) {
    const int strand_correction = t->strand == STRAND_FWD ? 1 : -1;

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);
        int64_t isoform_position;

        const bool donor_lost = ncr.donor != NULL && ncr.donor->posterior_probability == 0;
        const bool acceptor_lost = ncr.acceptor != NULL && ncr.acceptor->posterior_probability == 0;
        if (donor_lost || acceptor_lost) {
            const int64_t position = Intron_donor_base(ncr.canonical_coords, t->strand);
            const int offset = -range_len(ncr.canonical_coords); // INTRON_RETENTION

            const int coding = isoform_position_of(position, ncrs, hap_in, t, gene, &isoform_position);
            if (donor_lost) add_splice_shifts(ncr.donor, coding, 0, isoform_position, offset, offset, out);
            if (acceptor_lost) add_splice_shifts(ncr.acceptor, coding, 0, isoform_position, offset, offset, out);
            continue;
        }

        if (ncr.donor != NULL) {
            const int modified_bases = strand_correction * (ncr.donor->position - ncr.donor->canonical_position);
            const int shift = strand_correction * (Intron_donor_base(ncr.coords, t->strand) - Intron_donor_base(ncr.canonical_coords, t->strand));
            const int64_t position = modified_bases < 0 ? Intron_donor_base(ncr.canonical_coords, t->strand) : Intron_donor_base(ncr.coords, t->strand);

            const int coding = isoform_position_of(position, ncrs, hap_in, t, gene, &isoform_position);
            add_splice_shifts(ncr.donor, coding, !ncr.uses_competitor, isoform_position, shift, modified_bases, out);
        }

        if (ncr.acceptor != NULL) {
            const int modified_bases = strand_correction * (ncr.acceptor->position - ncr.acceptor->canonical_position);
            const int shift = strand_correction * (Intron_acceptor_base(ncr.coords, t->strand) - Intron_acceptor_base(ncr.canonical_coords, t->strand));
            const int64_t position = Intron_acceptor_base(ncr.coords, t->strand);

            const int coding = isoform_position_of(position, ncrs, hap_in, t, gene, &isoform_position);
            add_splice_shifts(ncr.acceptor, coding, !ncr.uses_competitor, isoform_position, shift, modified_bases, out);
        }
    }
}

void collect_isoform_shifts(const IsoformContext *ctx, const NonCodingRegions *ncrs, IsoformVariants *shifts) {
    const Transcript *t = ctx->t;

    for (size_t i = 0; i < kv_size(*ctx->hap_in); i++) {
        Record *hap_r = kv_A(*ctx->hap_in, i);

        int64_t isoform_position;
        int indel_offset = 0;
        const int coding = !get_range_in_isoform(&hap_r->alleles, ncrs, ctx->hap_in, t->strand, ctx->gene->coords, &isoform_position, &indel_offset);

        const int offset = coding ? indel_offset : 0; // Coding bases deleted minus bases inserted into the coding sequence

        IsoformVariant isoshift = { coding, 1, INDEL, isoform_position, offset, offset, hap_r };
        kv_push(IsoformVariant, *shifts, isoshift);
    }

    for (size_t i = 0; i < kv_size(*ncrs); i++) {
        const NonCodingRegion ncr = kv_A(*ncrs, i);

        for (size_t j = 0; j < kv_size(ncr.records); j++) {
            Record *utr_r = kv_A(ncr.records, j);
            enum UTRType utr_type = is_utr_variant(&utr_r->alleles, ncrs);

            if (utr_type == UTR5) set_record_flag(utr_r, ctx->ti, ctx->phase, UTR5_VARIANT);
            if (utr_type == UTR3) set_record_flag(utr_r, ctx->ti, ctx->phase, UTR3_VARIANT);
        }
    }

    collect_splice_shifts(ncrs, ctx->hap_in, t, ctx->gene, shifts);

    if (kv_size(*shifts) > 0) ks_introsort_frameshift_sort(kv_size(*shifts), shifts->a);

    size_t kept = 0;
    for (size_t i = 0; i < kv_size(*shifts); i++) {
        const IsoformVariant shift = kv_A(*shifts, i);
        if (shift.r == NULL) continue;

        if ((shift.shift_type == SPLICE && shift.shift != 0) || (shift.shift_type == INDEL && shift.coding)) {
            kv_A(*shifts, kept++) = shift;
        }
    }
    shifts->n = kept;
}

enum RescueType shift_rescue_type(const IsoformVariant *shift) {
    if (shift->partial) return shift->shift < 0 ? PARTIAL_INTRON_RETENTION : PARTIAL_EXON_SKIPPING;
    return shift->shift < 0 ? INTRON_RETENTION : EXON_SKIPPING;
}

// For every shift of the record that leaves a PTV on the isoform, flags the other records with a shift from there on in the
// sorted order, and a trigger on the same transcript, with UPSTREAM_PTV
static void flag_upstream_ptvs(const IsoformContext *ctx, const IsoformVariants *shifts, const int rescued, const enum PTVType ptv_type) {
    const TscriptCsqs *csq = ctx->csq;
    const int stop_trigger = Trigger_value(csq, STOP_TRIGGER, ctx->phase);

    for (size_t i = 0; i < kv_size(*shifts); i++) {
        const IsoformVariant iso_shift = kv_A(*shifts, i);
        if (iso_shift.r != ctx->r) continue;

        if (
            (csq->triggers[SPLICE_TRIGGER] == AF_SPLICE_VALUE || stop_trigger == STOP_GAIN_VALUE) &&
            (csq->triggers[SPLICE_TRIGGER] != CANONICAL_SPLICE_VALUE && !csq->triggers[FRAMESHIFT_TRIGGER] && stop_trigger != STOP_RESCUE_VALUE) &&
            (ptv_type == NO_PTV || ptv_type == LONG_INFRAME_INDEL)
        ) continue;
        if (rescued && !(ptv_type == STOP_GAINED || ptv_type == FRAMESHIFT)) continue;

        for (size_t j = i; j < kv_size(*shifts); j++) {
            const IsoformVariant upstream_shift = kv_A(*shifts, j);
            if (upstream_shift.r == iso_shift.r) continue;

            for (int ti2 = 0; ti2 < upstream_shift.r->tscripts->n; ti2++) {
                Transcript *tscript = kv_A(*upstream_shift.r->tscripts, ti2);
                if (tscript == NULL || tscript->id != ctx->t->id) continue;

                TscriptCsqs *upstream_csq = kv_A(*upstream_shift.r->vep, ti2);
                if (!(upstream_csq->triggers[SPLICE_TRIGGER] || upstream_csq->triggers[FRAMESHIFT_TRIGGER] || Trigger_value(upstream_csq, STOP_TRIGGER, ctx->phase))) continue;

                set_record_flag(upstream_shift.r, ti2, ctx->phase, UPSTREAM_PTV);
            }
        }
    }
}

// The earlier of two positions in the isoform sequence, where -1 is no position
static int64_t earlier_stop(const int64_t a, const int64_t b) {
    if (a < 0) return b;
    if (b < 0) return a;

    return a < b ? a : b;
}

// The escape region (NMD_START_PROXIMAL, NMD_LONG_EXON or NMD_DISTAL, 0 for none) of the first base of the new stop codon (-1 for
// none): the escape regions of the reference transcript are put in the coordinates of the isoform, which hap_ncrs (the non-coding
// regions in haplotype coordinates) and the records of the haplotype tell
static uint32_t nmd_flag_of_stop(const IsoformContext *ctx, const NonCodingRegions *hap_ncrs, const int64_t stop) {
    const Transcript *t = ctx->t;

    Ranges segments = isoform_segments(hap_ncrs, ctx->gene->coords.start, t->strand);
    NmdRegions regions = Transcript_nmd_regions(t);
    NmdRegions escape;
    nmd_escape_ranges_in_isoform(&regions, ctx->hap_in, ctx->gene->coords, &segments, t->strand, &escape);

    const uint32_t flag = nmd_flag_for_stop(stop, &escape);

    NmdRegions_destroy(&escape);
    NmdRegions_destroy(&regions);
    kv_destroy(segments);

    return flag;
}

uint32_t nmd_flag_of_isoforms(const Isoforms *isoforms) {
    const Isoform *best = NULL;
    for (size_t k = 0; k < kv_size(isoforms->isoforms); k++) {
        const Isoform *isoform = &kv_A(isoforms->isoforms, k);
        if (!isoform->nmd_ptv) continue;
        if (best == NULL || isoform->posterior_probability > best->posterior_probability) best = isoform;
    }

    return best == NULL ? 0 : best->nmd_flag;
}

Isoform evaluate_isoform(const IsoformContext *ctx, const NonCodingRegions *ncrs) {
    Record *r = ctx->r;
    const Transcript *t = ctx->t;
    const TscriptCsqs *csq = ctx->csq;
    const int ti = ctx->ti;
    const int phase = ctx->phase;
    const int stop_trigger = Trigger_value(csq, STOP_TRIGGER, phase); // Includes the reading frame conflicts of this phase

    log_debug("Gene length is: %i", ctx->hap_seq->l);

    NonCodingRegions hap_ncrs; kv_init(hap_ncrs);
    for (size_t k = 0; k < kv_size(*ncrs); k++) {
        NonCodingRegion ncr = kv_A(*ncrs, k);
        ncr.coords = region_to_haplotype(ctx->hap_in, ctx->gene->coords, ncr.coords, ncr.is_utr != 0);
        ncr.canonical_coords = region_to_haplotype(ctx->hap_in, ctx->gene->coords, ncr.canonical_coords, ncr.is_utr != 0);
        kv_push(NonCodingRegion, hap_ncrs, ncr);
    }

    // The isoform sequence, 5' to 3'
    kstring_t spliced_seq = {0};
    splice_haplotype_sequence(ctx->hap_seq->s, &hap_ncrs, ctx->gene->coords.start, t->strand, &spliced_seq);
    if (t->strand == STRAND_REV) reverse_complement(spliced_seq.s, spliced_seq.l);

    IsoformVariants shifts; kv_init(shifts);
    collect_isoform_shifts(ctx, ncrs, &shifts);

    log_info("%s:%li (isoform_length=%i)", bcf_hdr_id2name(ctx->hdr, r->v->rid), r->v->pos, spliced_seq.l);

    int affected_bases = 0;
    int rescued = 0;
    int unrescued_stop = 0, frameshift_checked = 0;
    enum RescueType rescue_type = NO_RESCUE;
    enum PTVType ptv_type = NO_PTV;
    int touched = 0;

    // Where the new stop of an unrescued or induced PTV is, for the NMD flag: indices into the isoform sequence, -1 for none yet
    int64_t vep_stop = -1;                 // The stop of a VEP stop that the isoform keeps
    int64_t gained_stop = -1;              // The earliest stop of the induced stops that set STOP_GAINED
    int64_t splice_frameshift_from = -1;   // Where the splice shift that set FRAMESHIFT is
    int64_t vep_frameshift_from = -1;      // Where the frameshift of the record itself is, which only the frameshift trigger looks at

    for (size_t ri = 0; ri < kv_size(shifts); ri++) {
        const IsoformVariant iso_shift = kv_A(shifts, ri);

        if (iso_shift.r != r) continue;
        if (unrescued_stop) continue;

        // Check if constitutive and in frame
        for (size_t ei = 0; ei < kv_size(t->exons); ei++) {
            const Exon exon = *kv_A(t->exons, ei);
            if (exon.cds == NULL) continue;
            if (!variant_in_cds_exon(&iso_shift.r->alleles, &exon)) continue;

            if (!exon.constitutive) set_record_flag(r, ti, phase, NON_CONSTITUTIVE);
            if (exon.in_frame) set_record_flag(r, ti, phase, INFRAME);
            break;
        }

        touched = 1;
        if (stop_trigger) {
            // Check if the codons of the variant code for a stop
            // The position is the index into the coding sequence. A deletion of the first or last coding bases places it just
            // outside the sequence
            int64_t coding_index = iso_shift.seq_pos;
            if (coding_index >= (int64_t) spliced_seq.l) coding_index = (int64_t) spliced_seq.l - 1;
            if (coding_index < 0) coding_index = 0;
            const Alleles *a = &iso_shift.r->alleles;
            const int scanned = (a->prefix_len > 0) + a->alt_len; // The attached base of an indel and the ALT bases

            // From the start of the codon of the first base up to the last base of the change. The window is measured from the
            // first base: a change that begins late in a codon and runs into the next one has the next codon looked at too
            const Range window = codon_window(coding_index, scanned);
            const bool is_stop = has_stop_codon(spliced_seq.s, spliced_seq.l, window.start, range_len(window));
            const int64_t found_stop = first_stop_codon(spliced_seq.s, spliced_seq.l, window.start, range_len(window), true);

            if (!is_stop && stop_trigger == STOP_RESCUE_VALUE) {
                if (rescue_type == NO_RESCUE) {
                    rescue_type = COMPLETE_RESCUE;
                    rescued = 1;
                }
            } else if (stop_trigger == STOP_RESCUE_VALUE) {
                unrescued_stop = 1;
                vep_stop = earlier_stop(vep_stop, found_stop);
            }

            if (is_stop && (stop_trigger == STOP_GAIN_VALUE)) {
                ptv_type = STOP_GAINED;
                gained_stop = earlier_stop(gained_stop, found_stop);
            }
        }

        int should_check_induced = 0;
        Range possible_induced_sequence = {0};

        if (csq->triggers[SPLICE_TRIGGER] && !unrescued_stop && iso_shift.shift_type == SPLICE) {
            if (iso_shift.shift == 0) {
                rescue_type = COMPLETE_RESCUE;
                rescued = 1;
                continue;
            }

            rescue_type = shift_rescue_type(&iso_shift);
            affected_bases = abs(iso_shift.shift);

            const int splice_insertion_length = iso_shift.modified_bases < 0 ? abs(iso_shift.modified_bases) : 0;
            if (iso_shift.shift % 3 == 0 && abs(iso_shift.shift) < MAX_SPLICE_RADIUS) {
                rescued = 1;
                should_check_induced = 1;
                possible_induced_sequence = codon_window_within(iso_shift.seq_pos, splice_insertion_length + 1, spliced_seq.l);
            } else if (csq->triggers[SPLICE_TRIGGER] && iso_shift.shift % 3 != 0) {
                rescued = 0;
                ptv_type = FRAMESHIFT;
                splice_frameshift_from = iso_shift.seq_pos;
            } else if (csq->triggers[SPLICE_TRIGGER] && abs(iso_shift.shift) > MAX_SPLICE_RADIUS) {
                ptv_type = LONG_INFRAME_INDEL;
                should_check_induced = 1;
                possible_induced_sequence = codon_window_within(iso_shift.seq_pos, splice_insertion_length + 1, spliced_seq.l);
            }
        }

        if ((csq->triggers[FRAMESHIFT_TRIGGER] || csq->triggers[SPLICE_TRIGGER]) && (ptv_type == FRAMESHIFT || ptv_type == NO_PTV || rescued == 1) && !unrescued_stop) {
            if (iso_shift.shift == 0) continue;
            frameshift_checked = 1;

            if (iso_shift.shift % 3 == 0 && abs(iso_shift.shift) < MAX_FRAMESHIFT_RADIUS) {
                rescued = 1;
                rescue_type = shift_rescue_type(&iso_shift);
                affected_bases = abs(iso_shift.shift);
            } else {
                // The change of the record itself, which nothing may restore: its stop is looked for from here
                if (iso_shift.shift_type != SPLICE && vep_frameshift_from < 0) vep_frameshift_from = iso_shift.seq_pos;

                int32_t summed_shift = 0;
                int64_t shift_start = 0;
                int included_current_shift = 0;
                int induced_frameshift_rescue = 0;

                IsoformVariants unique_shifts;
                unique_isoform_shifts(&shifts, &unique_shifts);

                for (size_t i = 0; i < kv_size(unique_shifts); i++) {
                    const IsoformVariant shift = kv_A(unique_shifts, i);
                    if (shift.shift == 0) continue; // Rescued splice, do nothing
                    if (abs(shift.shift) > MAX_FRAMESHIFT_RADIUS) continue;

                    // If outside of rescue radius
                    if (llabs(shift_start - shift.seq_pos) > MAX_FRAMESHIFT_RADIUS) {
                        summed_shift = shift.shift;
                        shift_start = shift.seq_pos;
                        included_current_shift = 0 | (shift.r == r);
                        induced_frameshift_rescue = shift.shift_type == SPLICE;
                        continue;
                    }

                    // Keep track of whether the set of frameshifts includes the target
                    if (!included_current_shift && shift.r == r) included_current_shift = 1;

                    summed_shift += shift.shift;

                    induced_frameshift_rescue = shift.shift_type == SPLICE && shift.shift % 3 != 0;
                    if (summed_shift % 3 != 0) continue;

                    if (included_current_shift && shift.seq_pos + MAX_FRAMESHIFT_RADIUS > iso_shift.seq_pos && shift.seq_pos - MAX_FRAMESHIFT_RADIUS < iso_shift.seq_pos) {
                        rescued = 1;
                        rescue_type = (rescue_type == NO_RESCUE) ? FRAME_RESTORING_INDEL : rescue_type; // Only set to rescue type
                        int64_t shift_end = shift.seq_pos;
                        affected_bases = llabs(shift_end - shift_start) + (shift.modified_bases < 0 ? abs(shift.modified_bases) : 0);

                        if (induced_frameshift_rescue) set_record_flag(r, ti, phase, RESCUED_BY_SPLICE_FRAMESHIFT);

                        if (ptv_type == FRAMESHIFT) ptv_type = NO_PTV;

                        should_check_induced = 1;
                        possible_induced_sequence = codon_window_within(shift_start, affected_bases, spliced_seq.l);

                        affected_bases = range_len(possible_induced_sequence);
                    }

                    summed_shift = 0;
                    shift_start = 0;
                    included_current_shift = 0;
                }

                kv_destroy(unique_shifts);
            }
        }

        if (!should_check_induced) continue;

        // Check if the induced sequence codes for a stop
        if (has_stop_codon(spliced_seq.s, spliced_seq.l, possible_induced_sequence.start, range_len(possible_induced_sequence))) {
            ptv_type = STOP_GAINED;
            gained_stop = earlier_stop(gained_stop, first_stop_codon(spliced_seq.s, spliced_seq.l, possible_induced_sequence.start,
                                                                      range_len(possible_induced_sequence), true));
        }
    }

    if (!touched && !unrescued_stop && !frameshift_checked && rescue_type == NO_RESCUE) {
        rescued = 1;
        rescue_type = COMPLETE_RESCUE;
    }

    // The NMD flag of a PTV that the isoform keeps (a VEP stop that is not rescued, or a frameshift that nothing restores) or
    // induces (a stop, a frameshift or a long in-frame change): the escape region of the first stop the ribosome meets. A stop
    // that is found where the change is, or by reading on from it for a frameshift. Without one there is no flag
    bool has_ptv = unrescued_stop || ptv_type != NO_PTV;
    int64_t first_stop = -1;
    if (unrescued_stop) first_stop = earlier_stop(first_stop, vep_stop);
    if (ptv_type == STOP_GAINED) first_stop = earlier_stop(first_stop, gained_stop);
    if (ptv_type == FRAMESHIFT) first_stop = earlier_stop(first_stop, ptv_stop_downstream(spliced_seq.s, spliced_seq.l, splice_frameshift_from));
    if (ptv_type == NO_PTV && !rescued && !unrescued_stop && csq->triggers[FRAMESHIFT_TRIGGER] && vep_frameshift_from >= 0) {
        has_ptv = true;
        first_stop = earlier_stop(first_stop, ptv_stop_downstream(spliced_seq.s, spliced_seq.l, vep_frameshift_from));
    }
    const uint32_t nmd_flag = has_ptv ? nmd_flag_of_stop(ctx, &hap_ncrs, first_stop) : 0;

    flag_upstream_ptvs(ctx, &shifts, rescued, ptv_type);

    const Isoform isoform = {
        .prior_probability = isoform_prior(&hap_ncrs),
        .posterior_probability = 0.0,
        .rescue = rescued,
        .rescue_type = rescue_type,
        .affected_bases = affected_bases,
        .ptv = (ptv_type == STOP_GAINED || ptv_type == FRAMESHIFT || ptv_type == LONG_INFRAME_INDEL),
        .ptv_type = ptv_type,
        .nmd_ptv = has_ptv,
        .nmd_flag = nmd_flag
    };

    kv_destroy(shifts);
    kv_destroy(hap_ncrs);
    free(spliced_seq.s);

    return isoform;
}

void build_isoforms_for_phase(const Session *s, const int phase, const Records *records) {
    for (int i = 0; i < records->n; i++) {
        Record *r = kv_A(*records, i);

        if (!variant_fits_window(&r->alleles, s->prediction_radius)) continue;

        for (int ti = 0; ti < r->tscripts->n; ti++) {
            Transcript *t = kv_A(*r->tscripts, ti);
            Gene *gene = kv_A(*r->transcript_genes, ti);
            TscriptCsqs *csq = kv_A(*r->vep, ti);

            if (t == NULL || gene == NULL) continue;
            if (!csq->triggers[SPLICE_TRIGGER] && !csq->triggers[FRAMESHIFT_TRIGGER] && !Trigger_value(csq, STOP_TRIGGER, phase)) continue;
            if (!record_overlaps_gene(r, gene)) continue;

            Records hap_in;
            find_haplotype_records(records, r, gene, false, &hap_in, NULL);

            NonCodingRegionSets nc_region_groups;
            if (!build_noncoding_regions(t, gene, r, &nc_region_groups)) { // Transcript is not protein-coding after all!
                free_noncoding_region_groups(&nc_region_groups);
                kv_destroy(hap_in);
                continue;
            }

            // One isoform per combination of one non-coding region from each group
            NonCodingRegionSets intron_sets;
            enumerate_region_combinations(nc_region_groups.a, kv_size(nc_region_groups), &intron_sets);

            kstring_t hap_seq;
            get_haplotype_sequence(gene->seq, &hap_in, gene->coords, &hap_seq);
            log_info("%s:%li (num_isoforms=%zu)", bcf_hdr_id2name(s->hdr, r->v->rid), r->v->pos, kv_size(intron_sets));

            Isoforms isoforms = { .valid = 1, .m_rescue_types = 1, .ptv_type = NO_PTV };
            kv_init(isoforms.isoforms);

            // TODO: On rework the order of this should be changed, ALL splice effects should be resolved before haplotype effects get applied.
            // Right now the back and forth on the coordinates makes a mess of things.
            const IsoformContext ctx = { s->hdr, phase, r, ti, t, gene, csq, &hap_in, &hap_seq };
            for (size_t k = 0; k < kv_size(intron_sets); k++) {
                log_debug("Processing isoform: %zu", k);
                kv_push(Isoform, isoforms.isoforms, evaluate_isoform(&ctx, &kv_A(intron_sets, k)));
            }

            // Calculate isoforms posterior
            float probability_sum = 0.0;
            for (size_t k = 0; k < kv_size(isoforms.isoforms); k++) probability_sum += kv_A(isoforms.isoforms, k).prior_probability;
            for (size_t k = 0; k < kv_size(isoforms.isoforms); k++) kv_A(isoforms.isoforms, k).posterior_probability = kv_A(isoforms.isoforms, k).prior_probability / probability_sum;

            calculate_rescue_probabilities(&isoforms);
            calculate_ptv_probabilities(&isoforms);
            set_record_flag(r, ti, phase, nmd_flag_of_isoforms(&isoforms));

            r->isoforms[ti * NUM_PHASES + phase] = isoforms;

            NonCodingRegionSets_destroy(&intron_sets);
            free_noncoding_region_groups(&nc_region_groups);
            kv_destroy(hap_in);
            free(hap_seq.s);
        }
    }
}
