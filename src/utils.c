#include "utils.h"
#include "region.h"
#include <string.h>

int get_allele_offset(int64_t pos, int rlen, int alen, int boundary_start, int boundary_end) {
    if (pos - boundary_start < 0) { // Insertion overlapping start
        return 1 - (rlen + (pos - boundary_start));
    } else if (pos + rlen - boundary_end > 0) { // Insertion overlapping end
        return 1 - (rlen - (pos + rlen - boundary_end));
    }
    return alen - rlen; // WARN: Only supports a single alternate allele per variant
}

int get_allele_offset_variant(bcf1_t *v, int boundary_start, int boundary_end) {
    return get_allele_offset(v->pos, v->rlen, strlen(v->d.allele[1]), boundary_start, boundary_end);
}

int64_t get_position_in_haplotype(int64_t pos, Records *rs, int64_t start, int64_t stop) {
    int offset = 0;
    for (int i = 0; i < rs->n; i++) {
        Record *r = kv_A(*rs, i);
        if (r->v->pos >= pos) break;

        offset += get_allele_offset_variant(r->v, start, stop);
    }

    return pos + offset;
}

int get_haplotype_sequence(kstring_t ref, Records *rs, int seq_offset, int pos_offset, kstring_t *out_hap) {
    // Prepare kstring
    kstring_t seq;
    seq.m = ref.l + seq_offset + 1; 
    seq.l = 0;
    seq.s = malloc(seq.m);

    // Current state tracking
    int prev_cpy_end = 0;
    int c_offset = 0;

    IRegion del_region = {-1, -1, -1};
    int prev_alt_allele_len = 0;
    for (int i = 0; i < rs->n; i++) {
        // Get current variant info
        Record *r = kv_A(*rs, i);

        int var_idx = r->v->pos - pos_offset;
        if (var_idx < 0) { // If a deletion overlaps the transcript start
            prev_cpy_end -= var_idx;
            continue;
        }

        if (var_idx + r->v->rlen > ref.l) { // If a deletion overlaps the transcript end
            ref.l -= (var_idx + r->v->rlen) - ref.l;
            continue;
        }

        char *dst = seq.s + prev_cpy_end + c_offset;
        char *src = ref.s + prev_cpy_end;
        int segment_len = var_idx - prev_cpy_end;
        kputsn(src, segment_len, &seq);

        // Copy over variant
        int alt_allele_len = strlen(r->v->d.allele[1]);
        kputsn(r->v->d.allele[1], alt_allele_len, &seq);

        // Update state for next iteration
        prev_cpy_end = var_idx + r->v->rlen;
        c_offset += alt_allele_len - r->v->rlen;

        if (alt_allele_len == 1 && r->v->rlen > 1) del_region = variant_to_iregion(r->v);
    }

    // Copy over remaining reference sequence
    kputsn(ref.s + prev_cpy_end, ref.l - prev_cpy_end, &seq);

    *out_hap = seq;

    return 0;
}

