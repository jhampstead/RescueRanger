#include "region.h"


Region get_region(bcf_hdr_t *hdr, bcf1_t *r) {
    return (Region) { bcf_hdr_id2name(hdr, r->rid), r->pos, r->pos + r->rlen }; // NOTE: Should this be +1?
}

void format_region(Region *region, kstring_t *s) {
    kputs(region->chr, s);
    kputc(':', s);
    kputl(region->start, s);
    kputc('-', s);
    kputl(region->end, s);
}

Region iregion_to_region(IRegion r, bcf_hdr_t *hdr) {
    return (Region) { bcf_hdr_id2name(hdr, r.rid), r.start, r.end };
}

IRegion variant_to_iregion(bcf1_t *v) {
    return (IRegion) {v->rid, v->pos, v->pos + v->rlen};
}

int iregion_intersects(IRegion r1, IRegion r2) {
    return r1.start <= r2.end-1 && r2.start <= r1.end-1 && r1.rid == r2.rid;
}
