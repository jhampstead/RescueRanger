#include "region.h"

IRegion variant_to_iregion(bcf1_t *v) {
    return (IRegion) {v->rid, v->pos, v->pos + v->rlen};
}

IRegion iregion_from_bases(int rid, int64_t first, int64_t last) {
    return (IRegion) { rid, first, last + 1 };
}

int64_t iregion_last_base(IRegion r) {
    return r.end > r.start ? r.end - 1 : r.start;
}

int iregion_intersects(IRegion r1, IRegion r2) {
    return r1.start <= r2.end-1 && r2.start <= r1.end-1 && r1.rid == r2.rid;
}
