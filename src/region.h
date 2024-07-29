#include <stdint.h>

#include <klib/kstring.h>
#include <htslib/vcf.h>

typedef struct Region {
    const char *chr;
    int64_t start;
    int64_t end;
} Region;

typedef struct IRegion {
    int rid;
    int64_t start;
    int64_t end;
} IRegion;


Region get_region(bcf_hdr_t *hdr, bcf1_t *r);

void format_region(Region *region, kstring_t *s);

Region iregion_to_region(IRegion r, bcf_hdr_t *hdr);

IRegion variant_to_iregion(bcf1_t *v);

int iregion_intersects(IRegion r1, IRegion r2);

