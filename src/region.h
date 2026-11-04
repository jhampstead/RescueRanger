#ifndef REGION_HEADER
#define REGION_HEADER

#include <stdint.h>

#include <htslib/vcf.h>

// A region of a sequence by its id, half-open: the bases start to end - 1.
typedef struct IRegion {
    int rid;
    int64_t start;
    int64_t end;
} IRegion;

// The bases of the REF allele of v (including the base an indel is attached to)
IRegion variant_to_iregion(bcf1_t *v);

// The region of the bases first to last, both inclusive, as the GFF and the region index give them
IRegion iregion_from_bases(int rid, int64_t first, int64_t last);

// The last base of the region, which is the end a region index query takes (inclusive). An empty region gives its start.
int64_t iregion_last_base(IRegion r);

// Whether the regions are on the same sequence and share at least one base
int iregion_intersects(IRegion r1, IRegion r2);

#endif
