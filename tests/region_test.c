#include "unity/unity.h"
#include "unity/unity_internals.h"

#include <region.h>


void setUp(void) {}
void tearDown(void) {}

void IRegionIntersects_OverlappingRegions_ShouldIntersect(void) {
    IRegion r1 = {0, 0, 7};
    IRegion r2 = {0, 5, 10};

    TEST_ASSERT_TRUE(iregion_intersects(r1, r2));
    TEST_ASSERT_TRUE(iregion_intersects(r2, r1));

    r1 = (IRegion) {0, 0, 5};
    TEST_ASSERT_TRUE(iregion_intersects(r1, r2));
}

void IRegionIntersects_TouchingRegions_ShouldNotIntersect(void) {
    IRegion r1 = {0, 0, 5};
    IRegion r2 = {0, 6, 10};

    TEST_ASSERT_FALSE(iregion_intersects(r1, r2));
}

void IRegionIntersects_DifferentRIDs_ShouldNotIntersect(void) {
    IRegion r1 = {0, 0, 7};
    IRegion r2 = {1, 0, 7};

    TEST_ASSERT_FALSE(iregion_intersects(r1, r2));
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(IRegionIntersects_OverlappingRegions_ShouldIntersect);
    RUN_TEST(IRegionIntersects_TouchingRegions_ShouldNotIntersect);
    RUN_TEST(IRegionIntersects_DifferentRIDs_ShouldNotIntersect);

    return UNITY_END();
}
