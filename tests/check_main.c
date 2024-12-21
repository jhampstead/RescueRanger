#include <stdlib.h>
#include <check.h>

#include "check_utils.h"
#include "check_haplotype.h"

int main(void) {
    Suite *ms = suite_create("Unit test");
    SRunner *sr = srunner_create(ms);

    srunner_add_suite(sr, utils_suite());
    srunner_add_suite(sr, haplotype_suite());

    srunner_run_all(sr, CK_NORMAL);
    int number_failed = srunner_ntests_failed(sr);

    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
