#include <stdlib.h>
#include <check.h>

#include "check_utils.h"
#include "check_haplotype.h"
#include "check_phase.h"
#include "check_record.h"
#include "check_triggers.h"
#include "check_scores.h"
#include "check_isoform_position.h"
#include "check_region.h"
#include "check_tags.h"
#include "check_splice_prediction.h"
#include "check_isoforms.h"
#include "check_record_block.h"
#include "check_reference_scores.h"

int main(void) {
    Suite *ms = suite_create("Unit test");
    SRunner *sr = srunner_create(ms);

    srunner_add_suite(sr, utils_suite());
    srunner_add_suite(sr, haplotype_suite());
    srunner_add_suite(sr, record_suite());
    srunner_add_suite(sr, phase_suite());
    srunner_add_suite(sr, triggers_suite());
    srunner_add_suite(sr, scores_suite());
    srunner_add_suite(sr, isoform_position_suite());
    srunner_add_suite(sr, region_suite());
    srunner_add_suite(sr, tags_suite());
    srunner_add_suite(sr, splice_prediction_suite());
    srunner_add_suite(sr, isoforms_suite());
    srunner_add_suite(sr, record_block_suite());
    srunner_add_suite(sr, reference_scores_suite());

    srunner_run_all(sr, CK_NORMAL);
    int number_failed = srunner_ntests_failed(sr);

    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
