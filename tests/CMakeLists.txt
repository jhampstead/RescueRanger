include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

find_package(Check REQUIRED)
include_directories(${CHECK_INCLUDE_DIRS})
link_directories(${CHECK_LIBRARY_DIRS})


set(TEST_SOURCES
    mocks.c
    check_haplotype.c
    check_utils.c
    check_record.c
    check_phase.c
    check_triggers.c
    check_scores.c
    check_isoform_position.c
    check_region.c
    check_tags.c
    check_splice_prediction.c
    check_isoforms.c
    check_record_block.c
    check_reference_scores.c
    check_main.c
)

add_executable(unit_test ${TEST_SOURCES})

target_link_libraries(unit_test ${CHECK_LIBRARIES})

target_link_libraries(unit_test RecordBlock)
target_link_libraries(unit_test Session)
target_link_libraries(unit_test ReferenceScores)
target_link_libraries(unit_test Isoforms)
target_link_libraries(unit_test SplicePrediction)
target_link_libraries(unit_test Tags)
target_link_libraries(unit_test IsoformPosition)
target_link_libraries(unit_test Scores)
target_link_libraries(unit_test Triggers)
target_link_libraries(unit_test Phase)
target_link_libraries(unit_test Record)
target_link_libraries(unit_test Utils)
target_link_libraries(unit_test Haplotype)
target_link_libraries(unit_test Region)

# cpliceai and htslib come in through the libraries above, after bcftools. Given here they would come before those libraries, and
# when pkg-config finds cpliceai its libraries include htslib (-lhts), whose region index has the same function names as
# bcftools/regidx.c: the region index would be the one of htslib.
