include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

find_package(Check REQUIRED)
include_directories(${CHECK_INCLUDE_DIRS})
link_directories(${CHECK_LIBRARY_DIRS})


set(TEST_SOURCES
    mocks.c
    check_haplotype.c
    check_utils.c
    check_main.c
)

add_executable(unit_test ${TEST_SOURCES})
add_dependencies(unit_test rescueranger)

target_link_libraries(unit_test ${CHECK_LIBRARIES})

target_link_libraries(unit_test ${HTSLIB_LIBRARIES})
target_link_libraries(unit_test ${CPLICEAI_LIBRARIES})


target_link_libraries(unit_test Record)
target_link_libraries(unit_test Utils)
target_link_libraries(unit_test Haplotype)

