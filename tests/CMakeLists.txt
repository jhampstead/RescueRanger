add_subdirectory(unity)

# Runs all tests
add_executable(region_test region_test.c)
add_executable(utils_test utils_test.c)

target_link_libraries(region_test PRIVATE Unity Region)
target_include_directories(region_test PRIVATE ${CMAKE_SOURCE_DIR}/src unity)

target_link_libraries(utils_test PRIVATE Unity Utils Region Logging ${HTSLIB_LIBRARIES})
target_include_directories(utils_test PRIVATE ${CMAKE_SOURCE_DIR}/src unity)

# Add tests
add_test(region region_test)
add_test(utils utils_test)

