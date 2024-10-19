set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_C_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# Create executable
add_executable(rescueranger rescueranger.c)

# Link external libraries
# target_link_libraries(rescueranger asan)
target_link_libraries(rescueranger ${HTSLIB_LIBRARIES})
target_link_libraries(rescueranger ${CPLICEAI_LIBRARIES})

# Add logging
add_library(Logging STATIC logging/log.c)
target_compile_definitions(Logging PRIVATE -DLOG_USE_COLOR)
target_include_directories(rescueranger PRIVATE Logging)
target_link_libraries(rescueranger Logging)

# Add region
add_library(Region STATIC region.c)
target_link_libraries(rescueranger Region)

# Add range
#add_library(Range STATIC range.c)
#target_link_libraries(rescueranger Range)

# Add record
add_library(Record STATIC record.c)
target_link_libraries(rescueranger Record)

# Add utils
add_library(Utils STATIC utils.c)
target_link_libraries(rescueranger Utils)

# Add BCFTools
add_library(bcftools STATIC bcftools/regidx.c bcftools/gff.c)
target_include_directories(rescueranger PRIVATE bcftools)
target_link_libraries(rescueranger bcftools)

