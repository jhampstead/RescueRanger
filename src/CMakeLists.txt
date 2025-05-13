set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_C_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3")
set(CMAKE_C_FLAGS_RELEASE_INIT "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

#find_package(Sanitizers)

# Create executable
add_executable(rescueranger rescueranger.c)
#add_sanitizers(rescueranger)

# Link external libraries
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

# Add record
add_library(Record STATIC record.c)
target_link_libraries(rescueranger Record)

# Add haplotype
add_library(Haplotype STATIC haplotype.c)
target_link_libraries(rescueranger Haplotype)

# Add (misc) utils
add_library(Utils STATIC utils.c)
target_link_libraries(rescueranger Utils)

# Add
add_library(Scores STATIC scores.c)
target_link_libraries(rescueranger Scores)

# Add BCFTools
add_library(bcftools STATIC bcftools/regidx.c bcftools/gff.c)
target_include_directories(rescueranger PRIVATE bcftools)
target_link_libraries(rescueranger bcftools)

