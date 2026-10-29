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
target_link_libraries(Haplotype Utils)
target_link_libraries(rescueranger Haplotype)

# Add (misc) utils
add_library(Utils STATIC utils.c)
target_link_libraries(rescueranger Utils)
target_link_libraries(Record Utils) # bcf1_variant_in_region

# Add
add_library(Scores STATIC scores.c)
target_link_libraries(Scores Utils)
target_link_libraries(rescueranger Scores)

# Add BCFTools
add_library(bcftools STATIC bcftools/regidx.c bcftools/gff.c)
target_include_directories(rescueranger PRIVATE bcftools)
target_link_libraries(rescueranger bcftools)

# Add the precomputed reference scores of cpliceai. Optional at run time, and a stub when the cpliceai it is built against has none
add_library(ReferenceScores STATIC reference_scores.c)
target_link_libraries(ReferenceScores Logging ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger ReferenceScores)

# Add session
add_library(Session STATIC session.c)
target_link_libraries(Session Record ReferenceScores bcftools ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger Session)

# Add annotation
add_library(Annotation STATIC annotation.c)
target_link_libraries(Annotation Record Region Logging bcftools ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger Annotation)

# Add phase
add_library(Phase STATIC phase.c)
target_link_libraries(Phase Record Utils Logging ${HTSLIB_LIBRARIES})
target_link_libraries(rescueranger Phase)

# Add triggers
add_library(Triggers STATIC triggers.c)
target_link_libraries(Triggers Record Phase Utils ${HTSLIB_LIBRARIES})
target_link_libraries(rescueranger Triggers)

# Add isoform positions
add_library(IsoformPosition STATIC isoform_position.c)
target_link_libraries(IsoformPosition Record Utils Logging)
target_link_libraries(rescueranger IsoformPosition)

# Add the splice prediction stages (initialising the splice sites and predicting the splice changes of a phase)
add_library(SplicePrediction STATIC splice_prediction.c)
target_link_libraries(SplicePrediction Record Scores Phase Haplotype Utils Logging ReferenceScores bcftools ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger SplicePrediction)

# Add the isoform stage (building the isoforms of a phase)
add_library(Isoforms STATIC isoforms.c)
target_link_libraries(Isoforms Record IsoformPosition Haplotype Triggers Utils Logging ${HTSLIB_LIBRARIES})
target_link_libraries(rescueranger Isoforms)

# Add the output tags
add_library(Tags STATIC tags.c)
target_link_libraries(Tags Record Triggers Utils Logging ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger Tags)

# Add record block
add_library(RecordBlock STATIC record_block.c)
target_link_libraries(RecordBlock Annotation Record Region Utils Logging bcftools ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger RecordBlock)

