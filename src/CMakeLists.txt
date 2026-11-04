set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_C_FLAGS_DEBUG_INIT "-Wall -g -DLOG_USE_COLOR")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3")
set(CMAKE_C_FLAGS_RELEASE_INIT "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

#find_package(Sanitizers)

# Create executable
add_executable(rescueranger rescueranger.c)
#add_sanitizers(rescueranger)

# cpliceai and htslib are not given to the program itself, but come in through the libraries below, which give them after
# bcftools. CMake puts what is given to the program before what the libraries bring in, so given here they would come before
# objects that are only taken from the libraries later. When pkg-config finds cpliceai its libraries include htslib (-lhts), and
# htslib has a region index with the same function names as bcftools/regidx.c, which would then be the one used.

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
add_library(bcftools STATIC bcftools/gff.c bcftools/regidx.c bcftools/error.c)
target_link_libraries(bcftools ${HTSLIB_LIBRARIES})
set_target_properties(bcftools PROPERTIES C_VISIBILITY_PRESET hidden) # Not exported by the programs, so libraries that use htslib's region index keep it
target_include_directories(rescueranger PRIVATE bcftools)
target_link_libraries(rescueranger bcftools)

# Add the precomputed reference scores of cpliceai. Optional at run time, and a stub when the cpliceai it is built against has none
add_library(ReferenceScores STATIC reference_scores.c)
target_link_libraries(ReferenceScores Logging ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger ReferenceScores)

# Add session
add_library(Session STATIC session.c)
target_link_libraries(Session Record ReferenceScores Tags Phase Logging bcftools ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
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

# Add the public interface (include/rescueranger/rescueranger.h): the options, opening a session and running the stages
add_library(Api STATIC api.c)
target_link_libraries(Api Session RecordBlock Phase Triggers SplicePrediction Isoforms Tags Utils Logging ${HTSLIB_LIBRARIES} ${CPLICEAI_LIBRARIES})
target_link_libraries(rescueranger Api)

# The shared library that is installed for other programs, such as the bcftools plugin. It is built from the same sources as the
# libraries above, which are not installed. Only the functions of include/rescueranger/rescueranger.h are exported (RR_EXPORT).
add_library(rescueranger_shared SHARED
    api.c
    annotation.c
    haplotype.c
    isoform_position.c
    isoforms.c
    phase.c
    record.c
    record_block.c
    reference_scores.c
    region.c
    scores.c
    session.c
    splice_prediction.c
    tags.c
    triggers.c
    utils.c
    logging/log.c
    bcftools/error.c
    bcftools/gff.c
    bcftools/regidx.c
)
set_target_properties(rescueranger_shared PROPERTIES
    OUTPUT_NAME rescueranger
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    C_VISIBILITY_PRESET hidden)
target_compile_definitions(rescueranger_shared PRIVATE -DLOG_USE_COLOR)
target_link_libraries(rescueranger_shared ${CPLICEAI_LIBRARIES} ${HTSLIB_LIBRARIES} m)

install(TARGETS rescueranger rescueranger_shared
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
