cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 99)

# The version of the library: RR_VERSION in rescueranger/version.h, the version of librescueranger.so and of rescueranger.pc
project(RescueRanger VERSION 1.0.0 LANGUAGES C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/externals/sanitizers-cmake/cmake" ${CMAKE_MODULE_PATH})

include(GNUInstallDirs)

include_directories(include)

# rescueranger/version.h is generated, next to the headers in include/
configure_file(include/rescueranger/version.h.in "${PROJECT_BINARY_DIR}/include/rescueranger/version.h" @ONLY)
include_directories("${PROJECT_BINARY_DIR}/include")

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
pkg_check_modules(CHECK REQUIRED check)

# cpliceai. By default it is the one pkg-config finds, which is how it has always been found. To build against a particular
# version, give the prefix it is installed in (its headers are in <prefix>/include/cpliceai and libcpliceai in <prefix>/lib):
#
#   cmake -S . -B build-old -DCPLICEAI_ROOT=/opt/cpliceai/old
#   cmake -S . -B build-new -DCPLICEAI_ROOT=/opt/cpliceai/new
#
# Use a build directory for each choice. The CMake cache keeps the choice, and the binaries of the build directories do not
# share anything. What the version offers is detected from its headers (see the definitions below), so the same sources build
# against the older cpliceai and against the newer one.
set(CPLICEAI_ROOT "" CACHE PATH "Install prefix of the cpliceai to build against (default: the one pkg-config finds)")

if(CPLICEAI_ROOT)
    find_path(CPLICEAI_INCLUDE_DIR cpliceai/predict.h HINTS "${CPLICEAI_ROOT}/include" NO_DEFAULT_PATH)
    find_library(CPLICEAI_LIBRARY cpliceai HINTS "${CPLICEAI_ROOT}/lib" "${CPLICEAI_ROOT}/lib64" NO_DEFAULT_PATH)
    if(NOT CPLICEAI_INCLUDE_DIR OR NOT CPLICEAI_LIBRARY)
        message(FATAL_ERROR "cpliceai not found in CPLICEAI_ROOT=${CPLICEAI_ROOT}: expected include/cpliceai/predict.h and "
                            "lib/libcpliceai.so there (found: include ${CPLICEAI_INCLUDE_DIR}, library ${CPLICEAI_LIBRARY})")
    endif()

    set(CPLICEAI_HEADER_ROOT "${CPLICEAI_INCLUDE_DIR}")
    set(CPLICEAI_LIBRARIES "${CPLICEAI_LIBRARY}")
    include_directories("${CPLICEAI_INCLUDE_DIR}")

    # Every cpliceai has the soname libcpliceai.so.1, and the runtime loader takes the first it finds. The directory of the one
    # that was linked is put in the binaries as DT_RPATH, which is searched before LD_LIBRARY_PATH (DT_RUNPATH is not), so a
    # cpliceai on LD_LIBRARY_PATH or in /usr/local/lib cannot take its place without anyone noticing.
    #
    # The binaries in the build directory get the same RPATH as the installed ones. Otherwise CMake ends the RPATH of the build
    # directory with an empty entry, to rewrite it at install, and the runtime loader reads an empty entry as the current directory:
    # the libraries that are not in the directories before it, such as libhts.so.3, would be looked for in the directory the
    # program is run from. Nothing in the build directory links a shared library of the build directory, so nothing else is needed.
    get_filename_component(CPLICEAI_LIBRARY_DIR "${CPLICEAI_LIBRARY}" DIRECTORY)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
    set(CMAKE_INSTALL_RPATH "${CPLICEAI_LIBRARY_DIR}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")

    # The newer cpliceai is linked against ONNX Runtime. The linker checks that what libcpliceai.so refers to is defined, and for
    # that it has to find the libraries libcpliceai.so needs: -rpath-link says where. The same directories go into the binaries as
    # RPATH, so they are found when the program runs too. Where is the directory of onnxruntime for pkg-config, or what is given:
    #   -DCPLICEAI_DEPENDENCY_DIRS=/opt/onnxruntime/lib
    set(CPLICEAI_DEPENDENCY_DIRS "" CACHE STRING "Directories with the libraries cpliceai needs, such as ONNX Runtime (default: where pkg-config finds onnxruntime)")
    set(CPLICEAI_DEPENDENCY_DIRS_USED ${CPLICEAI_DEPENDENCY_DIRS})
    if(NOT CPLICEAI_DEPENDENCY_DIRS_USED)
        pkg_check_modules(ONNXRUNTIME QUIET onnxruntime)
        if(ONNXRUNTIME_FOUND)
            set(CPLICEAI_DEPENDENCY_DIRS_USED ${ONNXRUNTIME_LIBRARY_DIRS})
        endif()
    endif()
    foreach(CPLICEAI_DEPENDENCY_DIR IN LISTS CPLICEAI_DEPENDENCY_DIRS_USED)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath-link,${CPLICEAI_DEPENDENCY_DIR}")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath-link,${CPLICEAI_DEPENDENCY_DIR}")
        list(APPEND CMAKE_INSTALL_RPATH "${CPLICEAI_DEPENDENCY_DIR}")
    endforeach()

    message(STATUS "cpliceai: ${CPLICEAI_LIBRARY} (headers in ${CPLICEAI_INCLUDE_DIR}/cpliceai), its libraries in: ${CPLICEAI_DEPENDENCY_DIRS_USED}")
else()
    pkg_check_modules(CPLICEAI REQUIRED cpliceai)

    # The sources include <cpliceai/range.h> and the like, so the directory that holds the cpliceai directory is what matters
    if(CPLICEAI_PREFIX)
        set(CPLICEAI_HEADER_ROOT "${CPLICEAI_PREFIX}/include")
    else()
        set(CPLICEAI_HEADER_ROOT "/usr/local/include")
    endif()

    message(STATUS "cpliceai: found by pkg-config (${CPLICEAI_LIBRARIES}), headers in ${CPLICEAI_HEADER_ROOT}/cpliceai")
endif()

# What this cpliceai offers, from its installed headers
set(CPLICEAI_HAS_PREDICT_PADDED_SEQUENCE OFF)
if(EXISTS "${CPLICEAI_HEADER_ROOT}/cpliceai/predict.h")
    file(READ "${CPLICEAI_HEADER_ROOT}/cpliceai/predict.h" CPLICEAI_PREDICT_HEADER)
    string(FIND "${CPLICEAI_PREDICT_HEADER}" "predict_padded_sequence" CPLICEAI_PREDICT_PADDED_POSITION)
    if(NOT CPLICEAI_PREDICT_PADDED_POSITION EQUAL -1)
        set(CPLICEAI_HAS_PREDICT_PADDED_SEQUENCE ON)
    endif()
endif()

set(CPLICEAI_HAS_REFERENCE_SCORES OFF)
if(EXISTS "${CPLICEAI_HEADER_ROOT}/cpliceai/reference.h" AND EXISTS "${CPLICEAI_HEADER_ROOT}/cpliceai/gene_reference.h" AND EXISTS "${CPLICEAI_HEADER_ROOT}/cpliceai/gene_regions.h")
    set(CPLICEAI_HAS_REFERENCE_SCORES ON)
endif()

# predict_padded_sequence: encode, predict and, on the reverse strand, turn the result around. The older cpliceai has none, and its
# one_hot_encode allocates the encoding itself, where the newer one fills a buffer, so predict_for_seq has a version for each.
if(CPLICEAI_HAS_PREDICT_PADDED_SEQUENCE)
    add_definitions(-DCPLICEAI_HAS_PREDICT_PADDED_SEQUENCE)
endif()

# Precomputed reference scores (reference.h, gene_reference.h, gene_regions.h). Without them the option to use them is refused at
# run time with a message, and the reference predictions are computed from the sequence as they have always been.
if(CPLICEAI_HAS_REFERENCE_SCORES)
    add_definitions(-DCPLICEAI_HAS_REFERENCE_SCORES)
endif()

message(STATUS "cpliceai predict_padded_sequence: ${CPLICEAI_HAS_PREDICT_PADDED_SEQUENCE}, precomputed reference scores: ${CPLICEAI_HAS_REFERENCE_SCORES}")

add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
add_test(NAME unit_test COMMAND unit_test)

# Installed: the shared library and the program (see src/CMakeLists.txt), the headers and rescueranger.pc
install(FILES "${PROJECT_SOURCE_DIR}/include/rescueranger/rescueranger.h" "${PROJECT_BINARY_DIR}/include/rescueranger/version.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/rescueranger")
configure_file(rescueranger.pc.in "${PROJECT_BINARY_DIR}/rescueranger.pc" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/rescueranger.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

