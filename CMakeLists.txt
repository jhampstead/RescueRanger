cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_STANDARD 99)

project(RescueRanger C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(include)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
pkg_check_modules(CPLICEAI REQUIRED cpliceai)

add_subdirectory(src)

include(CTest)
add_subdirectory(tests)
