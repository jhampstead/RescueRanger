# Add Unity as a static library
add_library(Unity STATIC unity.c)

# Specify include directories for Unity
target_include_directories(Unity PUBLIC ${CMAKE_SOURCE_DIR}/tests/unity)
