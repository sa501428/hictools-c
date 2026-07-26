cmake_minimum_required(VERSION 3.13)
project(HiCTools VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address")

# Required libraries
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# Optional: zstd (for future V10 support)
find_library(ZSTD_LIB zstd)
if(ZSTD_LIB)
  message(STATUS "zstd found: ${ZSTD_LIB}")
  add_compile_definitions(HAVE_ZSTD)
else()
  message(STATUS "zstd not found - skipping V10 support")
endif()

# Common sources shared between both executables
set(COMMON_SOURCES
    common/genome.cpp
    common/pair_parser.cpp
    common/block_writer.cpp
    common/expected_value.cpp
    common/hic_writer.cpp
)

set(COMMON_HEADERS
    common/hic_file_def.h
    common/genome.h
    common/pair_parser.h
    common/block_writer.h
    common/expected_value.h
    common/hic_writer.h
    common/thread_pool.h
    common/little_endian.h
)

# ---- hic_pre executable ----
add_executable(hic_pre
    ${COMMON_SOURCES}
    hic_pre/preprocessor.cpp
    hic_pre/main.cpp
)
target_include_directories(hic_pre PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hic_pre PRIVATE ZLIB::ZLIB Threads::Threads)
if(ZSTD_LIB)
  target_link_libraries(hic_pre PRIVATE ${ZSTD_LIB})
endif()

# ---- hic_addnorm executable ----
add_executable(hic_addnorm
    ${COMMON_SOURCES}
    hic_addnorm/norm_updater.cpp
    hic_addnorm/scale_norm.cpp
    hic_addnorm/main.cpp
)
target_include_directories(hic_addnorm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hic_addnorm PRIVATE ZLIB::ZLIB Threads::Threads)
if(ZSTD_LIB)
  target_link_libraries(hic_addnorm PRIVATE ${ZSTD_LIB})
endif()

# Math library (needed for log2 in SCALE norm)
target_link_libraries(hic_pre PRIVATE m)
target_link_libraries(hic_addnorm PRIVATE m)

# Install targets
install(TARGETS hic_pre hic_addnorm DESTINATION bin)

include(CTest)
if(BUILD_TESTING AND UNIX)
  add_test(
    NAME hic_pre_contiguous_pair_streaming
    COMMAND bash
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_contiguous_pair_streaming.sh
            $<TARGET_FILE:hic_pre>
  )
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
