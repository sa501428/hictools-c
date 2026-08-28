# Explicit opt-in; existing V9 executables never link this implementation.
find_path(V10_ZSTD_INCLUDE_DIR zstd.h HINTS /opt/homebrew/include)
find_library(V10_ZSTD_LIBRARY NAMES zstd HINTS /opt/homebrew/lib)
if(NOT V10_ZSTD_INCLUDE_DIR OR NOT V10_ZSTD_LIBRARY)
  message(FATAL_ERROR "BUILD_V10 requires the zstd headers and library")
endif()
add_executable(hic_v10
    main.cpp writer.cpp convert.cpp pre.cpp
    ../common/genome.cpp ../common/pair_parser.cpp ../common/expected_value.cpp
)
target_include_directories(hic_v10 PRIVATE ${PROJECT_SOURCE_DIR} ${V10_ZSTD_INCLUDE_DIR})
target_link_libraries(hic_v10 PRIVATE ${V10_ZSTD_LIBRARY} ZLIB::ZLIB Threads::Threads)
set_target_properties(hic_v10 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
install(TARGETS hic_v10 DESTINATION bin)

if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  add_executable(test_v10_block_geometry tests/test_geometry.cpp writer.cpp)
  target_include_directories(test_v10_block_geometry PRIVATE
      ${PROJECT_SOURCE_DIR} ${V10_ZSTD_INCLUDE_DIR})
  target_link_libraries(test_v10_block_geometry PRIVATE ${V10_ZSTD_LIBRARY})
  add_test(NAME v10_block_geometry COMMAND test_v10_block_geometry)

  set(STRAW_TEST_EXECUTABLE "" CACHE FILEPATH "Optional straw executable for V9/V10 integration tests")
  if(Python3_Interpreter_FOUND)
    add_test(NAME v10_writer_roundtrip
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_writer.py
              $<TARGET_FILE:hic_v10> $<TARGET_FILE:hic_pre> $<TARGET_FILE:hic_addnorm>
              ${STRAW_TEST_EXECUTABLE})
  endif()
endif()
