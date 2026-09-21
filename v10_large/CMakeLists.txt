# Large-data V10 builder.  This is intentionally a separate executable and
# implementation from v10/pre.cpp: its on-disk stages are designed to be
# resumable and schedulable across cluster workers.
add_executable(hic_v10_large
    main.cpp
    io.cpp
    manifest.cpp
    stage.cpp
    sort.cpp
    build.cpp
    writer.cpp
    vectors.cpp
    normalize.cpp
    scale.cpp
)
target_include_directories(hic_v10_large PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${V10_ZSTD_INCLUDE_DIR}
)
target_link_libraries(hic_v10_large PRIVATE
    ${V10_ZSTD_LIBRARY}
    ZLIB::ZLIB
    Threads::Threads
    m
)
set_target_properties(hic_v10_large PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
install(TARGETS hic_v10_large DESTINATION bin)

if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME v10_large_stage_and_rollup
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_large.py
              $<TARGET_FILE:hic_v10_large> $<TARGET_FILE:hic_v10>)
    set_tests_properties(v10_large_stage_and_rollup PROPERTIES
      ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1")
  endif()
endif()
