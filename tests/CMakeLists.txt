enable_testing()
include(GoogleTest)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

file(GLOB_RECURSE TEST_SOURCES test_*.cc)

foreach (TEST_FILE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_FILE_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_FILE_NAME} gtest gtest_main ${PROJECT_NAME})
    set_target_properties(${TEST_FILE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    gtest_discover_tests(${TEST_FILE_NAME})
endforeach ()
