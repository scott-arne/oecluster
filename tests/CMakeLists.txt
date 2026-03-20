include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(oecluster_tests
    cpp/test_main.cpp
)

target_link_libraries(oecluster_tests
    PRIVATE
        oecluster
        GTest::gtest
)

include(GoogleTest)
gtest_discover_tests(oecluster_tests)
