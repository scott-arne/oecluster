include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(oecluster_tests
    cpp/test_main.cpp
    cpp/test_storage.cpp
    cpp/test_thread_pool.cpp
    cpp/test_pdist.cpp
    cpp/test_fingerprint_metric.cpp
)

target_link_libraries(oecluster_tests
    PRIVATE
        oecluster
        GTest::gtest
)

if(OpenEye_GraphSim_FOUND)
    target_link_libraries(oecluster_tests PRIVATE OpenEye::OEGraphSim)
endif()

include(GoogleTest)
gtest_discover_tests(oecluster_tests)
