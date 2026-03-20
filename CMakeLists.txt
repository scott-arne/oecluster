cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oecluster VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(OECLUSTER_BUILD_TESTS "Build tests" ON)
option(OECLUSTER_BUILD_PYTHON "Build Python bindings" ON)
option(OECLUSTER_UNIVERSAL2 "Build universal2 binary for macOS" OFF)

if(APPLE AND OECLUSTER_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version" FORCE)
endif()

include(FetchContent)

FetchContent_Declare(
    cmake_openeye
    GIT_REPOSITORY https://github.com/scott-arne/cmake-openeye.git
    GIT_TAG v1.0.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cmake_openeye)
list(APPEND CMAKE_MODULE_PATH "${cmake_openeye_SOURCE_DIR}")

if(OECLUSTER_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(SWIG 4.0)
endif()

if(OPENEYE_USE_SHARED)
    find_package(OpenEyePython REQUIRED)
endif()

find_package(OpenEye REQUIRED)

set(OECLUSTER_SOURCES
    src/StorageBackend.cpp
    src/ThreadPool.cpp
    src/PDist.cpp
    src/DistanceMatrix.cpp
)

add_library(oecluster ${OECLUSTER_SOURCES})

set_target_properties(oecluster PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(oecluster
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(oecluster
    PUBLIC
        OpenEye::OEChem
)

if(OECLUSTER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(OECLUSTER_BUILD_PYTHON AND SWIG_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/swig/CMakeLists.txt")
    add_subdirectory(swig)
endif()

install(TARGETS oecluster
    EXPORT oecluster-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/oecluster DESTINATION include)
install(EXPORT oecluster-targets
    FILE OEClusterTargets.cmake
    NAMESPACE OECluster::
    DESTINATION lib/cmake/oecluster
)
