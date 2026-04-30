cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oecluster VERSION 3.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(OECLUSTER_BUILD_TESTS "Build tests" ON)
option(OECLUSTER_BUILD_PYTHON "Build Python bindings" ON)
option(OECLUSTER_UNIVERSAL2 "Build universal2 binary for macOS" OFF)
option(OECLUSTER_USE_STABLE_ABI "Use Python stable ABI" ON)

if(APPLE)
    add_link_options("LINKER:-no_warn_duplicate_libraries")
endif()

if(APPLE AND OECLUSTER_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
elseif(APPLE)
    unset(CMAKE_OSX_ARCHITECTURES CACHE)
endif()

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
endif()

if(APPLE)
    execute_process(
        COMMAND brew --prefix
        OUTPUT_VARIABLE HOMEBREW_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(HOMEBREW_PREFIX)
        link_directories(${HOMEBREW_PREFIX}/lib)
        # ICU is keg-only on macOS and not linked into the main homebrew lib dir
        execute_process(
            COMMAND brew --prefix icu4c
            OUTPUT_VARIABLE ICU_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(ICU_PREFIX)
            link_directories(${ICU_PREFIX}/lib)
        endif()
    endif()
endif()

include(FetchContent)

FetchContent_Declare(
    cmake_openeye
    GIT_REPOSITORY https://github.com/scott-arne/cmake-openeye.git
    GIT_TAG v1.2.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cmake_openeye)
list(APPEND CMAKE_MODULE_PATH "${cmake_openeye_SOURCE_DIR}")

FetchContent_Declare(
    oeselect
    GIT_REPOSITORY https://github.com/scott-arne/oeselect.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
set(OESELECT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(OESELECT_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(oeselect)

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
    src/CDist.cpp
    src/DistanceMatrix.cpp
    src/metrics/FingerprintMetric.cpp
    src/metrics/ROCSMetric.cpp
    src/metrics/SuperposeMetric.cpp
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
        OpenEye::OEGraphSim
        OpenEye::OEShape
        OpenEye::OEBio
        OpenEye::OESpruce
    PRIVATE
        oeselect
)

if(OECLUSTER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(OECLUSTER_BUILD_PYTHON AND SWIG_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/swig/CMakeLists.txt")
    add_subdirectory(swig)
endif()

option(OECLUSTER_BUILD_TOOLS "Build CLI tools" ON)

if(OECLUSTER_BUILD_TOOLS)
    add_subdirectory(tools)
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
