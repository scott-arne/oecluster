cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oecluster VERSION 0.1.1 LANGUAGES CXX)

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
    GIT_TAG v1.0.0
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

# OEShape and its transitive dependencies are not yet in cmake-openeye
# FindOpenEye.cmake; discover them manually.
set(_OE_LIB_SEARCH_PATHS
    ${OPENEYE_LIB_DIR}
    ${OPENEYE_ROOT}/lib
    $ENV{OPENEYE_ROOT}/lib
    $ENV{OE_DIR}/lib
    /opt/openeye/lib
    /usr/local/openeye/lib
)

macro(_find_oe_dep VAR_NAME LIB_NAME TARGET_NAME DEPS)
    if(NOT TARGET ${TARGET_NAME})
        find_library(${VAR_NAME} NAMES ${LIB_NAME} PATHS ${_OE_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
        if(${VAR_NAME})
            add_library(${TARGET_NAME} UNKNOWN IMPORTED)
            set_target_properties(${TARGET_NAME} PROPERTIES
                IMPORTED_LOCATION "${${VAR_NAME}}"
                INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            )
            if(NOT "${DEPS}" STREQUAL "")
                set_property(TARGET ${TARGET_NAME} PROPERTY INTERFACE_LINK_LIBRARIES ${DEPS})
            endif()
        endif()
    endif()
endmacro()

_find_oe_dep(OEOPT_LIBRARY oeopt OpenEye::OEOpt "OpenEye::OESystem")
_find_oe_dep(OEMOLPOTENTIAL_LIBRARY oemolpotential OpenEye::OEMolPotential "OpenEye::OEChem;OpenEye::OEOpt")
_find_oe_dep(OEHERMITE_LIBRARY oehermite OpenEye::OEHermite "OpenEye::OEChem;OpenEye::OEOpt")

if(NOT TARGET OpenEye::OEShape)
    find_library(OESHAPE_LIBRARY NAMES oeshape PATHS ${_OE_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(OESHAPE_LIBRARY)
        add_library(OpenEye::OEShape UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEShape PROPERTIES
            IMPORTED_LOCATION "${OESHAPE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem;OpenEye::OEGrid;OpenEye::OEOpt;OpenEye::OEMolPotential;OpenEye::OEHermite"
        )
        message(STATUS "OpenEye: Found OEShape: ${OESHAPE_LIBRARY}")
    endif()
endif()

if(NOT TARGET OpenEye::OEBio)
    find_library(OEBIO_LIBRARY NAMES oebio PATHS ${_OE_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(OEBIO_LIBRARY)
        add_library(OpenEye::OEBio UNKNOWN IMPORTED)
        set_target_properties(OpenEye::OEBio PROPERTIES
            IMPORTED_LOCATION "${OEBIO_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENEYE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "OpenEye::OEChem"
        )
        message(STATUS "OpenEye: Found OEBio: ${OEBIO_LIBRARY}")
    endif()
endif()

_find_oe_dep(OEZAP_LIBRARY oezap OpenEye::OEZap "OpenEye::OEChem;OpenEye::OEGrid")
_find_oe_dep(OESPICOLI_LIBRARY oespicoli OpenEye::OESpicoli "OpenEye::OEChem;OpenEye::OEZap")
_find_oe_dep(OESITEHOPPER_LIBRARY oesitehopper OpenEye::OESiteHopper "OpenEye::OEChem;OpenEye::OEShape;OpenEye::OESpicoli;${CMAKE_DL_LIBS}")
_find_oe_dep(OESPRUCE_LIBRARY oespruce OpenEye::OESpruce "OpenEye::OEChem;OpenEye::OEBio;OpenEye::OESiteHopper")

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
