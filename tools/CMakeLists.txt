include(FetchContent)

FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cli11)

FetchContent_Declare(
    indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(indicators)

FetchContent_Declare(
    oemaestro
    GIT_REPOSITORY https://github.com/scott-arne/oemaestro.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
set(OEMAESTRO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(OEMAESTRO_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(MAEPARSER_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(oemaestro)

if(TARGET maeparser)
    # Boost.Spirit in maeparser triggers libc++'s char_traits<unsigned int>
    # deprecation on recent Clang toolchains; keep that third-party warning
    # out of oecluster builds without changing project-wide warning behavior.
    target_compile_options(maeparser PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Wno-deprecated-declarations>
    )
endif()

add_executable(oepdist oepdist.cpp OutputWriter.cpp MolReader.cpp)

target_link_libraries(oepdist PRIVATE oecluster oemaestro CLI11::CLI11 indicators::indicators)

install(TARGETS oepdist RUNTIME DESTINATION bin)

if(OECLUSTER_BUILD_PYTHON)
    install(TARGETS oepdist RUNTIME DESTINATION oecluster/_bin)
    install(FILES ${CMAKE_SOURCE_DIR}/python/_oepdist_main.py DESTINATION .)
endif()

if(OECLUSTER_BUILD_TESTS)
    add_test(NAME oepdist_fp_metadata
             COMMAND ${CMAKE_COMMAND}
                     -DOEPDIST=$<TARGET_FILE:oepdist>
                     -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/oepdist_fp_metadata
                     -P ${CMAKE_SOURCE_DIR}/tests/cli/oepdist_fp_metadata.cmake)
    add_test(NAME oepdist_fp_behavior
             COMMAND ${CMAKE_COMMAND}
                     -DOEPDIST=$<TARGET_FILE:oepdist>
                     -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/oepdist_fp_behavior
                     -P ${CMAKE_SOURCE_DIR}/tests/cli/oepdist_fp_behavior.cmake)
    add_test(NAME oepdist_fp_explicitness
             COMMAND ${CMAKE_COMMAND}
                     -DOEPDIST=$<TARGET_FILE:oepdist>
                     -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/oepdist_fp_explicitness
                     -P ${CMAKE_SOURCE_DIR}/tests/cli/oepdist_fp_explicitness.cmake)
endif()
