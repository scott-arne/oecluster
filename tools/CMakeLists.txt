include(FetchContent)

FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cli11)

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

add_executable(oepdist oepdist.cpp OutputWriter.cpp MolReader.cpp)

target_link_libraries(oepdist PRIVATE oecluster oemaestro CLI11::CLI11)

install(TARGETS oepdist RUNTIME DESTINATION bin)

if(OECLUSTER_BUILD_PYTHON)
    install(TARGETS oepdist RUNTIME DESTINATION oecluster/_bin)
endif()
