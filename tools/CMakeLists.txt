include(FetchContent)

FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.4.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cli11)

add_executable(oepdist oepdist.cpp OutputWriter.cpp MolReader.cpp)

target_link_libraries(oepdist PRIVATE oecluster CLI11::CLI11)

if(OpenEye_GraphSim_FOUND)
    target_link_libraries(oepdist PRIVATE OpenEye::OEGraphSim)
endif()

if(OpenEye_Shape_FOUND)
    target_link_libraries(oepdist PRIVATE OpenEye::OEShape)
endif()

if(OpenEye_Bio_FOUND)
    target_link_libraries(oepdist PRIVATE OpenEye::OEBio)
endif()

install(TARGETS oepdist RUNTIME DESTINATION bin)
