include(OpenEyeSWIG)

set(_OECLUSTER_EXPECTED_LIB_VARS
    OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
    OEMATH_LIBRARY OEZSTD_LIBRARY
    OEGRAPHSIM_LIBRARY
    OESHAPE_LIBRARY OEHERMITE_LIBRARY OEOPT_LIBRARY OEMOLPOTENTIAL_LIBRARY
    OEBIO_LIBRARY
    OEGRID_LIBRARY OEZAP_LIBRARY OESPICOLI_LIBRARY OESITEHOPPER_LIBRARY
    OESPRUCE_LIBRARY
)

set(_OECLUSTER_SWIG_DEPENDS
    ${CMAKE_SOURCE_DIR}/include/oecluster/oecluster.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/Error.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/GateFacts.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/PairwiseComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/StorageBackend.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/ThreadPool.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/PDist.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/CDist.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/DistanceMatrix.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/DescriptorStatistics.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/FingerprintComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/ROCSComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/SuperposeComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/DescriptorComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/RMSDComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/MCSComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/ClusterTypes.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Butina.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Representative.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/ClusterReport.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/PartitionAgreement.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/SARCoherence.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/DiversitySelection.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/SetDiversity.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/SphereExclusion.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/KNNGraph.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/JarvisPatrick.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Leiden.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/DBSCAN.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/HDBSCAN.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Agglomerative.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/BitBirch.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/KMedoids.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/MurckoScaffold.h
)

set_source_files_properties(oecluster.i PROPERTIES
    DEPENDS "${_OECLUSTER_SWIG_DEPENDS}"
)

openeye_add_swig_module(
    NAME oecluster
    SWIG_FILE oecluster.i
    LINK_LIBS oecluster
    PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oecluster
    STABLE_ABI ${OECLUSTER_USE_STABLE_ABI}
    SWIG_FLAGS "-doxygen"
        "-I${CMAKE_SOURCE_DIR}/include"
        "-I${OECLUSTER_OEFP_RESOLVED_SOURCE_DIR}/include"
        "-I${OPENEYE_INCLUDE_DIR}"
    EXPECTED_LIB_VARS ${_OECLUSTER_EXPECTED_LIB_VARS}
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oecluster/__init__.py
)

# -----------------------------------------------------------------------------
# Editable dev RPATH fix for pyarrow
# -----------------------------------------------------------------------------
# The extension links libarrow/libparquet out of the active interpreter's
# pyarrow wheel, but those live in site-packages rather than anywhere on the
# default loader search path. openeye_add_swig_module() sets an RPATH covering
# @loader_path and the OpenEye runtime only, so the copy it drops in
# PYTHON_OUTPUT_DIR cannot resolve Arrow and the editable install fails to
# import with "Library not loaded: @rpath/libarrow.<soversion>.dylib".
#
# Patching the copy here keeps the fix tied to the build that produced it:
# pyarrow bumps its soversion on every major release, so a hand-applied
# install_name_tool goes stale silently the next time the environment moves.
#
# Wheels are excluded deliberately. delocate/auditwheel vendor the Arrow
# libraries next to the extension, where the existing @loader_path/$ORIGIN
# entry already finds them, and baking a build machine's absolute
# site-packages path into a redistributable artifact would be wrong.
if(NOT SKBUILD)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c
            "import os, pyarrow; print(os.path.dirname(pyarrow.__file__))"
        OUTPUT_VARIABLE _OECLUSTER_DEV_PYARROW_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )

    if(_OECLUSTER_DEV_PYARROW_DIR AND IS_DIRECTORY "${_OECLUSTER_DEV_PYARROW_DIR}")
        set(_OECLUSTER_DEV_SO
            ${CMAKE_SOURCE_DIR}/python/oecluster/$<TARGET_FILE_NAME:oecluster_python>)

        if(APPLE)
            find_program(INSTALL_NAME_TOOL install_name_tool)
            if(INSTALL_NAME_TOOL)
                # Adding an RPATH invalidates the signature, so re-sign ad hoc.
                add_custom_command(TARGET oecluster_python POST_BUILD
                    COMMAND ${INSTALL_NAME_TOOL} -add_rpath
                        "${_OECLUSTER_DEV_PYARROW_DIR}" ${_OECLUSTER_DEV_SO}
                    COMMAND codesign -f -s - ${_OECLUSTER_DEV_SO}
                    COMMENT "Adding pyarrow RPATH for editable dev install (oecluster)"
                )
            endif()
        elseif(UNIX)
            find_program(PATCHELF patchelf)
            if(PATCHELF)
                add_custom_command(TARGET oecluster_python POST_BUILD
                    COMMAND ${PATCHELF} --add-rpath
                        "${_OECLUSTER_DEV_PYARROW_DIR}" ${_OECLUSTER_DEV_SO}
                    COMMENT "Adding pyarrow RPATH for editable dev install (oecluster)"
                )
            endif()
        endif()
    else()
        message(WARNING
            "pyarrow not importable from ${Python3_EXECUTABLE}; the editable "
            "install will not be able to resolve libarrow at import time.")
    endif()
endif()
