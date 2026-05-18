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
    ${CMAKE_SOURCE_DIR}/include/oecluster/PairwiseComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/StorageBackend.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/ThreadPool.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/PDist.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/CDist.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/DistanceMatrix.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/FingerprintComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/ROCSComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/comparisons/SuperposeComparison.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/ClusterTypes.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Butina.h
    ${CMAKE_SOURCE_DIR}/include/oecluster/clustering/Centroid.h
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
        "-I${OPENEYE_INCLUDE_DIR}"
    EXPECTED_LIB_VARS ${_OECLUSTER_EXPECTED_LIB_VARS}
    INIT_PY ${CMAKE_SOURCE_DIR}/python/oecluster/__init__.py
)
