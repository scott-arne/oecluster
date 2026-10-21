# A homebrew SWIG upgrade deletes the versioned Cellar directory that
# CMakeCache.txt recorded as SWIG_DIR (and that the generated build rules use
# as SWIG_LIB). FindSWIG trusts a cached SWIG_DIR or SWIG_EXECUTABLE without
# checking that it still exists, and recomputes SWIG_VERSION only when it is
# unset, so a vanished path otherwise breaks the build until the cache is
# edited by hand. When either cached path has vanished, drop all three so
# find_package(SWIG) detects SWIG afresh. A cached path that still exists is
# left alone.
set(_oecluster_swig_stale FALSE)
foreach(_oecluster_swig_var SWIG_EXECUTABLE SWIG_DIR)
    if(DEFINED CACHE{${_oecluster_swig_var}})
        set(_oecluster_swig_path "$CACHE{${_oecluster_swig_var}}")
        if(_oecluster_swig_path AND NOT _oecluster_swig_path MATCHES "-NOTFOUND$"
                AND NOT EXISTS "${_oecluster_swig_path}")
            message(STATUS "Dropping cached ${_oecluster_swig_var}, which names "
                           "a path that no longer exists: ${_oecluster_swig_path}")
            set(_oecluster_swig_stale TRUE)
        endif()
    endif()
endforeach()
if(_oecluster_swig_stale)
    unset(SWIG_EXECUTABLE CACHE)
    unset(SWIG_DIR CACHE)
    unset(SWIG_VERSION CACHE)
endif()
unset(_oecluster_swig_stale)
unset(_oecluster_swig_var)
unset(_oecluster_swig_path)
