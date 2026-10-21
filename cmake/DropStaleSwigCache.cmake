# FindSWIG trusts a cached SWIG_EXECUTABLE or SWIG_DIR without checking that
# the path still exists, and recomputes SWIG_VERSION only when it is unset. A
# Homebrew SWIG upgrade deletes the versioned Cellar directory that
# CMakeCache.txt recorded as SWIG_DIR (and that the generated build rules pass
# to SWIG as SWIG_LIB), so without this guard the build stays broken until the
# cache is edited or deleted by hand.
#
# When the cached SWIG_EXECUTABLE names a path that no longer exists, drop it
# together with SWIG_DIR and SWIG_VERSION: the library directory belongs to
# that executable, and keeping it would pair a newly found executable with an
# old library. Otherwise, when the cached SWIG_DIR names a path that no longer
# exists, drop SWIG_DIR and SWIG_VERSION but keep the executable, so FindSWIG
# asks that executable for its library directory again (swig -swiglib). A live
# SWIG_DIR is dropped only along with a vanished executable, and a live
# executable is never dropped.
#
# The guard runs only when CMake configures. CMake does not treat the SWIG
# installation as a regeneration input, so after a SWIG upgrade a plain
# `cmake --build` does not re-run CMake and keeps the stale SWIG_LIB in the
# generated rules. Re-run CMake on the build directory (`cmake <build-dir>`)
# to apply the guard.
set(_oecluster_swig_executable "$CACHE{SWIG_EXECUTABLE}")
set(_oecluster_swig_dir "$CACHE{SWIG_DIR}")
set(_oecluster_swig_executable_gone FALSE)
set(_oecluster_swig_dir_gone FALSE)
if(_oecluster_swig_executable
        AND NOT _oecluster_swig_executable MATCHES "-NOTFOUND$"
        AND NOT EXISTS "${_oecluster_swig_executable}")
    set(_oecluster_swig_executable_gone TRUE)
endif()
if(_oecluster_swig_dir
        AND NOT _oecluster_swig_dir MATCHES "-NOTFOUND$"
        AND NOT EXISTS "${_oecluster_swig_dir}")
    set(_oecluster_swig_dir_gone TRUE)
endif()

if(_oecluster_swig_executable_gone)
    message(STATUS "Dropping cached SWIG_EXECUTABLE, which names a path that "
                   "no longer exists: ${_oecluster_swig_executable}")
    unset(SWIG_EXECUTABLE CACHE)
endif()
if(_oecluster_swig_dir_gone)
    message(STATUS "Dropping cached SWIG_DIR, which names a path that "
                   "no longer exists: ${_oecluster_swig_dir}")
elseif(_oecluster_swig_executable_gone AND _oecluster_swig_dir
        AND NOT _oecluster_swig_dir MATCHES "-NOTFOUND$")
    message(STATUS "Dropping cached SWIG_DIR, which belongs to the dropped "
                   "SWIG_EXECUTABLE: ${_oecluster_swig_dir}")
endif()
if(_oecluster_swig_executable_gone OR _oecluster_swig_dir_gone)
    unset(SWIG_DIR CACHE)
    unset(SWIG_VERSION CACHE)
endif()

unset(_oecluster_swig_executable)
unset(_oecluster_swig_dir)
unset(_oecluster_swig_executable_gone)
unset(_oecluster_swig_dir_gone)
