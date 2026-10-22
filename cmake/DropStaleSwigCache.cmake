# FindSWIG trusts the cache. It asks the executable for its library directory
# (swig -swiglib) only when SWIG_DIR is false in CMake's sense (such as unset,
# empty or a -NOTFOUND marker), and for its version (swig -version) only when
# SWIG_VERSION is false in the same sense. A cached SWIG_DIR or SWIG_VERSION
# therefore survives even when it no longer describes the cached executable:
# a library directory that a Homebrew SWIG upgrade deleted, the library of a
# different SWIG, or the version of a SWIG that is no longer the one found.
# UseSWIG passes SWIG_DIR to SWIG as SWIG_LIB, so SWIG generates the wrapper
# from whatever library SWIG_DIR names.
#
# This guard drops SWIG_DIR and SWIG_VERSION from the cache on every
# configure, and find_package(SWIG) derives both again: SWIG_VERSION from the
# executable, and SWIG_DIR by searching for swig.swg, with the directories the
# executable's -swiglib answer names as the last candidates. That costs two
# short SWIG runs per configure, and it means a SWIG_DIR in the cache,
# including one passed with -DSWIG_DIR, does not survive a configure. Only the
# cache entry is dropped: a normal variable named SWIG_DIR, such as a
# toolchain file or a parent project can set, is not touched and still wins.
#
# The guard does not change FindSWIG's search policy. That search looks in
# CMAKE_INCLUDE_PATH and CMake's other search locations before the -swiglib
# directories, so a SWIG library placed in one of those locations is chosen on
# any configure, fresh or not.
#
# SWIG_EXECUTABLE is dropped only when it names a path that no longer exists,
# so that find_package(SWIG) searches for SWIG again. A live executable,
# including one passed with -DSWIG_EXECUTABLE, is always kept.
#
# The guard runs only when CMake configures. CMake does not treat the SWIG
# installation as a regeneration input, so after a SWIG upgrade a plain
# `cmake --build` does not re-run CMake and keeps the stale SWIG_LIB in the
# generated rules. Re-run CMake on the build directory (`cmake <build-dir>`)
# to apply the guard.
set(_oecluster_swig_executable "$CACHE{SWIG_EXECUTABLE}")
if(_oecluster_swig_executable
        AND NOT _oecluster_swig_executable MATCHES "-NOTFOUND$"
        AND NOT EXISTS "${_oecluster_swig_executable}")
    message(STATUS "Dropping cached SWIG_EXECUTABLE, which names a path that "
                   "no longer exists: ${_oecluster_swig_executable}")
    unset(SWIG_EXECUTABLE CACHE)
endif()
unset(_oecluster_swig_executable)
unset(SWIG_DIR CACHE)
unset(SWIG_VERSION CACHE)
