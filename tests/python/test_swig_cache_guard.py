"""The build re-derives SWIG's library path and version on every configure.

FindSWIG trusts a cached SWIG_DIR and SWIG_VERSION whatever executable they
came from, and a cached SWIG_EXECUTABLE whether or not it still exists, so a
homebrew SWIG bump that deletes the versioned Cellar directory recorded as
SWIG_DIR breaks the build. cmake/DropStaleSwigCache.cmake drops SWIG_DIR and
SWIG_VERSION before every find_package(SWIG), so FindSWIG derives both again,
and drops SWIG_EXECUTABLE only when its path has vanished.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GUARD = REPO_ROOT / "cmake" / "DropStaleSwigCache.cmake"
# The probe project asks for the version the stable-ABI build needs.
MIN_SWIG = (4, 4, 1)
CACHE_TYPES = {"SWIG_EXECUTABLE": "FILEPATH", "SWIG_DIR": "PATH", "SWIG_VERSION": "STRING"}


def _swig_version(swig):
    """The version ``swig -version`` reports, read the way FindSWIG reads it."""
    result = subprocess.run([str(swig), "-version"], capture_output=True,
                            text=True, timeout=60, check=False)
    match = re.search(r"SWIG Version[^0-9.]*([0-9.]+)", result.stdout + result.stderr)
    return match.group(1) if match else None


def _too_old(version):
    if version is None:
        return True
    return tuple(int(part) for part in version.split(".") if part) < MIN_SWIG


SWIG_ON_PATH = shutil.which("swig")
FOUND_SWIG = _swig_version(SWIG_ON_PATH) if SWIG_ON_PATH else None

pytestmark = [
    pytest.mark.skipif(
        shutil.which("cmake") is None or SWIG_ON_PATH is None,
        reason="needs cmake and swig on PATH"),
    pytest.mark.skipif(
        SWIG_ON_PATH is not None and _too_old(FOUND_SWIG),
        reason=f"needs SWIG {'.'.join(map(str, MIN_SWIG))} or newer on PATH, "
               f"found {FOUND_SWIG or 'an unreadable version'}"),
]


def _configure(tmp_path, **cache):
    project = tmp_path / "project"
    project.mkdir()
    (project / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.16)\n"
        "project(guard_probe NONE)\n"
        f'include("{GUARD.as_posix()}")\n'
        "find_package(SWIG 4.4.1 REQUIRED)\n")
    build = tmp_path / "build"
    args = ["cmake", "-S", str(project), "-B", str(build)]
    # Typed, as CMakeCache.txt holds them. An untyped -D does not reproduce a
    # cached -NOTFOUND marker: find_program leaves an untyped
    # -DSWIG_EXECUTABLE=SWIG_EXECUTABLE-NOTFOUND unresolved.
    args += [f"-D{key}:{CACHE_TYPES[key]}={value}" for key, value in cache.items()]
    result = subprocess.run(args, capture_output=True, text=True, timeout=120,
                            check=False)
    cached = {}
    cache_file = build / "CMakeCache.txt"
    if cache_file.exists():
        for line in cache_file.read_text().splitlines():
            match = re.match(r"^(SWIG_EXECUTABLE|SWIG_DIR|SWIG_VERSION):[A-Z]+=(.*)$", line)
            if match:
                cached[match.group(1)] = match.group(2)
    return result, cached


def _dropped(result):
    """Map each variable the guard reported dropping to the path it named."""
    dropped = {}
    for line in result.stdout.splitlines():
        match = re.match(r"^-- Dropping cached (SWIG_\w+), .*: (.+)$", line)
        if match:
            dropped[match.group(1)] = Path(match.group(2))
    return dropped


def _swiglib(swig):
    """The library directories ``swig -swiglib`` reports, as FindSWIG reads them."""
    output = subprocess.run([str(swig), "-swiglib"], capture_output=True,
                            text=True, timeout=60, check=True).stdout
    return {Path(line.strip()) for line in output.splitlines() if line.strip()}


def _assert_derived_from_executable(cached):
    """SWIG_DIR and SWIG_VERSION describe the cached SWIG_EXECUTABLE."""
    swig = cached["SWIG_EXECUTABLE"]
    assert Path(swig).exists()
    assert Path(cached["SWIG_DIR"]) in _swiglib(swig)
    assert cached["SWIG_VERSION"] == _swig_version(swig)


def _user_swig(tmp_path):
    """A live SWIG executable at a path find_program would never choose, the
    way a user passes -DSWIG_EXECUTABLE for a SWIG outside PATH."""
    found = Path(SWIG_ON_PATH)
    user_swig = tmp_path / "user-bin" / found.name
    user_swig.parent.mkdir()
    try:
        user_swig.symlink_to(found)
    except OSError:
        # Creating a symlink can need a privilege on Windows. A copy is still a
        # live executable outside PATH, which is all these cases need.
        shutil.copy2(found, user_swig)
    return user_swig


def test_both_paths_vanished(tmp_path):
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="/nonexistent/swig",
        SWIG_DIR="/nonexistent/swig-lib", SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    _assert_derived_from_executable(cached)
    assert _dropped(result) == {"SWIG_EXECUTABLE": Path("/nonexistent/swig")}


def test_only_the_library_directory_vanished(tmp_path):
    """A live executable with only the library directory gone, the shape a
    Homebrew SWIG upgrade leaves in a cache. FindSWIG skips rediscovery while
    SWIG_DIR is set, so a guard that checked only the executable fails here,
    and the executable the user chose must survive the repair."""
    user_swig = _user_swig(tmp_path)
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE=user_swig.as_posix(),
        SWIG_DIR="/nonexistent/swig-lib", SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_EXECUTABLE"]) == user_swig
    _assert_derived_from_executable(cached)
    assert _dropped(result) == {}


def test_only_the_executable_vanished(tmp_path):
    """A vanished executable does not leave its library directory behind,
    even one that still exists: the re-detected executable gets its own
    library. The cached directory is a live decoy so that keeping it shows."""
    old_library = tmp_path / "old-swig-lib"
    old_library.mkdir()
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="/nonexistent/swig",
        SWIG_DIR=old_library.as_posix(), SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    _assert_derived_from_executable(cached)
    assert _dropped(result) == {"SWIG_EXECUTABLE": Path("/nonexistent/swig")}


def test_live_library_directory_of_another_swig(tmp_path):
    """Both cached paths exist, but SWIG_DIR is not the executable's library,
    as when the cache names an older SWIG's library. Existence alone cannot
    tell, so both are dropped and derived again, here from the executable,
    which the user chose and which survives."""
    user_swig = _user_swig(tmp_path)
    other_library = tmp_path / "other-swig-lib"
    other_library.mkdir()
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE=user_swig.as_posix(),
        SWIG_DIR=other_library.as_posix(), SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_EXECUTABLE"]) == user_swig
    _assert_derived_from_executable(cached)
    assert _dropped(result) == {}


def test_notfound_markers_with_a_stale_version(tmp_path):
    """FindSWIG searches again for a -NOTFOUND path but keeps a cached
    SWIG_VERSION, so the version check would judge the SWIG it finds by a
    version it never reported."""
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="SWIG_EXECUTABLE-NOTFOUND",
        SWIG_DIR="SWIG_DIR-NOTFOUND", SWIG_VERSION="1.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    _assert_derived_from_executable(cached)
    assert _dropped(result) == {}
