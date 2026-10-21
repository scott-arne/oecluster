"""The build drops a cached SWIG path that no longer exists.

A homebrew SWIG bump deletes the versioned Cellar directory CMakeCache.txt
recorded as SWIG_DIR, and FindSWIG neither re-checks a cached path nor
recomputes a cached SWIG_VERSION. cmake/DropStaleSwigCache.cmake drops a
vanished SWIG_EXECUTABLE together with SWIG_DIR and SWIG_VERSION, and a
vanished SWIG_DIR together with SWIG_VERSION while keeping a live executable,
so find_package(SWIG) detects afresh what it needs to.
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


def _swig_version():
    swig = shutil.which("swig")
    if swig is None:
        return None
    output = subprocess.run([swig, "-version"], capture_output=True, text=True,
                            timeout=60, check=False).stdout
    match = re.search(r"SWIG Version (\d+)\.(\d+)\.(\d+)", output)
    return tuple(int(part) for part in match.groups()) if match else (0, 0, 0)


FOUND_SWIG = _swig_version()

pytestmark = [
    pytest.mark.skipif(
        shutil.which("cmake") is None or FOUND_SWIG is None,
        reason="needs cmake and swig on PATH"),
    pytest.mark.skipif(
        FOUND_SWIG is not None and FOUND_SWIG < MIN_SWIG,
        reason=f"needs SWIG {'.'.join(map(str, MIN_SWIG))} or newer on PATH, "
               f"found {'.'.join(map(str, FOUND_SWIG or ()))}"),
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
    args += [f"-D{key}={value}" for key, value in cache.items()]
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
        match = re.match(r"^-- Dropping cached (SWIG_EXECUTABLE|SWIG_DIR), .*: (.+)$", line)
        if match:
            dropped[match.group(1)] = Path(match.group(2))
    return dropped


def _swiglib(swig):
    """The library directories ``swig -swiglib`` reports, as FindSWIG reads them."""
    output = subprocess.run([swig, "-swiglib"], capture_output=True, text=True,
                            timeout=60, check=True).stdout
    return {Path(line.strip()) for line in output.splitlines() if line.strip()}


def _user_swig(tmp_path):
    """A live SWIG executable at a path find_program would never choose, the
    way a user passes -DSWIG_EXECUTABLE for a SWIG outside PATH."""
    found = Path(shutil.which("swig"))
    user_swig = tmp_path / "user-bin" / found.name
    user_swig.parent.mkdir()
    try:
        user_swig.symlink_to(found)
    except OSError:
        # Creating a symlink can need a privilege on Windows. A copy is still a
        # live executable outside PATH, which is all this case needs.
        shutil.copy2(found, user_swig)
    return user_swig


def test_both_paths_vanished(tmp_path):
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="/nonexistent/swig",
        SWIG_DIR="/nonexistent/swig-lib", SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_EXECUTABLE"]).exists()
    assert Path(cached["SWIG_DIR"]).is_dir()
    assert cached["SWIG_VERSION"] != "99.0.0"
    assert _dropped(result) == {
        "SWIG_EXECUTABLE": Path("/nonexistent/swig"),
        "SWIG_DIR": Path("/nonexistent/swig-lib")}


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
    assert Path(cached["SWIG_DIR"]) in _swiglib(user_swig)
    assert cached["SWIG_VERSION"] != "99.0.0"
    assert _dropped(result) == {"SWIG_DIR": Path("/nonexistent/swig-lib")}


def test_only_the_executable_vanished(tmp_path):
    """A vanished executable takes its library directory with it, even one
    that still exists: the re-detected executable gets its own library, not
    the old one. The cached directory is a live decoy so that keeping it
    shows."""
    old_library = tmp_path / "old-swig-lib"
    old_library.mkdir()
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="/nonexistent/swig",
        SWIG_DIR=old_library.as_posix(), SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_EXECUTABLE"]).exists()
    assert Path(cached["SWIG_DIR"]) in _swiglib(cached["SWIG_EXECUTABLE"])
    assert cached["SWIG_VERSION"] != "99.0.0"
    assert _dropped(result) == {
        "SWIG_EXECUTABLE": Path("/nonexistent/swig"),
        "SWIG_DIR": old_library}
