"""The build drops a cached SWIG path that no longer exists.

A homebrew SWIG bump deletes the versioned Cellar directory CMakeCache.txt
recorded as SWIG_DIR, and FindSWIG neither re-checks a cached path nor
recomputes a cached SWIG_VERSION. cmake/DropStaleSwigCache.cmake clears all
three when either path has vanished, so find_package(SWIG) detects afresh.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
GUARD = REPO_ROOT / "cmake" / "DropStaleSwigCache.cmake"

pytestmark = pytest.mark.skipif(
    shutil.which("cmake") is None or shutil.which("swig") is None,
    reason="needs cmake and swig on PATH")


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


def test_both_paths_vanished(tmp_path):
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE="/nonexistent/swig",
        SWIG_DIR="/nonexistent/swig-lib", SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_EXECUTABLE"]).exists()
    assert Path(cached["SWIG_DIR"]).is_dir()
    assert cached["SWIG_VERSION"] != "99.0.0"
    assert "/nonexistent/swig" in result.stdout


def test_only_the_library_directory_vanished(tmp_path):
    """The broken scikit-build cache in this repository: a live executable,
    a vanished SWIG_DIR and a stale version. FindSWIG skips rediscovery while
    SWIG_DIR is set, so a guard that checked only the executable fails here."""
    result, cached = _configure(
        tmp_path, SWIG_EXECUTABLE=shutil.which("swig"),
        SWIG_DIR="/nonexistent/swig-lib", SWIG_VERSION="99.0.0")
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(cached["SWIG_DIR"]).is_dir()
    assert cached["SWIG_VERSION"] != "99.0.0"
