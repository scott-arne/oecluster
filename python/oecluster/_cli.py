"""Wrapper to invoke the bundled oepdist binary."""

import os
import sys


def _binary_path():
    name = "oepdist.exe" if sys.platform == "win32" else "oepdist"

    # Installed location (wheel)
    pkg_dir = os.path.dirname(__file__)
    installed = os.path.join(pkg_dir, "_bin", name)
    if os.path.isfile(installed):
        return installed

    # Development build directories (editable install)
    project_root = os.path.dirname(os.path.dirname(pkg_dir))
    for build_dir in ("build", "cmake-build-debug", "cmake-build-release"):
        dev_path = os.path.join(project_root, build_dir, "tools", name)
        if os.path.isfile(dev_path):
            return dev_path

    return installed


def _setup_library_env():
    """Set LD_LIBRARY_PATH so the binary can find OpenEye shared libraries.

    The SWIG module gets symlinks created by _ensure_library_compat() at import
    time, but the oepdist binary runs via execv and needs LD_LIBRARY_PATH to
    locate version-mismatched shared libraries.
    """
    pkg_dir = os.path.dirname(__file__)
    extra_paths = [pkg_dir]

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
        if os.path.isdir(oe_lib_dir):
            extra_paths.append(oe_lib_dir)
    except (ImportError, Exception):
        pass

    if sys.platform == "darwin":
        env_var = "DYLD_LIBRARY_PATH"
    else:
        env_var = "LD_LIBRARY_PATH"

    existing = os.environ.get(env_var, "")
    new_paths = os.pathsep.join(extra_paths)
    if existing:
        os.environ[env_var] = new_paths + os.pathsep + existing
    else:
        os.environ[env_var] = new_paths


def main():
    binary = _binary_path()
    if not os.path.isfile(binary):
        print(f"Error: oepdist binary not found. Build with CMake first, or install the wheel.", file=sys.stderr)
        sys.exit(1)

    # Ensure compat symlinks exist and set library search paths
    try:
        from oecluster import _ensure_library_compat
        _ensure_library_compat()
    except ImportError:
        pass
    _setup_library_env()

    os.execv(binary, [binary] + sys.argv[1:])


if __name__ == "__main__":
    main()
