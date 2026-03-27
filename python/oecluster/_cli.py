"""Wrapper to invoke the bundled oepdist binary."""

import os
import re
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


def _ensure_compat_symlinks():
    """Create compatibility symlinks for version-mismatched OpenEye libraries.

    This is a lightweight version of oecluster._ensure_library_compat() that
    avoids importing the main package (which would trigger the SWIG module
    import and potentially segfault on version mismatch).
    """
    try:
        from oecluster import _build_info
    except ImportError:
        return

    if getattr(_build_info, 'OPENEYE_LIBRARY_TYPE', 'STATIC') != 'SHARED':
        return

    expected_libs = getattr(_build_info, 'OPENEYE_EXPECTED_LIBS', [])
    if not expected_libs:
        return

    try:
        from openeye import libs
        oe_lib_dir = libs.FindOpenEyeDLLSDirectory()
    except (ImportError, Exception):
        return

    if not os.path.isdir(oe_lib_dir):
        return

    pkg_dir = os.path.dirname(__file__)
    for expected_name in expected_libs:
        if os.path.exists(os.path.join(oe_lib_dir, expected_name)):
            continue

        symlink_path = os.path.join(pkg_dir, expected_name)
        if os.path.islink(symlink_path):
            if os.path.exists(symlink_path):
                continue
            try:
                os.unlink(symlink_path)
            except OSError:
                continue
        elif os.path.exists(symlink_path):
            continue

        match = re.match(r'(lib\w+?)(-[\d.]+)?(\.[\d.]*\w+)$', expected_name)
        if not match:
            continue

        base_name = match.group(1)
        actual_path = None
        for f in os.listdir(oe_lib_dir):
            if f.startswith(base_name + '-') or f.startswith(base_name + '.'):
                actual_path = os.path.join(oe_lib_dir, f)
                break

        if actual_path:
            try:
                os.symlink(actual_path, os.path.join(pkg_dir, expected_name))
            except OSError:
                pass


def _setup_library_env():
    """Set LD_LIBRARY_PATH so the binary can find OpenEye shared libraries."""
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

    _ensure_compat_symlinks()
    _setup_library_env()
    os.execv(binary, [binary] + sys.argv[1:])


if __name__ == "__main__":
    main()
