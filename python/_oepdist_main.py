"""Entry point for the oepdist CLI.

This module lives outside the oecluster package so that the console_scripts
entry point can invoke it without importing oecluster (which loads the SWIG
extension and can segfault on OpenEye version mismatch).
"""

import os
import re
import sys


def _find_binary():
    name = "oepdist.exe" if sys.platform == "win32" else "oepdist"

    # Installed location (wheel): _oepdist_main.py sits next to oecluster/
    pkg_dir = os.path.join(os.path.dirname(__file__), "oecluster")
    installed = os.path.join(pkg_dir, "_bin", name)
    if os.path.isfile(installed):
        return installed, pkg_dir

    # Development build directories
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for build_dir in ("build", "cmake-build-debug", "cmake-build-release"):
        dev_path = os.path.join(project_root, build_dir, "tools", name)
        if os.path.isfile(dev_path):
            return dev_path, pkg_dir

    return installed, pkg_dir


def _ensure_compat_symlinks(pkg_dir):
    """Create compatibility symlinks for version-mismatched OpenEye libraries."""
    try:
        # Import only _build_info (a plain data module) via importlib to avoid
        # triggering oecluster/__init__.py and the SWIG extension import.
        import importlib.util
        bi_path = os.path.join(pkg_dir, "_build_info.py")
        if not os.path.isfile(bi_path):
            return
        spec = importlib.util.spec_from_file_location("_build_info", bi_path)
        _build_info = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_build_info)
    except Exception:
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


def _setup_library_env(pkg_dir):
    """Set LD_LIBRARY_PATH so the binary can find OpenEye shared libraries."""
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
    binary, pkg_dir = _find_binary()
    if not os.path.isfile(binary):
        print("Error: oepdist binary not found. Build with CMake first, or install the wheel.", file=sys.stderr)
        sys.exit(1)

    _ensure_compat_symlinks(pkg_dir)
    _setup_library_env(pkg_dir)
    os.execv(binary, [binary] + sys.argv[1:])


if __name__ == "__main__":
    main()
