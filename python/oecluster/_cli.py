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


def main():
    binary = _binary_path()
    if not os.path.isfile(binary):
        print(f"Error: oepdist binary not found. Build with CMake first, or install the wheel.", file=sys.stderr)
        sys.exit(1)
    os.execv(binary, [binary] + sys.argv[1:])


if __name__ == "__main__":
    main()
