"""Wrapper to invoke the bundled oepdist binary."""

import os
import sys


def _binary_path():
    pkg_dir = os.path.dirname(__file__)
    name = "oepdist.exe" if sys.platform == "win32" else "oepdist"
    return os.path.join(pkg_dir, "_bin", name)


def main():
    binary = _binary_path()
    if not os.path.isfile(binary):
        print(f"Error: oepdist binary not found at {binary}", file=sys.stderr)
        sys.exit(1)
    os.execv(binary, [binary] + sys.argv[1:])


if __name__ == "__main__":
    main()
