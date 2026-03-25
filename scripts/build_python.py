#!/usr/bin/env python3
"""
Build and package oecluster for PyPI distribution.

Usage:
    python scripts/build_python.py [options]

Options:
    --openeye-root PATH    Path to OpenEye C++ SDK (headers)
    --python PATH          Python executable to use
    --clean                Clean dist/ before building
    --upload               Upload to PyPI after building (requires twine)
    --test-upload          Upload to TestPyPI instead of PyPI
    --verbose              Verbose output
"""

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    NC = '\033[0m'

    @classmethod
    def disable(cls):
        cls.RED = cls.GREEN = cls.YELLOW = cls.BLUE = cls.NC = ''


def print_header(msg):
    print(f"\n{Colors.GREEN}{'=' * 60}{Colors.NC}")
    print(f"{Colors.GREEN}{msg}{Colors.NC}")
    print(f"{Colors.GREEN}{'=' * 60}{Colors.NC}\n")

def print_step(msg):
    print(f"{Colors.YELLOW}>>> {msg}{Colors.NC}")

def print_error(msg):
    print(f"{Colors.RED}ERROR: {msg}{Colors.NC}", file=sys.stderr)

def print_success(msg):
    print(f"{Colors.GREEN}{msg}{Colors.NC}")

def run_command(cmd, cwd=None, check=True, capture_output=False, verbose=False):
    if verbose:
        print(f"{Colors.BLUE}Running: {' '.join(cmd)}{Colors.NC}")
    return subprocess.run(cmd, cwd=cwd, check=check, capture_output=capture_output, text=True)

def get_openeye_info(python_exe):
    code = """
from openeye import libs, oechem
import os
dll_dir = libs.FindOpenEyeDLLSDirectory()
version = oechem.OEToolkitsGetRelease()
print(f'VERSION:{version}')
print(f'LIB_DIR:{dll_dir}')
print(f'PLATFORM:{os.path.basename(dll_dir)}')
"""
    try:
        result = run_command([python_exe, '-c', code], capture_output=True, check=True)
        info = {}
        for line in result.stdout.strip().split('\n'):
            if ':' in line:
                key, value = line.split(':', 1)
                info[key] = value
        return info
    except subprocess.CalledProcessError as e:
        if e.stderr:
            print_step(f"openeye import error: {e.stderr.strip()}")
        return None

def verify_openeye_root(openeye_root):
    include_dir = Path(openeye_root) / 'include'
    if not include_dir.exists():
        print_error(f"OpenEye include directory not found: {include_dir}")
        return False
    if not (include_dir / 'oechem.h').exists():
        print_error(f"oechem.h not found in {include_dir}")
        return False
    return True

def get_openeye_root_from_cmake_presets(project_dir):
    all_presets = {}
    for filename in ("CMakePresets.json", "CMakeUserPresets.json"):
        filepath = project_dir / filename
        if not filepath.exists():
            continue
        try:
            data = json.loads(filepath.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        for preset in data.get("configurePresets", []):
            name = preset.get("name")
            if name:
                all_presets[name] = preset

    def resolve_cache_variables(preset_name, visited=None):
        if visited is None:
            visited = set()
        if preset_name in visited or preset_name not in all_presets:
            return {}
        visited.add(preset_name)
        entry = all_presets[preset_name]
        inherits = entry.get("inherits", [])
        if isinstance(inherits, str):
            inherits = [inherits]
        merged = {}
        for parent_name in inherits:
            merged.update(resolve_cache_variables(parent_name, visited))
        merged.update(entry.get("cacheVariables", {}))
        return merged

    for hidden in (False, True):
        for name, preset in all_presets.items():
            if preset.get("hidden", False) != hidden:
                continue
            cache_vars = resolve_cache_variables(name)
            value = cache_vars.get("OPENEYE_ROOT")
            if value and isinstance(value, str):
                return value
    return None

def get_version_from_pyproject(pyproject_path):
    if not pyproject_path.exists():
        return None
    content = pyproject_path.read_text()
    match = re.search(r'version\s*=\s*"([^"]+)"', content)
    return match.group(1) if match else None

def get_lib_version(project_dir):
    lib_pyproject = Path(project_dir) / 'pyproject.toml'
    lib_version = get_version_from_pyproject(lib_pyproject)
    if lib_version is None:
        print_error(f"Could not extract version from {lib_pyproject}")
    return lib_version

def find_missing_libraries(lib_dir, verbose=False):
    """Find OpenEye libraries that are missing from the library directory."""
    expected_missing = {
        'liboebio',
        'liboehermite',
        'liboeopt',
        'liboesitehopper',
        'liboeshape',
        'liboegraphsim',
    }
    system = platform.system()
    if system == 'Darwin':
        ext = '.dylib'
    elif system == 'Linux':
        ext = '.so'
    elif system == 'Windows':
        ext = '.dll'
    else:
        print_error(f"Unsupported platform: {system}")
        return []

    missing = []
    lib_path = Path(lib_dir)
    if not lib_path.exists():
        print_error(f"Library directory does not exist: {lib_dir}")
        return []

    lib_files = [f.name for f in lib_path.iterdir() if ext in f.suffixes or f.name.endswith(ext)]
    if verbose:
        print_step(f"Found {len(lib_files)} libraries in {lib_dir}")

    for lib_name in expected_missing:
        if not any(f.startswith(lib_name + '-') or f.startswith(lib_name + '.') for f in lib_files):
            missing.append(lib_name)

    return missing

def build_oecluster(python_exe, openeye_root, project_dir, verbose=False):
    """Build oecluster wheel using pip."""
    print_step("Building oecluster wheel")

    cmd = [
        str(python_exe), '-m', 'pip', 'wheel',
        '--no-build-isolation',
        '--no-deps',
        '--wheel-dir', 'dist',
        '-C', 'cmake.define.OECLUSTER_BUILD_TESTS=OFF',
        '-C', f'cmake.define.OPENEYE_ROOT={openeye_root}',
        '.'
    ]

    run_command(cmd, cwd=project_dir, verbose=verbose)
    print_success("Build completed successfully")

def fix_rpath_and_sign(wheel_path, verbose=False):
    """Fix RPATH on macOS and sign the binary."""
    system = platform.system()
    if system != 'Darwin':
        return

    print_step("Fixing RPATH and signing binary on macOS")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Unzip wheel
        run_command(['unzip', '-q', str(wheel_path)], cwd=tmpdir, verbose=verbose)

        # Find the .so file
        so_file = tmpdir / 'oecluster' / '_oecluster.so'
        if not so_file.exists():
            print_error(f"Binary not found: {so_file}")
            return

        # Get current RPATH
        result = run_command(
            ['otool', '-l', str(so_file)],
            capture_output=True,
            verbose=verbose
        )

        # Fix RPATH if needed
        if '@loader_path' not in result.stdout:
            print_step("Adding @loader_path to RPATH")
            run_command(
                ['install_name_tool', '-add_rpath', '@loader_path', str(so_file)],
                verbose=verbose
            )

        # Sign the binary
        print_step("Signing binary")
        run_command(
            ['codesign', '--force', '--sign', '-', str(so_file)],
            verbose=verbose
        )

        # Repackage wheel
        print_step("Repackaging wheel")
        wheel_path.unlink()

        # Get wheel contents
        wheel_contents = list(tmpdir.glob('*'))

        # Create new wheel
        run_command(
            ['zip', '-q', '-r', str(wheel_path)] + [c.name for c in wheel_contents],
            cwd=tmpdir,
            verbose=verbose
        )

    print_success("RPATH and signing completed")

def upload_to_pypi(dist_dir, test_pypi=False, verbose=False):
    """Upload wheels to PyPI using twine."""
    repo = 'testpypi' if test_pypi else 'pypi'
    print_step(f"Uploading to {repo}")

    cmd = ['twine', 'upload', '--repository', repo, str(dist_dir / '*.whl')]
    run_command(cmd, verbose=verbose)
    print_success(f"Upload to {repo} completed")

def main():
    parser = argparse.ArgumentParser(
        description='Build and package oecluster for PyPI distribution',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--openeye-root', type=str, help='Path to OpenEye C++ SDK (headers)')
    parser.add_argument('--python', type=str, default=sys.executable, help='Python executable to use')
    parser.add_argument('--clean', action='store_true', help='Clean dist/ before building')
    parser.add_argument('--upload', action='store_true', help='Upload to PyPI after building')
    parser.add_argument('--test-upload', action='store_true', help='Upload to TestPyPI instead of PyPI')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--no-color', action='store_true', help='Disable colored output')

    args = parser.parse_args()

    if args.no_color:
        Colors.disable()

    project_dir = Path(__file__).parent.parent.absolute()
    dist_dir = project_dir / 'dist'

    print_header("oecluster Build Script")

    # Get Python executable
    python_exe = Path(args.python).absolute()
    if not python_exe.exists():
        print_error(f"Python executable not found: {python_exe}")
        return 1

    print_step(f"Using Python: {python_exe}")

    # Get library version
    lib_version = get_lib_version(project_dir)
    if lib_version:
        print_step(f"Library version: {lib_version}")

    # Get OpenEye information
    print_step("Detecting OpenEye configuration")
    oe_info = get_openeye_info(python_exe)

    if oe_info:
        print_step(f"OpenEye version: {oe_info['VERSION']}")
        print_step(f"OpenEye library directory: {oe_info['LIB_DIR']}")
        print_step(f"Platform: {oe_info['PLATFORM']}")

        # Check for missing libraries
        missing = find_missing_libraries(oe_info['LIB_DIR'], verbose=args.verbose)
        if missing:
            print_error(f"Missing expected libraries: {', '.join(missing)}")
    else:
        print_error("Could not detect OpenEye configuration")
        print_step("Make sure openeye-toolkits is installed in the Python environment")
        return 1

    # Get OpenEye root
    openeye_root = args.openeye_root
    if not openeye_root:
        openeye_root = get_openeye_root_from_cmake_presets(project_dir)
        if openeye_root:
            print_step(f"Found OPENEYE_ROOT from CMake presets: {openeye_root}")

    if not openeye_root:
        print_error("Could not determine OPENEYE_ROOT")
        print_step("Please specify --openeye-root or configure CMakeUserPresets.json")
        return 1

    if not verify_openeye_root(openeye_root):
        return 1

    print_step(f"Using OPENEYE_ROOT: {openeye_root}")

    # Clean dist directory if requested
    if args.clean and dist_dir.exists():
        print_step("Cleaning dist/ directory")
        shutil.rmtree(dist_dir)

    # Create dist directory
    dist_dir.mkdir(exist_ok=True)

    # Build wheel
    try:
        build_oecluster(python_exe, openeye_root, project_dir, verbose=args.verbose)
    except subprocess.CalledProcessError as e:
        print_error("Build failed")
        return 1

    # Find the built wheel
    wheels = list(dist_dir.glob('oecluster-*.whl'))
    if not wheels:
        print_error("No wheel found in dist/ directory")
        return 1

    wheel_path = wheels[0]
    print_success(f"Built wheel: {wheel_path.name}")

    # Fix RPATH and sign on macOS
    if platform.system() == 'Darwin':
        try:
            fix_rpath_and_sign(wheel_path, verbose=args.verbose)
        except subprocess.CalledProcessError as e:
            print_error("Failed to fix RPATH or sign binary")
            return 1

    # Upload if requested
    if args.upload or args.test_upload:
        try:
            upload_to_pypi(dist_dir, test_pypi=args.test_upload, verbose=args.verbose)
        except subprocess.CalledProcessError as e:
            print_error("Upload failed")
            return 1

    print_header("Build Complete")
    return 0


if __name__ == '__main__':
    sys.exit(main())
