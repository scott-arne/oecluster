"""Release-gate checks for the version strings and the OEFP dependency pin.

``vrzn`` keeps the version locations in sync during a bump; these assertions
catch a hand-edit that bypassed it, and the SWIG module dunder in particular,
which drifted to 4.2.0 while the package was at 4.2.3. The OEFP pin is not
managed by ``vrzn`` at all, so the pin sites are checked against one literal.
"""

import re
import tomllib
from pathlib import Path

import oecluster

REPO_ROOT = Path(__file__).resolve().parents[2]

# Not a bare floor: oecluster compiles OEFP's core into its own extension and
# exchanges batch pointers with the installed wheel, so the two must share a
# minor series and the upper bound is part of the contract.
OEFP_REQUIREMENT = "oefp>=0.3.0,<0.4"
OEFP_SERIES = "0.3"

# The floor of the requirement is the tag that gets compiled in, and its minor
# series is what the README documents. Deriving both from the one requirement
# string keeps the four sites from drifting apart.
OEFP_FLOOR = OEFP_REQUIREMENT.split(">=", 1)[1].split(",", 1)[0]


def _pyproject_oefp_requirement(relative_path: str) -> str:
    """Return the single ``oefp`` entry from a pyproject's dependencies.

    :param relative_path: Path to the pyproject file, relative to the
        repository root.
    :returns: The requirement string with surrounding whitespace stripped.
    :raises AssertionError: If the file does not declare exactly one ``oefp``
        dependency.
    """
    with (REPO_ROOT / relative_path).open("rb") as handle:
        data = tomllib.load(handle)
    dependencies = data["project"]["dependencies"]
    matches = [dep.strip() for dep in dependencies if dep.strip().startswith("oefp")]
    assert len(matches) == 1, f"{relative_path} declares {len(matches)} oefp deps"
    return matches[0]


def _c_macro_version(relative_path: str) -> str:
    """Return the dotted version spelled by a file's ``OECLUSTER_VERSION_*`` macros.

    :param relative_path: Path to the file, relative to the repository root.
    :returns: The major, minor and patch macro values joined with dots.
    :raises AssertionError: If any of the three macros is missing.
    """
    text = (REPO_ROOT / relative_path).read_text()
    components = []
    for component in ("MAJOR", "MINOR", "PATCH"):
        match = re.search(rf"#define OECLUSTER_VERSION_{component}\s+(\d+)", text)
        assert match is not None, f"no OECLUSTER_VERSION_{component} in {relative_path}"
        components.append(match.group(1))
    return ".".join(components)


class TestOEFPPin:
    """The four OEFP pin sites must agree on one requirement."""

    def test_root_pyproject(self):
        """The distribution requires the pinned OEFP range."""
        assert _pyproject_oefp_requirement("pyproject.toml") == OEFP_REQUIREMENT

    def test_python_pyproject(self):
        """The Python-only package requires the same range as the root project."""
        assert _pyproject_oefp_requirement("python/pyproject.toml") == OEFP_REQUIREMENT

    def test_cmake_fetchcontent_tag(self):
        """The vendored OEFP source is fetched at the pinned floor."""
        text = (REPO_ROOT / "CMakeLists.txt").read_text()
        match = re.search(
            r"GIT_REPOSITORY\s+\S*oefp\.git\s+GIT_TAG\s+v(\d+\.\d+\.\d+)", text
        )
        assert match is not None, "no oefp FetchContent_Declare in CMakeLists.txt"
        assert match.group(1) == OEFP_FLOOR

    def test_readme_requirements_line(self):
        """The documented OEFP series matches the pinned minor series."""
        assert OEFP_SERIES == OEFP_FLOOR.rsplit(".", 1)[0]
        text = (REPO_ROOT / "README.md").read_text()
        match = re.search(r"\*\*OEFP\*\* (\d+\.\d+)\.x", text)
        assert match is not None, "no OEFP requirements line in README.md"
        assert match.group(1) == OEFP_SERIES


class TestVersionPins:
    """Every version location must agree with ``oecluster.__version__``."""

    def test_version_info_tuple(self):
        """The tuple form spells the same version as the string form."""
        expected = tuple(int(part) for part in oecluster.__version__.split("."))
        assert oecluster.__version_info__ == expected

    def test_cmake_project_version(self):
        """The CMake project version matches the package version."""
        text = (REPO_ROOT / "CMakeLists.txt").read_text()
        match = re.search(r"project\(oecluster VERSION (\d+\.\d+\.\d+)", text)
        assert match is not None, "no project(oecluster VERSION ...) in CMakeLists.txt"
        assert match.group(1) == oecluster.__version__

    def test_umbrella_header_macros(self):
        """The umbrella header macros match the package version."""
        assert _c_macro_version("include/oecluster/oecluster.h") == oecluster.__version__

    def test_swig_interface_macros(self):
        """The SWIG interface macros match the package version."""
        assert _c_macro_version("swig/oecluster.i") == oecluster.__version__

    def test_swig_module_dunder(self):
        """The generated extension module reports the package version.

        This location drifted to 4.2.0 while the package was at 4.2.3 because it
        was registered with ``vrzn`` only as a C macro block, which does not see
        the Python string below it. Reads the ``.i`` source rather than the
        generated ``oecluster.py``, so it holds without a rebuild.
        """
        text = (REPO_ROOT / "swig/oecluster.i").read_text()
        match = re.search(r'__version__ = "(\d+\.\d+\.\d+)"', text)
        assert match is not None, "no __version__ assignment in swig/oecluster.i"
        assert match.group(1) == oecluster.__version__

    def test_changelog_documents_this_version(self):
        """The release notes name the version being released."""
        text = (REPO_ROOT / "CHANGELOG.md").read_text()
        assert f"## [{oecluster.__version__}]" in text
