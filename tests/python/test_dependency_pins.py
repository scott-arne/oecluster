"""Release-gate checks for the version strings and the OEFP dependency pin.

``vrzn`` keeps the version locations in sync during a bump; these assertions
catch a hand-edit that bypassed it, the SWIG module dunder in particular, which
drifted to 4.2.0 while the package was at 4.2.3. ``TestVersionPins`` records
which locations it checks and which are covered elsewhere. The OEFP pin is not
managed by ``vrzn`` at all, so the pin sites are checked against one literal.
``TestManifestAgreement`` covers the rest of the requirement metadata, which the
two pyproject files must spell identically.
"""

import re
import tomllib
from pathlib import Path

import oecluster

REPO_ROOT = Path(__file__).resolve().parents[2]

# An exact pin, not a range: oecluster compiles OEFP's core into its own
# extension and exchanges raw fingerprint-batch pointers with the separately
# compiled wheel, so the compiled-against and installed versions must be
# identical. The extension enforces that itself the first time such a pointer
# crosses -- see tests/python/test_oefp_abi_guard.py -- and a range pin would
# admit a wheel that check then refuses.
OEFP_REQUIREMENT = "oefp==0.3.0"

# The pinned version is the tag that gets compiled in and the version the
# README names. Deriving it from the one requirement string keeps the four
# sites from drifting apart.
OEFP_VERSION = OEFP_REQUIREMENT.split("==", 1)[1]

# Publication metadata the root manifest alone is permitted, because the root
# project builds the distributed wheel through scikit-build-core and so owns it;
# python/pyproject.toml is the pure-Python manifest and declares only what a
# consumer resolves against. These names are excused on the root and nowhere
# else: the consumer manifest carrying one is itself a defect and fails the
# comparison. A field that starts appearing on the root alone belongs here with
# its own reason, deliberately admitted, rather than slipping past unnoticed.
ROOT_ONLY_PROJECT_FIELDS = frozenset(
    {"authors", "classifiers", "keywords", "license", "license-files"}
)

# The documentation toolchain, permitted on the root manifest alone because the
# docs build runs from the repository root. Excused there and nowhere else: the
# consumer manifest declaring it is a defect and fails the comparison. A new
# root-only extra belongs here with its own reason, deliberately admitted,
# rather than slipping past unnoticed.
ROOT_ONLY_EXTRAS = frozenset({"docs"})


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


def _pyproject_version(relative_path: str) -> str:
    """Return a pyproject's declared ``project.version``.

    :param relative_path: Path to the pyproject file, relative to the
        repository root.
    :returns: The version string as written in the file.
    """
    with (REPO_ROOT / relative_path).open("rb") as handle:
        data = tomllib.load(handle)
    return data["project"]["version"]


def _pyproject_dependencies(relative_path: str) -> list[str]:
    """Return a pyproject's declared ``project.dependencies``.

    :param relative_path: Path to the pyproject file, relative to the
        repository root.
    :returns: The runtime requirement strings, in declaration order.
    """
    with (REPO_ROOT / relative_path).open("rb") as handle:
        data = tomllib.load(handle)
    return data["project"]["dependencies"]


def _pyproject_dev_extra(relative_path: str) -> list[str]:
    """Return a pyproject's ``dev`` optional-dependency group.

    :param relative_path: Path to the pyproject file, relative to the
        repository root.
    :returns: The ``dev`` extra's requirement strings, in declaration order.
    :raises KeyError: If the file declares no ``dev`` extra.
    """
    with (REPO_ROOT / relative_path).open("rb") as handle:
        data = tomllib.load(handle)
    return data["project"]["optional-dependencies"]["dev"]


def _pyproject_project_table(
    relative_path: str, *, drop_root_only: bool
) -> dict[str, object]:
    """Return a pyproject's ``project`` table, ready to compare with the other's.

    The reduction is a subtraction rather than an enumeration: everything the
    two manifests are expected to spell identically is kept. ``optional-
    dependencies`` is normalized to always be present, so a manifest that
    declares no extras still compares against one that does.

    :param relative_path: Path to the pyproject file, relative to the
        repository root.
    :param drop_root_only: Whether to subtract ``ROOT_ONLY_PROJECT_FIELDS`` and
        ``ROOT_ONLY_EXTRAS``. True for the root manifest, the only one those
        names are excused on; False for the consumer manifest, so that it
        acquiring one of them surfaces as a difference rather than being
        ignored on both sides.
    :returns: A fresh dictionary holding the table. The parsed file is re-read
        on each call, so callers cannot affect one another.
    """
    with (REPO_ROOT / relative_path).open("rb") as handle:
        data = tomllib.load(handle)
    project = dict(data["project"])
    extras = dict(project.setdefault("optional-dependencies", {}))
    if drop_root_only:
        for field in ROOT_ONLY_PROJECT_FIELDS:
            project.pop(field, None)
        extras = {
            name: requirements
            for name, requirements in extras.items()
            if name not in ROOT_ONLY_EXTRAS
        }
    project["optional-dependencies"] = extras
    return project


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
        """The vendored OEFP source is fetched at the pinned version."""
        text = (REPO_ROOT / "CMakeLists.txt").read_text()
        match = re.search(
            r"GIT_REPOSITORY\s+\S*oefp\.git\s+GIT_TAG\s+v(\d+\.\d+\.\d+)", text
        )
        assert match is not None, "no oefp FetchContent_Declare in CMakeLists.txt"
        assert match.group(1) == OEFP_VERSION

    def test_readme_requirements_line(self):
        """The README names the pinned version, all three components of it."""
        text = (REPO_ROOT / "README.md").read_text()
        match = re.search(r"\*\*OEFP\*\* (\d+\.\d+\.\d+)", text)
        assert match is not None, "no OEFP requirements line in README.md"
        assert match.group(1) == OEFP_VERSION


class TestVersionPins:
    """Each version location checked here agrees with ``oecluster.__version__``.

    Not every ``vrzn`` location is covered directly. ``oecluster.__version__``
    is the reference the rest are compared against, so a hand-edit there fails
    these tests from the other side rather than being asserted on its own, and
    the two ``tests/python/test_oecluster.py`` locations are asserted by that
    file's own ``TestVersion``.
    """

    def test_version_info_tuple(self):
        """The tuple form spells the same version as the string form."""
        expected = tuple(int(part) for part in oecluster.__version__.split("."))
        assert oecluster.__version_info__ == expected

    def test_root_pyproject_version(self):
        """The distribution metadata declares the package version."""
        assert _pyproject_version("pyproject.toml") == oecluster.__version__

    def test_python_pyproject_version(self):
        """The wheel metadata declares the package version.

        ``python/pyproject.toml`` is what a consumer resolves against, so a
        hand-edit here ships a wrong version without touching any other site.
        """
        assert _pyproject_version("python/pyproject.toml") == oecluster.__version__

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
        """The version literal in the SWIG interface matches the package version.

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


class TestManifestAgreement:
    """The two manifests declare the same ``project`` table, bar the known gaps.

    Both files describe the same distribution, and ``python/pyproject.toml`` is
    what a consumer resolves against, so a requirement added to one and not the
    other ships metadata that installs a different environment than the one CI
    tests. ``scikit-learn`` drifted exactly that way: it was declared in the
    root ``dev`` extra alone, leaving an install of the Python-only package's
    extra short a dependency that several test modules import. What is compared
    is the whole ``project`` table, so a field nobody thought to enumerate
    cannot drift quietly. The rule is asymmetric: the names in
    ``ROOT_ONLY_PROJECT_FIELDS`` and ``ROOT_ONLY_EXTRAS`` are subtracted from
    the root table only, never from the consumer's, so the consumer manifest
    acquiring one of them fails rather than being excused along with the root.
    """

    def test_dependencies_agree(self):
        """Both manifests declare the same runtime dependencies."""
        assert _pyproject_dependencies("python/pyproject.toml") == (
            _pyproject_dependencies("pyproject.toml")
        )

    def test_dev_extra_agrees(self):
        """Both manifests declare the same ``dev`` extra."""
        assert _pyproject_dev_extra("python/pyproject.toml") == (
            _pyproject_dev_extra("pyproject.toml")
        )

    def test_project_tables_agree(self):
        """The two ``project`` tables are identical under the asymmetric rule.

        The root-only names are subtracted from the root table only, so the
        comparison bites in both directions: the consumer manifest missing
        something the root declares, and the consumer manifest acquiring a
        field or an extra that only the root is permitted.

        This subsumes ``test_dependencies_agree`` and ``test_dev_extra_agrees``:
        both compare fields already inside this table. They are kept because
        they name the drift that actually shipped and report it more sharply
        than a whole-table diff does. This case earns its place by covering the
        fields nobody enumerated -- ``requires-python``, ``description``, an
        entry point, or an extra invented after this was written.
        """
        assert _pyproject_project_table(
            "python/pyproject.toml", drop_root_only=False
        ) == _pyproject_project_table("pyproject.toml", drop_root_only=True)
