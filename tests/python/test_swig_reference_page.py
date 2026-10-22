"""docs/api/swig-classes.rst matches what the current build generates."""

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_the_swig_class_reference_is_current():
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "docs" / "swig_reference.py"), "--check"],
        capture_output=True, text=True, timeout=120, cwd=REPO_ROOT, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
