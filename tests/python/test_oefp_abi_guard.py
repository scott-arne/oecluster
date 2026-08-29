"""The OEFPBatch typemap must refuse a mismatched OEFP build.

oecluster compiles its own copy of OEFP and hands raw batch pointers to the
separately compiled ``oefp`` wheel, so a version disagreement is a memory
hazard rather than an ordinary error. These tests drive the guard in a
subprocess: it caches its verdict on first use, so the mismatch case needs a
fresh interpreter with the reported version patched before the first call.
"""

import subprocess
import sys

import numpy as np
import pytest

oefp = pytest.importorskip("oefp")


_MISMATCH_SCRIPT = """
import sys
import numpy as np
import oefp._native

# Move the reported build version off whatever oecluster compiled against.
oefp._native.OEFP_VERSION_PATCH = oefp._native.OEFP_VERSION_PATCH + 97

import oecluster
import oefp

# Create a minimal batch to trigger the typemap.
bits = np.array([[1, 0, 1, 0]], dtype=np.uint8)
fingerprints = []
for row in bits:
    on_bits = np.flatnonzero(row).astype(int).tolist()
    fingerprints.append(oefp.OEFP.from_on_bits(bits.shape[1], on_bits))
batch = oefp.OEFPBatch.from_fingerprints(fingerprints)

# This should fail with the ABI guard message.
oecluster.bitbirch(batch, threshold=0.5, branching_factor=2)
"""


def test_matching_build_is_accepted():
    """The installed oefp matches this build, so a batch call must work."""
    import oecluster

    # Create a minimal batch and call bitbirch to exercise the typemap.
    bits = np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=np.uint8)
    fingerprints = []
    for row in bits:
        on_bits = np.flatnonzero(row).astype(int).tolist()
        fingerprints.append(oefp.OEFP.from_on_bits(bits.shape[1], on_bits))
    batch = oefp.OEFPBatch.from_fingerprints(fingerprints)

    result = oecluster.bitbirch(batch, threshold=0.5, branching_factor=2)
    assert isinstance(result, oecluster.BitBirchResult)


def test_mismatched_build_is_rejected(tmp_path):
    script = tmp_path / "mismatch.py"
    script.write_text(_MISMATCH_SCRIPT)
    result = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, check=False
    )
    assert result.returncode != 0
    assert "compiled against OEFP" in result.stderr
    assert "raw fingerprint-batch pointers" in result.stderr
