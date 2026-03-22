#!/usr/bin/env python3

import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run(cmd, cwd=ROOT):
    completed = subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        target = Path(tmpdir) / "site"
        run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                ".",
                "--target",
                str(target),
                "--no-deps",
                "--no-build-isolation",
            ]
        )

        package_dir = target / "satzilla_features"
        if not (package_dir / "libsatzilla_features.so").exists():
            raise AssertionError("installed package is missing libsatzilla_features.so")
        if not (package_dir / "liblpk5.so").exists():
            raise AssertionError("installed package is missing liblpk5.so")
        if not (package_dir / "satzilla_Solvers" / "sbva").exists():
            raise AssertionError("installed package is missing bundled solver binaries")

        smoke_script = f"""
import sys
sys.path.insert(0, {str(target)!r})
from satzilla_features import extract_features

class CNFShim:
    def __init__(self, clauses):
        self.clauses = clauses
        self.nv = max(abs(l) for clause in clauses for l in clause)

features = extract_features(
    CNFShim([[1, 2], [-1, 3], [-2, 3, 4], [-3, -4, 5], [2, -5], [-1, -2, 4]]),
    groups=["structure", "ncnf-graphs"],
    build_if_missing=False,
)
assert "variable_alpha" in features
assert "vg_al_weights_mean" in features
"""
        run([sys.executable, "-c", smoke_script])

    print("package install smoke test passed")


if __name__ == "__main__":
    main()
