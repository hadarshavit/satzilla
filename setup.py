from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


ROOT = Path(__file__).resolve().parent
FEATURES_DIR = ROOT / "src"
PACKAGE_NAME = "satzilla_features"


class build_py(_build_py):
    def run(self):
        subprocess.run(["make", "shared"], cwd=FEATURES_DIR, check=True)
        super().run()

        package_dir = Path(self.build_lib) / PACKAGE_NAME
        package_dir.mkdir(parents=True, exist_ok=True)

        for artifact in (
            FEATURES_DIR / "libsatzilla_features.so",
            FEATURES_DIR / "lp_solve_5.0" / "liblpk5.so",
        ):
            shutil.copy2(artifact, package_dir / artifact.name)

        solver_src = FEATURES_DIR / "satzilla_Solvers"
        solver_dst = package_dir / "satzilla_Solvers"
        if solver_dst.exists():
            shutil.rmtree(solver_dst)
        shutil.copytree(solver_src, solver_dst)


setup(
    cmdclass={"build_py": build_py},
    package_data={
        "satzilla_features": [
            "py.typed",
            "*.so",
            "satzilla_Solvers/*",
            "satzilla_Solvers/yalsat-03v/*",
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
