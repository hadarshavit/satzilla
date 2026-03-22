#!/usr/bin/env python3

import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from satzilla_features import (
    FeatureExtractionError,
    extract_features,
    extract_features_from_path,
)


class CNFShim:
    def __init__(self, clauses):
        self.clauses = clauses
        self.nv = max(abs(lit) for clause in clauses for lit in clause)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    smoke_path = ROOT / "SAT-features-competition2024" / "tests" / "package_smoke.cnf"
    cnf = CNFShim(
        [
            [1, 2],
            [-1, 3],
            [-2, 3, 4],
            [-3, -4, 5],
            [2, -5],
            [-1, -2, 4],
        ]
    )

    features = extract_features(
        cnf,
        groups=["structure", "ncnf-graphs", "ncnf-constraints", "ncnf-rwh"],
        build_if_missing=False,
    )
    require("variable_alpha" in features, "missing structure feature")
    require("vg_al_weights_mean" in features, "missing new-cnf graph feature")
    require("and_node_mean" in features, "missing new-cnf constraint feature")
    require("rwh_2_mean" in features, "missing new-cnf rwh feature")
    require("solved" in features, "missing final solved feature")

    with tempfile.TemporaryDirectory() as tmpdir:
        cnf_path = Path(tmpdir) / "shim.cnf"
        with cnf_path.open("w", encoding="ascii") as handle:
            handle.write("p cnf 5 6\n")
            for clause in cnf.clauses:
                handle.write(" ".join(str(lit) for lit in clause) + " 0\n")

        path_features = extract_features_from_path(
            cnf_path,
            groups=["structure"],
            build_if_missing=False,
        )

    require(
        abs(features["variable_alpha"] - path_features["variable_alpha"]) < 1e-9,
        "python CNF and path extraction disagree on variable_alpha",
    )

    second_features = extract_features(
        cnf,
        groups=["structure"],
        build_if_missing=False,
    )
    require(
        abs(second_features["variable_alpha"] - path_features["variable_alpha"]) < 1e-9,
        "repeated direct calls are inconsistent",
    )

    smoke_features = extract_features_from_path(
        smoke_path,
        groups=["all"],
        build_if_missing=False,
    )
    require("nvarsOrig" in smoke_features, "all-groups direct path extraction failed")

    try:
        extract_features_from_path(
            ROOT / "does-not-exist.cnf",
            groups=["structure"],
            build_if_missing=False,
        )
    except FeatureExtractionError:
        pass
    else:
        raise AssertionError("missing file did not raise FeatureExtractionError")

    try:
        from pysat.formula import CNF  # type: ignore
    except ModuleNotFoundError:
        pass
    else:
        pysat_cnf = CNF(from_clauses=cnf.clauses)
        pysat_features = extract_features(
            pysat_cnf,
            groups=["structure"],
            build_if_missing=False,
        )
        require("variable_alpha" in pysat_features, "real PySAT CNF extraction failed")

    print("python interface smoke test passed")


if __name__ == "__main__":
    main()
