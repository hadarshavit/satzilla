#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BINARY = ROOT / "features"
SMOKE_CNF = ROOT / "tests" / "package_smoke.cnf"


def run(cmd):
    completed = subprocess.run(
        cmd,
        cwd=ROOT,
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


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "features.csv"
        run([str(BINARY), "-all", str(SMOKE_CNF), str(output_path)])

        with output_path.open(newline="") as handle:
            rows = list(csv.reader(handle))

    require(len(rows) == 2, f"expected 2 csv rows, got {len(rows)}")
    header, values = rows
    require(len(header) == len(values), "header/value length mismatch")

    expected_features = {
        "nvarsOrig",
        "nclausesOrig",
        "vig_modularty",
        "vig_d_poly",
        "cvig_db_poly",
        "variable_alpha",
        "v_nd_p_node_mean",
        "vg_al_weights_mean",
        "cg_al_node_mean",
        "rg_weights_mean",
        "big_node_mean",
        "and_node_mean",
        "band_node_mean",
        "exo_node_mean",
        "rwh_0_mean",
        "rwh_1_mean",
        "rwh_2_mean",
        "structure-featuretime",
        "ncnf-graphs-featuretime",
        "ncnf-constraints-featuretime",
        "ncnf-rwh-featuretime",
        "solved",
    }

    header_set = set(header)
    missing = sorted(expected_features - header_set)
    require(not missing, f"missing expected features: {missing}")

    feature_map = dict(zip(header, values))
    require(feature_map["solved"] in {"0.000000000", "1.000000000", "2.000000000", "3.000000000", "4.000000000"}, "unexpected solved value")

    for feature_name in (
        "vig_modularty",
        "variable_alpha",
        "v_nd_p_node_mean",
        "and_node_mean",
        "rwh_0_mean",
    ):
        value = float(feature_map[feature_name])
        require(value > -512.0, f"{feature_name} was left at the reserved value")

    print("package smoke test passed")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
