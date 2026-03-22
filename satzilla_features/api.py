from __future__ import annotations

import ctypes
import subprocess
from pathlib import Path
from typing import Iterable, Mapping, Sequence


_PACKAGE_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _PACKAGE_DIR.parent
_FEATURES_DIR = _REPO_ROOT / "SAT-features-competition2024"
_SOURCE_SHARED_LIBRARY = _FEATURES_DIR / "libsatzilla_features.so"
_PACKAGE_SHARED_LIBRARY = _PACKAGE_DIR / "libsatzilla_features.so"

ALL_GROUPS = (
    "all",
    "base",
    "structure",
    "ncnf-graphs",
    "ncnf-constraints",
    "ncnf-rwh",
    "sp",
    "dia",
    "cl",
    "lp",
    "unit",
    "ls",
    "lobjois",
)

DEFAULT_GROUPS: tuple[str, ...] = ()


class FeatureExtractionError(RuntimeError):
    pass


class _COptions(ctypes.Structure):
    _fields_ = [
        ("do_base", ctypes.c_int),
        ("do_structure", ctypes.c_int),
        ("do_ncnf_graphs", ctypes.c_int),
        ("do_ncnf_constraints", ctypes.c_int),
        ("do_ncnf_rwh", ctypes.c_int),
        ("do_unit_probe", ctypes.c_int),
        ("do_ls_probe", ctypes.c_int),
        ("do_cl", ctypes.c_int),
        ("do_dia", ctypes.c_int),
        ("do_sp", ctypes.c_int),
        ("do_lobjois", ctypes.c_int),
        ("do_lp", ctypes.c_int),
        ("timeout_seconds", ctypes.c_int),
        ("solver_root", ctypes.c_char_p),
    ]


class _CResult(ctypes.Structure):
    _fields_ = [
        ("status_code", ctypes.c_int),
        ("feature_count", ctypes.c_size_t),
        ("feature_names", ctypes.POINTER(ctypes.c_char_p)),
        ("feature_values", ctypes.POINTER(ctypes.c_double)),
        ("error_message", ctypes.c_char_p),
    ]


_GROUP_ATTRS = {
    "base": "do_base",
    "structure": "do_structure",
    "ncnf-graphs": "do_ncnf_graphs",
    "ncnf-constraints": "do_ncnf_constraints",
    "ncnf-rwh": "do_ncnf_rwh",
    "unit": "do_unit_probe",
    "ls": "do_ls_probe",
    "cl": "do_cl",
    "dia": "do_dia",
    "sp": "do_sp",
    "lobjois": "do_lobjois",
    "lp": "do_lp",
}

_LIB = None


def _normalize_groups(groups: Iterable[str] | None) -> list[str]:
    if groups is None:
        return []

    normalized = []
    for group in groups:
        key = group.strip().lower()
        if key not in ALL_GROUPS:
            raise ValueError(
                f"unknown feature group {group!r}; expected one of {', '.join(ALL_GROUPS)}"
            )
        normalized.append(key)
    return normalized


def _candidate_library_paths() -> list[Path]:
    return [_PACKAGE_SHARED_LIBRARY, _SOURCE_SHARED_LIBRARY]


def _find_default_solver_root(lib_path: Path) -> Path:
    package_candidate = lib_path.parent
    if (package_candidate / "satzilla_Solvers").exists():
        return package_candidate
    return _FEATURES_DIR


def _ensure_shared_library(build_if_missing: bool) -> Path:
    for candidate in _candidate_library_paths():
        if candidate.exists():
            return candidate

    if not build_if_missing:
        raise FileNotFoundError(
            "shared feature extractor library not found. "
            "Build it with `make shared` in SAT-features-competition2024."
        )

    completed = subprocess.run(
        ["make", "shared"],
        cwd=_FEATURES_DIR,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise FeatureExtractionError(
            "failed to build shared feature extractor with `make shared`\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    for candidate in _candidate_library_paths():
        if candidate.exists():
            return candidate

    raise FileNotFoundError("`make shared` completed but no shared library was found")


def _load_library(build_if_missing: bool):
    global _LIB
    if _LIB is not None:
        return _LIB

    library_path = _ensure_shared_library(build_if_missing)
    lib = ctypes.CDLL(str(library_path))
    lib.satzilla_extract_from_dimacs_path.argtypes = [
        ctypes.c_char_p,
        ctypes.POINTER(_COptions),
        ctypes.POINTER(_CResult),
    ]
    lib.satzilla_extract_from_dimacs_path.restype = ctypes.c_int
    lib.satzilla_extract_from_clauses.argtypes = [
        ctypes.c_int,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(_COptions),
        ctypes.POINTER(_CResult),
    ]
    lib.satzilla_extract_from_clauses.restype = ctypes.c_int
    lib.satzilla_free_result.argtypes = [ctypes.POINTER(_CResult)]
    lib.satzilla_free_result.restype = None
    _LIB = lib
    return lib


def _extract_cnf_data(cnf: object) -> tuple[Sequence[Sequence[int]], int]:
    clauses = getattr(cnf, "clauses", None)
    nv = getattr(cnf, "nv", None)
    if clauses is None or nv is None:
        raise TypeError(
            "expected a PySAT CNF-like object exposing `clauses` and `nv`"
        )
    return clauses, int(nv)


def _build_options(
    groups: Sequence[str],
    timeout: int | None,
    solver_root: str | Path | None,
    build_if_missing: bool,
) -> _COptions:
    lib_path = _ensure_shared_library(build_if_missing)
    resolved_solver_root = (
        Path(solver_root) if solver_root is not None else _find_default_solver_root(lib_path)
    )

    options = _COptions()
    for group in groups:
        if group == "all":
            for attr in _GROUP_ATTRS.values():
                setattr(options, attr, 1)
        else:
            setattr(options, _GROUP_ATTRS[group], 1)

    options.timeout_seconds = -1 if timeout is None else int(timeout)
    options.solver_root = str(resolved_solver_root).encode("utf-8")
    return options


def _decode_result(result: _CResult) -> Mapping[str, float]:
    features: dict[str, float] = {}
    for idx in range(result.feature_count):
        name = result.feature_names[idx].decode("utf-8")
        features[name] = result.feature_values[idx]
    return features


def _raise_if_failed(status: int, result: _CResult) -> None:
    if status == 0 and result.status_code == 0:
        return

    if result.error_message:
        message = ctypes.cast(result.error_message, ctypes.c_char_p).value.decode("utf-8")
    else:
        message = f"feature extraction failed with status {status or result.status_code}"
    raise FeatureExtractionError(message)


def extract_features(
    cnf: object,
    *,
    groups: Iterable[str] | None = None,
    timeout: int | None = None,
    solver_root: str | Path | None = None,
    build_if_missing: bool = True,
) -> Mapping[str, float]:
    """
    Extract SAT features for a PySAT CNF-like object.

    PySAT documents `CNF` formulas in the `formula` module and exposes CNF data
    through the formula object interface used here (`clauses` and `nv`).
    Source: https://pysathq.github.io/docs/html/api/formula.html
    """

    clauses, nv = _extract_cnf_data(cnf)
    normalized_groups = _normalize_groups(groups)
    options = _build_options(normalized_groups, timeout, solver_root, build_if_missing)
    lib = _load_library(build_if_missing)

    flat_literals: list[int] = []
    clause_offsets = [0]
    for clause in clauses:
        flat_literals.extend(int(lit) for lit in clause)
        clause_offsets.append(len(flat_literals))

    literal_array = (ctypes.c_int * len(flat_literals))(*flat_literals) if flat_literals else None
    offset_array = (ctypes.c_size_t * len(clause_offsets))(*clause_offsets)
    result = _CResult()
    try:
        status = lib.satzilla_extract_from_clauses(
            int(nv),
            len(clauses),
            literal_array,
            offset_array,
            ctypes.byref(options),
            ctypes.byref(result),
        )
        _raise_if_failed(status, result)
        return _decode_result(result)
    finally:
        lib.satzilla_free_result(ctypes.byref(result))


def extract_features_from_path(
    cnf_path: str | Path,
    *,
    groups: Iterable[str] | None = None,
    timeout: int | None = None,
    solver_root: str | Path | None = None,
    build_if_missing: bool = True,
) -> Mapping[str, float]:
    normalized_groups = _normalize_groups(groups)
    options = _build_options(normalized_groups, timeout, solver_root, build_if_missing)
    lib = _load_library(build_if_missing)

    result = _CResult()
    try:
        status = lib.satzilla_extract_from_dimacs_path(
            str(cnf_path).encode("utf-8"),
            ctypes.byref(options),
            ctypes.byref(result),
        )
        _raise_if_failed(status, result)
        return _decode_result(result)
    finally:
        lib.satzilla_free_result(ctypes.byref(result))
