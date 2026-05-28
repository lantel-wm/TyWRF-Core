import json
import math
from pathlib import Path

import netCDF4

from tools.audit_pressure_column_runtime_probe import (
    audit_pressure_column_runtime_probe,
    main as audit_main,
    report_to_json,
)


FIELDS = "P,PB,P+PB,MU,MUB,MU+MUB,PH,PHB,PH+PHB,T,QVAPOR,HGT"
NOT_AVAILABLE = "ALB,C3F,C4F,C3H,C4H,P_TOP,theta_m"


def _write_candidate(path: Path, attrs: dict[str, object] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Time", None)
        for key, value in (attrs or {}).items():
            dataset.setncattr(key, value)


def _record(
    phase: str,
    i: int,
    j: int,
    k: int,
    *,
    p: float,
    pb: float = 99_700.0,
    mu: float = 4200.0,
    mub: float = 60_000.0,
) -> str:
    ph = 25_000.0 + k
    phb = 500.0 + k
    return (
        f"phase={phase};i={i};j={j};k={k};P={p};PB={pb};"
        f"P_PLUS_PB={p + pb};MU={mu};MUB={mub};MU_PLUS_MUB={mu + mub};"
        f"PH={ph};PHB={phb};PH_PLUS_PHB={ph + phb};T=1.5;"
        "QVAPOR=0.012;HGT=15"
    )


def _probe_attrs(*, values: str, record_count: int) -> dict[str, object]:
    return {
        "TYWRF_PRESSURE_COLUMN_PROBE_VERSION": "runtime_v0",
        "TYWRF_PRESSURE_COLUMN_PROBE_ENABLED": "true",
        "TYWRF_PRESSURE_COLUMN_PROBE_EVIDENCE_ONLY": "true",
        "TYWRF_PRESSURE_COLUMN_PROBE_INDEX_BASE": "zero_based_mass_grid",
        "TYWRF_PRESSURE_COLUMN_PROBE_COLUMN_COUNT": 1,
        "TYWRF_PRESSURE_COLUMN_PROBE_LEVEL_COUNT": 2,
        "TYWRF_PRESSURE_COLUMN_PROBE_PHASE_COUNT": 2,
        "TYWRF_PRESSURE_COLUMN_PROBE_RECORD_COUNT": record_count,
        "TYWRF_PRESSURE_COLUMN_PROBE_COLUMNS": "160,49",
        "TYWRF_PRESSURE_COLUMN_PROBE_LEVELS": "0,1",
        "TYWRF_PRESSURE_COLUMN_PROBE_PHASES": (
            "post_static_refresh,post_pressure_refresh"
        ),
        "TYWRF_PRESSURE_COLUMN_PROBE_FIELDS": FIELDS,
        "TYWRF_PRESSURE_COLUMN_PROBE_NOT_AVAILABLE": NOT_AVAILABLE,
        "TYWRF_PRESSURE_COLUMN_PROBE_VALUES": values,
    }


def _flag_codes(payload: dict[str, object]) -> set[str]:
    return {item["code"] for item in payload["risk_flags"]}


def test_runtime_probe_parses_two_phases_and_deltas(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    _write_candidate(candidate, _probe_attrs(values=values, record_count=4))

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert payload["usage_boundary"]["can_serve_as_gate_pass"] is False
    assert payload["usage_boundary"]["can_generate_candidate"] is False
    assert payload["columns"]["declared"] == [{"i": 160, "j": 49}]
    assert payload["levels"]["declared"] == [0, 1]
    assert payload["phases"]["declared"] == [
        "post_static_refresh",
        "post_pressure_refresh",
    ]
    assert payload["record_count"]["declared_matches_parsed"] is True
    assert payload["record_count"]["declared_matches_expected"] is True
    assert payload["record_count"]["expected_matches_parsed"] is True

    level0 = payload["observations"]["columns"][0]["levels"][0]
    assert level0["phases"]["post_static_refresh"]["P"] == 80.0
    assert level0["phases"]["post_pressure_refresh"]["P"] == 75.0
    delta = level0["deltas"]["post_pressure_refresh_minus_post_static_refresh"]
    assert math.isclose(delta["P"], -5.0)
    assert math.isclose(delta["P+PB"], -5.0)
    assert _flag_codes(payload) == {"formula_terms_unavailable"}


def test_runtime_probe_flags_mismatch_and_malformed_values(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            "phase=post_pressure_refresh;i=bad;j=49;k=0;P=79;broken",
        ]
    )
    _write_candidate(candidate, _probe_attrs(values=values, record_count=4))

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "malformed_values" in codes
    assert "record_count_mismatch" in codes
    assert payload["record_count"]["declared_matches_parsed"] is False
    assert payload["record_count"]["expected_matches_parsed"] is False
    assert payload["records"]["parsed"][0]["phase"] == "post_static_refresh"


def test_runtime_probe_flags_negative_post_refresh_and_large_drop(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.1953125),
            _record("post_static_refresh", 160, 49, 1, p=90.0),
            _record("post_pressure_refresh", 160, 49, 0, p=-4172.158691),
            _record("post_pressure_refresh", 160, 49, 1, p=70.0),
        ]
    )
    _write_candidate(candidate, _probe_attrs(values=values, record_count=4))

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "post_pressure_refresh_p_negative" in codes
    assert "large_p_drop_magnitude" in codes
    level0 = payload["observations"]["columns"][0]["levels"][0]
    delta = level0["deltas"]["post_pressure_refresh_minus_post_static_refresh"]
    assert math.isclose(delta["P"], -4252.3540035)


def test_runtime_probe_flags_missing_attrs(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    _write_candidate(candidate)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert payload["status"] == "computed_with_flags"
    assert "missing_probe_attrs" in codes
    assert "TYWRF_PRESSURE_COLUMN_PROBE_VALUES" in payload["attrs"]["missing_attrs"]
    assert payload["record_count"]["parsed"] == 0


def test_runtime_probe_cli_writes_pretty_json(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "runtime_probe.json"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    _write_candidate(candidate, _probe_attrs(values=values, record_count=4))

    exit_code = audit_main([str(candidate), "--output", str(output), "--pretty"])

    assert exit_code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["summary"]["parsed_record_count"] == 4
