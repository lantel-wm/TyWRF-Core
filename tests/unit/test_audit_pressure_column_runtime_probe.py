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
FORMULA_FIELDS = (
    "status,valid,i,j,k,mu_total,pfu,pfd,phm,log_ratio,phi_lower,phi_upper,"
    "delta_phi,ALB,PB,theta,alpha_total,alpha_perturbation,"
    "alpha_from_wrf_branch,pressure_base,total_pressure,perturbation_pressure_pa"
)


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


def _formula_record(
    *,
    i: int = 160,
    j: int = 49,
    k: int = 0,
    status: str = "recorded",
    valid: bool = True,
    perturbation_pressure_pa: float = 75.0,
    pb: float = 99_700.0,
    total_pressure: float | None = None,
) -> str:
    total_pressure = (
        pb + perturbation_pressure_pa
        if total_pressure is None
        else total_pressure
    )
    return (
        f"status={status};valid={str(valid).lower()};i={i};j={j};k={k};"
        "mu_total=64200;pfu=0.72;pfd=0.28;phm=25500;log_ratio=0.003;"
        f"phi_lower=25100;phi_upper=25230;delta_phi=130;ALB=0.92;PB={pb};"
        "theta=301.5;alpha_total=0.83;alpha_perturbation=0.014;"
        f"alpha_from_wrf_branch=0.014;pressure_base={pb};"
        f"total_pressure={total_pressure};"
        f"perturbation_pressure_pa={perturbation_pressure_pa}"
    )


def _formula_attrs(
    *,
    values: str,
    record_count: int,
    fields: str = FORMULA_FIELDS,
    valid_count: int = 1,
    invalid_count: int = 0,
    out_of_bounds_count: int = 0,
    outside_target_region_count: int = 0,
) -> dict[str, object]:
    return {
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_VERSION": "runtime_v0",
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_ENABLED": "true",
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_EVIDENCE_ONLY": "true",
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_INDEX_BASE": "zero_based_mass_grid",
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_REQUEST_COUNT": record_count,
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_RECORD_COUNT": record_count,
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALID_COUNT": valid_count,
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_INVALID_COUNT": invalid_count,
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUT_OF_BOUNDS_COUNT": (
            out_of_bounds_count
        ),
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_OUTSIDE_TARGET_REGION_COUNT": (
            outside_target_region_count
        ),
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_FIELDS": fields,
        "TYWRF_PRESSURE_FORMULA_OBSERVATION_VALUES": values,
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
    assert payload["formula_observation"]["present"] is False
    assert "pressure_budget" not in payload["formula_observation"]["correlation"]
    assert (
        "formula_sensitivity"
        not in payload["formula_observation"]["correlation"]
    )
    assert _flag_codes(payload) == {"formula_terms_unavailable"}


def test_formula_observation_correlates_with_probe_pressure(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(
        _formula_attrs(
            values=_formula_record(perturbation_pressure_pa=75.0),
            record_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_terms_unavailable" not in codes
    assert "probe_field_terms_unavailable" in codes
    assert "formula_pressure_probe_mismatch" not in codes
    assert "formula_observation_count_mismatch" not in codes
    assert "formula_observation_invalid_status" not in codes
    formula = payload["formula_observation"]
    assert formula["present"] is True
    assert formula["attrs"]["version"] == "runtime_v0"
    assert formula["record_count"]["valid_declared"] == 1
    assert formula["record_count"]["valid_observed"] == 1
    assert formula["record_count"]["valid_matches_declared"] is True
    assert {"ALB", "theta", "perturbation_pressure_pa"} <= set(
        formula["fields"]["supplied"]
    )
    assert "160,49,0" in formula["records_by_column_level"]["items"]
    match = formula["correlation"]["matched_records"][0]
    assert math.isclose(match["difference_pa"], 0.0)
    assert math.isclose(match["probe_delta_p"], -5.0)
    budget = formula["correlation"]["pressure_budget"]["records"][0]
    assert math.isclose(budget["total_pressure_minus_pb_pa"], 75.0)
    assert math.isclose(budget["perturbation_minus_total_minus_pb_pa"], 0.0)
    assert budget["formula_total_pressure_below_pb"] is False
    assert budget["post_refresh_p_matches_total_minus_pb"] is True
    assert math.isclose(budget["probe_delta_p"], -5.0)
    assert budget["large_drop_explained_by_formula_base_subtraction"] is False
    sensitivity = formula["correlation"]["pressure_budget"]["formula_sensitivity"]
    sensitivity_record = sensitivity["records"][0]
    assert math.isclose(sensitivity_record["total_pressure_gap_to_pb_pa"], -75.0)
    assert math.isclose(
        sensitivity_record[
            "pressure_change_needed_to_make_perturbation_nonnegative_pa"
        ],
        0.0,
    )
    assert math.isclose(
        sensitivity_record[
            "approximate_fractional_total_pressure_increase_needed"
        ],
        0.0,
    )
    assert sensitivity_record["available_drivers"]["theta"] == 301.5
    assert sensitivity["summary"][
        "records_requiring_total_pressure_increase_count"
    ] == 0


def test_formula_observation_flags_count_mismatch_and_malformed_values(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    malformed = _formula_record(perturbation_pressure_pa=75.0).replace(
        "i=160", "i=bad"
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(_formula_attrs(values=malformed, record_count=2, valid_count=2))
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "malformed_formula_observation_values" in codes
    assert "formula_observation_count_mismatch" in codes
    assert payload["formula_observation"]["record_count"]["parsed"] == 0
    assert payload["formula_observation"]["record_count"]["declared"] == 2


def test_formula_observation_flags_pressure_mismatch(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(
        _formula_attrs(
            values=_formula_record(perturbation_pressure_pa=75.25),
            record_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_pressure_probe_mismatch" in codes
    mismatch = payload["formula_observation"]["correlation"]["pressure_mismatches"][0]
    assert math.isclose(mismatch["difference_pa"], 0.25)
    budget = payload["formula_observation"]["correlation"]["pressure_budget"][
        "records"
    ][0]
    assert budget["post_refresh_p_matches_total_minus_pb"] is False
    assert math.isclose(budget["total_pressure_minus_pb_pa"], 75.25)


def test_formula_observation_pressure_budget_flags_below_base_pressure(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    pb = 99_711.8828
    total_pressure = 95_539.724
    perturbation_pressure_pa = total_pressure - pb
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.1953125, pb=pb),
            _record(
                "post_pressure_refresh",
                160,
                49,
                0,
                p=perturbation_pressure_pa,
                pb=pb,
            ),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=2)
    attrs["TYWRF_PRESSURE_COLUMN_PROBE_LEVEL_COUNT"] = 1
    attrs["TYWRF_PRESSURE_COLUMN_PROBE_LEVELS"] = "0"
    attrs.update(
        _formula_attrs(
            values=_formula_record(
                perturbation_pressure_pa=perturbation_pressure_pa,
                pb=pb,
                total_pressure=total_pressure,
            ),
            record_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_total_pressure_below_base_pressure" in codes
    assert "formula_pressure_drop_explained_by_base_subtraction" in codes
    assert "formula_total_pressure_gap_requires_staging_diagnosis" in codes
    budget = payload["formula_observation"]["correlation"]["pressure_budget"][
        "records"
    ][0]
    assert math.isclose(budget["total_pressure_minus_pb_pa"], perturbation_pressure_pa)
    assert math.isclose(budget["perturbation_minus_total_minus_pb_pa"], 0.0)
    assert budget["formula_total_pressure_below_pb"] is True
    assert budget["post_refresh_p_matches_total_minus_pb"] is True
    assert math.isclose(budget["probe_delta_p"], -4252.3541125)
    assert budget["large_drop_explained_by_formula_base_subtraction"] is True
    sensitivity = payload["formula_observation"]["correlation"]["pressure_budget"][
        "formula_sensitivity"
    ]
    sensitivity_record = sensitivity["records"][0]
    assert math.isclose(
        sensitivity_record["total_pressure_gap_to_pb_pa"],
        pb - total_pressure,
    )
    assert math.isclose(
        sensitivity_record[
            "pressure_change_needed_to_make_perturbation_nonnegative_pa"
        ],
        pb - total_pressure,
    )
    assert math.isclose(
        sensitivity_record[
            "approximate_fractional_total_pressure_increase_needed"
        ],
        (pb - total_pressure) / total_pressure,
    )
    assert sensitivity_record["exact_budget_directional_sensitivities"][
        "perturbation_pressure_pa"
    ]["with_respect_to_PB"] == -1.0
    assert sensitivity["summary"][
        "most_records_require_large_total_pressure_increase"
    ] is True
    assert sensitivity["summary"][
        "records_requiring_large_total_pressure_increase_count"
    ] == 1


def test_formula_sensitivity_preserves_records_missing_budget_terms(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
        ]
    )
    formula_values = (
        "status=recorded;valid=true;i=160;j=49;k=0;"
        "perturbation_pressure_pa=75"
    )
    formula_fields = "status,valid,i,j,k,perturbation_pressure_pa"
    attrs = _probe_attrs(values=values, record_count=2)
    attrs["TYWRF_PRESSURE_COLUMN_PROBE_LEVEL_COUNT"] = 1
    attrs["TYWRF_PRESSURE_COLUMN_PROBE_LEVELS"] = "0"
    attrs.update(
        _formula_attrs(
            values=formula_values,
            record_count=1,
            fields=formula_fields,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_observation_nonfinite_or_invalid_terms" not in codes
    assert "formula_total_pressure_gap_requires_staging_diagnosis" not in codes
    assert "formula_pressure_probe_mismatch" not in codes
    correlation = payload["formula_observation"]["correlation"]
    assert "pressure_budget" not in correlation
    assert correlation["matched_records"][0]["matches_within_tolerance"] is True


def test_formula_observation_flags_nonfinite_terms(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(
        _formula_attrs(
            values=_formula_record().replace("ALB=0.92", "ALB=nan"),
            record_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_observation_nonfinite_or_invalid_terms" in codes
    nonfinite = next(
        flag
        for flag in payload["risk_flags"]
        if flag["code"] == "formula_observation_nonfinite_or_invalid_terms"
    )
    assert nonfinite["evidence"]["errors"][0]["field"] == "ALB"


def test_formula_observation_accepts_cpp_skip_status_counts(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    formula_values = "|".join(
        [
            _formula_record(perturbation_pressure_pa=75.0),
            _formula_record(
                i=0,
                j=0,
                k=0,
                status="request_out_of_bounds",
                valid=False,
            ),
            _formula_record(
                i=1,
                j=1,
                k=0,
                status="request_outside_target_region",
                valid=False,
            ),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(
        _formula_attrs(
            values=formula_values,
            record_count=3,
            valid_count=1,
            invalid_count=0,
            out_of_bounds_count=1,
            outside_target_region_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_observation_count_mismatch" not in codes
    assert "formula_observation_invalid_status" not in codes
    record_count = payload["formula_observation"]["record_count"]
    assert record_count["valid_observed"] == 1
    assert record_count["valid_matches_declared"] is True
    assert record_count["invalid_observed"] == 0
    assert record_count["invalid_matches_declared"] is True
    assert record_count["out_of_bounds_observed"] == 1
    assert record_count["out_of_bounds_matches_declared"] is True
    assert record_count["outside_target_region_observed"] == 1
    assert record_count["outside_target_region_matches_declared"] is True


def test_formula_observation_flags_invalid_status(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs.update(
        _formula_attrs(
            values=_formula_record(status="invalid_mu_total", valid=False),
            record_count=1,
            valid_count=0,
            invalid_count=1,
        )
    )
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "formula_observation_invalid_status" in codes
    invalid = next(
        flag
        for flag in payload["risk_flags"]
        if flag["code"] == "formula_observation_invalid_status"
    )
    assert invalid["evidence"]["examples"][0]["status"] == "invalid_mu_total"


def test_formula_observation_flags_missing_when_enabled_claimed(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    values = "|".join(
        [
            _record("post_static_refresh", 160, 49, 0, p=80.0),
            _record("post_static_refresh", 160, 49, 1, p=120.0),
            _record("post_pressure_refresh", 160, 49, 0, p=75.0),
            _record("post_pressure_refresh", 160, 49, 1, p=118.0),
        ]
    )
    attrs = _probe_attrs(values=values, record_count=4)
    attrs["TYWRF_PRESSURE_COLUMN_PROBE_FORMULA_OBSERVATION_ENABLED"] = "true"
    _write_candidate(candidate, attrs)

    payload = json.loads(report_to_json(audit_pressure_column_runtime_probe(candidate)))
    codes = _flag_codes(payload)

    assert "missing_formula_observation" in codes
    assert payload["formula_observation"]["present"] is False


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
