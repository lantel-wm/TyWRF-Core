import json
import math
from pathlib import Path
import subprocess
import sys

import netCDF4
import numpy as np

from tools.audit_pressure_formula_column_probe import (
    audit_pressure_formula_column_probe,
    main as audit_main,
    report_to_json,
)


def _ensure_dim(dataset: netCDF4.Dataset, name: str, size: int | None) -> None:
    if name not in dataset.dimensions:
        dataset.createDimension(name, size)


def _variable_dims(name: str, data: np.ndarray) -> tuple[str, ...]:
    if data.ndim == 3:
        vertical = "bottom_top_stag" if name in {"PH", "PHB"} else "bottom_top"
        return ("Time", vertical, "south_north", "west_east")
    if data.ndim == 2:
        return ("Time", "south_north", "west_east")
    if data.ndim == 1:
        return ("Time", f"{name}_dim")
    return ("Time", *(f"{name}_dim_{index}" for index in range(data.ndim)))


def _write_wrfout(
    path: Path,
    variables: dict[str, np.ndarray],
    *,
    attrs: dict[str, object] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, "w") as dataset:
        _ensure_dim(dataset, "Time", None)
        for key, value in (attrs or {}).items():
            dataset.setncattr(key, value)

        for name, raw_values in variables.items():
            values = np.asarray(raw_values, dtype=np.float64)
            dims = _variable_dims(name, values)
            if values.ndim >= 2:
                _ensure_dim(dataset, "south_north", values.shape[-2])
                _ensure_dim(dataset, "west_east", values.shape[-1])
            if values.ndim == 3:
                _ensure_dim(dataset, dims[1], values.shape[-3])
            if values.ndim == 1:
                _ensure_dim(dataset, dims[1], values.shape[0])
            variable = dataset.createVariable(name, "f8", dims)
            variable[0, ...] = values


def _base_attrs() -> dict[str, object]:
    return {
        "TYWRF_DIAGNOSTIC_ONLY": "false",
        "TYWRF_GATE_CANDIDATE": "true",
        "TYWRF_INTEGRATOR_OUTPUT": "true",
        "TYWRF_VALIDATION_GATE_ONLY": "false",
        "TYWRF_CANDIDATE_KIND": "selected_field_integrator_v0",
        "TYWRF_FROM_PARENT_START": "0,0",
        "TYWRF_TO_PARENT_START": "1,0",
        "TYWRF_PARENT_GRID_RATIO": 1,
        "TYWRF_PRESSURE_REFRESH_APPLIED": "true",
        "TYWRF_PRESSURE_REFRESH_INTEGRATION_STATUS": "applied_to_candidate",
        "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 3,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 3,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_POINT_COUNT": 15,
        "TYWRF_PRESSURE_REFRESH_REFRESHED_P_POINTS": 15,
        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 15,
        "TYWRF_PROVIDER_TERRAIN_USES_MOVED_CANDIDATE_HGT": "true",
        "TYWRF_PRESSURE_REFRESH_TERRAIN_SOURCE": "moved_candidate_HGT",
    }


def _fields() -> dict[str, np.ndarray]:
    nz = 6
    ny = 3
    nx = 4
    return {
        "P": np.full((nz, ny, nx), 10_000.0),
        "PB": np.full((nz, ny, nx), 90_000.0),
        "MU": np.full((ny, nx), 900.0),
        "MUB": np.full((ny, nx), 10_000.0),
        "PH": np.full((nz + 1, ny, nx), 100.0),
        "PHB": np.full((nz + 1, ny, nx), 50_000.0),
        "T": np.full((nz, ny, nx), 5.0),
        "QVAPOR": np.full((nz, ny, nx), 0.01),
        "HGT": np.full((ny, nx), 50.0),
    }


def _write_source_tree(root: Path, *, complete_formula: bool = True) -> None:
    dynamics = root / "src" / "dynamics"
    include = root / "include" / "tywrf" / "dynamics"
    dynamics.mkdir(parents=True)
    include.mkdir(parents=True)
    if complete_formula:
        pressure_refresh = """
namespace tywrf::dynamics {
bool compute_krosa_pressure(const PressureRefreshInputs& inputs, const Options& options, int i, int j, int k, float& out) {
  const auto mu_total = inputs.mu(i, j) + inputs.mub(i, j);
  const auto pfu = inputs.c3f.values[k + 1] * mu_total + inputs.c4f.values[k + 1] + inputs.p_top_pa;
  const auto pfd = inputs.c3f.values[k] * mu_total + inputs.c4f.values[k] + inputs.p_top_pa;
  const auto phm = inputs.c3h.values[k] * mu_total + inputs.c4h.values[k] + inputs.p_top_pa;
  const auto phi = inputs.ph(i, j, k) + inputs.phb(i, j, k);
  const auto alb = inputs.alb(i, j, k);
  const auto theta = options.base_potential_temperature_k + inputs.t(i, j, k);
  const auto total_pressure = theta + phi + alb + pfu + pfd + phm;
  const auto pb = inputs.pb(i, j, k);
  out = total_pressure - pb;
  return true;
}
PressureRefreshReport refresh_krosa_moving_nest_pressure(const PressureRefreshInputs& inputs) {
  float refreshed_p = 0.0F;
  compute_krosa_pressure(inputs, refreshed_p);
  inputs.p(0, 0, 0) = refreshed_p;
  return {};
}
}
"""
    else:
        pressure_refresh = """
namespace tywrf::dynamics {
bool compute_krosa_pressure(const PressureRefreshInputs& inputs, int i, int j, int k, float& out) {
  out = inputs.p(i, j, k);
  return true;
}
}
"""
    (dynamics / "pressure_refresh.cpp").write_text(pressure_refresh, encoding="utf-8")
    (dynamics / "pressure_refresh_hook.cpp").write_text(
        """
namespace tywrf::dynamics {
Report apply_krosa_moving_nest_pressure_refresh_hook() {
  sync_exposed_3d(plan.mass, provider_views.pb, scratch_view.pb);
  sync_exposed_2d(plan.surface, provider_views.mub, scratch_view.mub);
  sync_exposed_3d(plan.w_full, provider_views.phb, scratch_view.phb);
  auto staging = make_krosa_pressure_refresh_inputs(scratch_view, provider_views.alb, metadata);
  return {};
}
}
""",
        encoding="utf-8",
    )
    (dynamics / "pressure_refresh_staging.cpp").write_text(
        """
namespace tywrf::dynamics {
Inputs make_krosa_pressure_refresh_inputs(StateView state, FieldView external_alb, Metadata metadata) {
  return {state.p, state.pb, state.t, external_alb, state.ph, state.phb, state.mu, state.mub, metadata.c3f, metadata.c4f, metadata.c3h, metadata.c4h, metadata.p_top_pa};
}
}
""",
        encoding="utf-8",
    )
    for name in (
        "pressure_refresh.hpp",
        "pressure_refresh_hook.hpp",
        "pressure_refresh_staging.hpp",
    ):
        (include / name).write_text("// test header\n", encoding="utf-8")


def _flag_codes(payload: dict[str, object]) -> set[str]:
    return {item["code"] for item in payload["risk_flags"]}


def test_column_probe_ranks_worst_target_columns_and_reports_low_levels(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference_end.nc"
    candidate = tmp_path / "candidate_end.nc"
    source = tmp_path / "source_start.nc"
    formula_json = tmp_path / "formula_inputs.json"
    source_root = tmp_path / "repo"

    ref_fields = _fields()
    cand_fields = {name: values.copy() for name, values in ref_fields.items()}
    source_fields = {name: values.copy() for name, values in ref_fields.items()}
    cand_fields["P"][:5, 2, 3] -= 2000.0
    cand_fields["P"][:5, 1, 3] -= 1000.0
    cand_fields["P"][:5, 0, 3] -= 100.0
    cand_fields["PB"][:5, :, 3] += 10.0
    cand_fields["MU"][:, 3] += 2.0
    cand_fields["MUB"][:, 3] += 1.0
    cand_fields["PHB"][:5, :, 3] += 1.0
    cand_fields["T"][:5, :, 3] += 0.1
    cand_fields["QVAPOR"][:5, :, 3] += 0.001

    _write_wrfout(reference, ref_fields)
    _write_wrfout(candidate, cand_fields, attrs=_base_attrs())
    _write_wrfout(source, source_fields)
    _write_source_tree(source_root)
    formula_json.write_text(
        json.dumps(
            {
                "status": "computed_with_flags",
                "formula_input_metrics": {
                    "target_region": {"status": "available"},
                    "pressure_refresh_metadata_counts": {
                        "TYWRF_PRESSURE_REFRESH_CHANGED_P_POINTS": 15
                    },
                },
                "risk_flags": [
                    {
                        "code": "low_level_target_p_bias",
                        "severity": "warning",
                        "message": "test risk",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    payload = json.loads(
        report_to_json(
            audit_pressure_formula_column_probe(
                candidate,
                reference_end_path=reference,
                source_start_path=source,
                formula_input_json_path=formula_json,
                source_root=source_root,
                max_columns=2,
            )
        )
    )

    assert payload["diagnostic_only"] is True
    assert payload["candidate_model_pass"] == "not_applicable"
    assert set(payload) >= {
        "column_selection",
        "column_formula_terms",
        "source_formula_inventory",
        "risk_flags",
        "summary",
    }

    selection = payload["column_selection"]
    assert selection["status"] == "ranked_by_reference"
    assert selection["target_region"]["inferred_target_column_count"] == 3
    selected = selection["selected_columns"]
    assert [(item["j"], item["i"]) for item in selected] == [(2, 3), (1, 3)]
    assert math.isclose(selected[0]["max_abs_low_level_p_error"], 2000.0)

    column = payload["column_formula_terms"]["columns"][0]
    level0 = column["levels"][0]
    p_record = level0["fields"]["P"]
    assert p_record["candidate"]["value"] == 8000.0
    assert p_record["reference"]["value"] == 10000.0
    assert p_record["source"]["value"] == 10000.0
    assert p_record["candidate_minus_reference"]["value"] == -2000.0
    assert p_record["candidate_minus_source"]["value"] == -2000.0

    full_pressure = level0["fields"]["P+PB"]
    assert full_pressure["candidate_minus_reference"]["value"] == -1990.0
    mu_total = level0["fields"]["MU+MUB"]
    assert mu_total["candidate_minus_reference"]["value"] == 3.0
    assert level0["fields"]["HGT"]["candidate"]["status"] == "available"
    assert level0["coefficient_terms"]["C3F/C4F"]["candidate"]["status"] == "not_available"
    assert level0["coefficient_terms"]["P_TOP"]["candidate"]["status"] == "not_available"
    assert level0["coefficient_terms"]["theta_m"]["candidate"]["status"] == "not_available"

    terms = payload["source_formula_inventory"]["compute_krosa_pressure"][
        "formula_terms"
    ]
    assert all(terms[term]["status"] == "present" for term in terms)

    codes = _flag_codes(payload)
    assert "low_level_column_p_error_large" in codes
    assert "p_bias_with_weak_column_companion_terms" in codes
    assert "column_candidate_source_delta_matches_p_error" in codes
    assert "coefficient_terms_missing_from_candidate_netcdf" in codes
    assert "formula_input_json_prior_risks_present" in codes
    assert payload["summary"]["strict_gate_status"] == "not_evaluated"
    json.dumps(payload, allow_nan=False)


def test_missing_reference_source_and_missing_source_terms_are_explicit(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "candidate.nc"
    source_root = tmp_path / "repo"
    _write_wrfout(
        candidate,
        {"P": np.full((2, 2, 4), 10_000.0)},
        attrs={
            **_base_attrs(),
            "TYWRF_PRESSURE_REFRESH_TARGET_COLUMN_COUNT": 2,
            "TYWRF_PRESSURE_REFRESH_REFRESHED_COLUMN_COUNT": 2,
        },
    )
    _write_source_tree(source_root, complete_formula=False)

    payload = json.loads(
        report_to_json(
            audit_pressure_formula_column_probe(
                candidate,
                source_root=source_root,
                max_columns=1,
            )
        )
    )

    assert payload["column_selection"]["status"] == "selected_without_reference"
    assert payload["column_selection"]["selected_columns"][0]["selection_basis"] == (
        "target_region_order_no_reference"
    )
    level0 = payload["column_formula_terms"]["columns"][0]["levels"][0]
    assert level0["fields"]["P"]["reference"]["status"] == "not_available"
    assert level0["fields"]["P"]["source"]["status"] == "not_available"
    assert level0["fields"]["PB"]["candidate"]["status"] == "not_available"
    assert level0["coefficient_terms"]["C3H/C4H"]["candidate"]["status"] == "not_available"

    compute = payload["source_formula_inventory"]["compute_krosa_pressure"]
    assert compute["missing_formula_terms"]

    codes = _flag_codes(payload)
    assert "reference_end_missing_ranked_by_target_order" in codes
    assert "source_start_missing" in codes
    assert "source_formula_terms_missing" in codes
    assert "coefficient_terms_missing_from_candidate_netcdf" in codes
    json.dumps(payload, allow_nan=False)


def test_cli_writes_json_and_direct_help_runs(tmp_path: Path) -> None:
    candidate = tmp_path / "candidate.nc"
    output = tmp_path / "column_probe.json"
    source_root = tmp_path / "repo"
    _write_wrfout(candidate, _fields(), attrs=_base_attrs())
    _write_source_tree(source_root)

    exit_code = audit_main(
        [
            "--candidate-end",
            str(candidate),
            "--source-root",
            str(source_root),
            "--max-columns",
            "1",
            "--output",
            str(output),
            "--pretty",
        ]
    )
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["diagnostic_only"] is True
    assert payload["column_formula_terms"]["columns"]

    result = subprocess.run(
        [sys.executable, "tools/audit_pressure_formula_column_probe.py", "--help"],
        cwd=Path(__file__).resolve().parents[2],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--candidate-end" in result.stdout
    assert "--formula-input-json" in result.stdout
