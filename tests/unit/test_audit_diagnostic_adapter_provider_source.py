import json
from pathlib import Path

import netCDF4

from tools.audit_diagnostic_adapter_provider_source import (
    audit_diagnostic_adapter_provider_source,
    main as audit_main,
    report_to_json,
)


PROVIDER = "TYWRF_DIAGNOSTIC_ADAPTER_PROVIDER_SOURCE_"
STAGING = "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_STAGING_"
DELTA = "TYWRF_DIAGNOSTIC_ADAPTER_SOURCE_CHILD_DELTA_"
DELTA_FIELDS = ("PHB", "MUB", "HT", "PB", "T_INIT", "ALB")


def _bool(value: bool) -> str:
    return "true" if value else "false"


def _provider_attrs() -> dict[str, object]:
    return {
        f"{PROVIDER}VERSION": "d77_provider_source_v0",
        f"{PROVIDER}ORIGIN": "base_state_reconstruction_provider+moved_candidate_HGT",
        f"{PROVIDER}SOURCE_ORIGIN": (
            "base_state_reconstruction_provider+moved_candidate_HGT"
        ),
        f"{PROVIDER}PROVIDER_SOURCE": "base_state_reconstruction_provider",
        f"{PROVIDER}TERRAIN_SOURCE": "moved_candidate_HGT",
        f"{PROVIDER}TERRAIN_PROVENANCE": "override:moved_candidate_HGT",
        f"{PROVIDER}TERRAIN_OVERRIDE_USED": "true",
        f"{PROVIDER}HT_SOURCE": "output_static.hgt",
        f"{PROVIDER}HT_PROVENANCE": "adapter_HT_from_output_static_HGT",
        f"{PROVIDER}PROVIDER_OK": "true",
        f"{PROVIDER}RESULT_MESSAGE": "",
        f"{PROVIDER}DIAGNOSTIC_ONLY": "true",
        f"{PROVIDER}GATE_CANDIDATE": "false",
        f"{PROVIDER}INTEGRATOR_OUTPUT": "false",
        f"{PROVIDER}WRITES_CANDIDATE": "false",
        f"{PROVIDER}WRITES_NETCDF": "false",
        f"{PROVIDER}NO_CANDIDATE_WRITE": "true",
        f"{PROVIDER}USES_REFERENCE_END_TRUTH": "false",
        f"{PROVIDER}NO_REFERENCE_END_TRUTH": "true",
        f"{PROVIDER}USES_DIRECT_P_SHORTCUT": "false",
        f"{PROVIDER}NO_DIRECT_P_SHORTCUT": "true",
        f"{PROVIDER}READS_DIRECT_P": "false",
        f"{PROVIDER}WROTE_PB": "true",
        f"{PROVIDER}WROTE_T_INIT": "true",
        f"{PROVIDER}WROTE_MUB": "true",
        f"{PROVIDER}WROTE_ALB": "true",
        f"{PROVIDER}WROTE_PHB": "true",
        f"{PROVIDER}PROVIDER_RECONSTRUCTED_PHB_NOT_WRF_REBALANCE_VALIDATED": (
            "true"
        ),
    }


def _staging_attrs() -> dict[str, object]:
    return {
        f"{STAGING}VERSION": "d75_provider_report_v0",
        f"{STAGING}PROVIDER_KIND": "BaseStateSourceStagingProvider",
        f"{STAGING}SOURCE": "explicit_base_state_source_staging_provider",
        f"{STAGING}DISPOSITION": "staging_report_only_no_gate_no_integrator",
        f"{STAGING}SOURCE_SHAPE": "active_nx=4 active_ny=3 mass_nz=2 full_nz=3",
        f"{STAGING}OK": "true",
        f"{STAGING}RESULT_MESSAGE": "",
        f"{STAGING}DIAGNOSTIC_ONLY": "true",
        f"{STAGING}GATE_CANDIDATE": "false",
        f"{STAGING}INTEGRATOR_OUTPUT": "false",
        f"{STAGING}WRITES_CANDIDATE": "false",
        f"{STAGING}WRITES_NETCDF": "false",
        f"{STAGING}CANDIDATE_BUFFERS_PRESERVED": "true",
        f"{STAGING}OWNS_STAGING_BUFFERS": "true",
        f"{STAGING}ALLOCATED_BUFFERS": "true",
        f"{STAGING}USES_REFERENCE_END_TRUTH": "false",
        f"{STAGING}USES_DIRECT_P_SHORTCUT": "false",
        f"{STAGING}READS_DIRECT_P": "false",
        f"{STAGING}ALIASES_CHILD": "false",
        f"{STAGING}EXPOSED_REGION_COUNT": 1.0,
        f"{STAGING}STAGED_VALUE_COUNT": 12.0,
        f"{STAGING}INVALID_EXPOSED_VALUE_COUNT": 0.0,
    }


def _delta_attrs(*, max_by_field: dict[str, float] | None = None) -> dict[str, object]:
    max_by_field = max_by_field or {}
    attrs: dict[str, object] = {
        f"{DELTA}VERSION": "a76_source_child_delta_v0",
        f"{DELTA}SOURCE": "BaseStateSourceStagingProvider_vs_child_staging_pre_adapter",
        f"{DELTA}SCOPE": "exposed_base_state_values_only",
        f"{DELTA}FIELDS": "PHB,MUB,HT,PB,T_INIT,ALB",
        f"{DELTA}DIAGNOSTIC_ONLY": "true",
        f"{DELTA}GATE_CANDIDATE": "false",
        f"{DELTA}INTEGRATOR_OUTPUT": "false",
        f"{DELTA}WRITES_CANDIDATE": "false",
        f"{DELTA}WRITES_NETCDF": "false",
        f"{DELTA}VALUES_IDENTICAL": "false",
        f"{DELTA}COMPARED_VALUE_COUNT": 12.0,
        f"{DELTA}DIFFERING_VALUE_COUNT": 6.0,
        f"{DELTA}MAX_ABS_DIFF": max(max_by_field.values(), default=6.0),
    }
    for index, field in enumerate(DELTA_FIELDS, start=1):
        max_abs_diff = max_by_field.get(field, float(index))
        attrs[f"{DELTA}{field}_COMPARED_VALUE_COUNT"] = 2.0
        attrs[f"{DELTA}{field}_DIFFERING_VALUE_COUNT"] = 1.0
        attrs[f"{DELTA}{field}_MAX_ABS_DIFF"] = max_abs_diff
    return attrs


def _write_netcdf(path: Path, attrs: dict[str, object]) -> None:
    with netCDF4.Dataset(path, "w") as dataset:
        for key, value in attrs.items():
            dataset.setncattr(key, value)


def _valid_attrs() -> dict[str, object]:
    return _provider_attrs() | _staging_attrs() | _delta_attrs()


def test_valid_diagnostic_metadata_passes_without_gate_evidence(tmp_path: Path) -> None:
    candidate = tmp_path / "diagnostic_adapter.nc"
    output = tmp_path / "audit.json"
    _write_netcdf(candidate, _valid_attrs())

    exit_code = audit_main([str(candidate), "--output", str(output), "--pretty"])
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert exit_code == 0
    assert payload["status"] == "passed"
    assert payload["gate_evidence"] is False
    assert payload["provider_source"]["source_origin"] == (
        "base_state_reconstruction_provider+moved_candidate_HGT"
    )
    assert payload["provider_source"]["terrain_source"] == "moved_candidate_HGT"
    assert payload["provider_source"]["provider_ok"] is True
    assert payload["provider_source"]["diagnostic_only"] is True
    assert payload["provider_source"]["gate_candidate"] is False
    assert payload["provider_source"]["integrator_output"] is False
    assert payload["provider_source"]["writes_candidate"] is False
    assert payload["provider_source"]["writes_netcdf"] is False
    assert payload["provider_source"]["no_candidate_write"] is True
    assert payload["provider_source"]["uses_reference_end_truth"] is False
    assert payload["provider_source"]["uses_direct_p_shortcut"] is False
    assert payload["provider_source"]["reads_direct_p"] is False
    assert payload["provider_source"]["wrote_pb"] is True
    assert payload["provider_source"]["wrote_t_init"] is True
    assert payload["provider_source"]["wrote_mub"] is True
    assert payload["provider_source"]["wrote_alb"] is True
    assert payload["provider_source"]["wrote_phb"] is True
    assert payload["provider_source"]["phb_diagnostic_only_marker"] is True
    assert payload["source_staging"]["ok"] is True
    assert payload["source_staging"]["diagnostic_only"] is True
    assert payload["source_staging"]["gate_candidate"] is False
    assert payload["source_staging"]["aliases_child"] is False
    assert payload["source_staging"]["counts"]["staged_value_count"] == 12.0
    assert payload["source_staging"]["counts"]["invalid_exposed_value_count"] == 0.0
    assert payload["source_staging"]["counts"]["exposed_region_count"] == 1.0
    assert payload["source_child_delta"]["values_identical"] is False
    assert payload["source_child_delta"]["compared_value_count"] == 12.0
    assert payload["source_child_delta"]["differing_value_count"] == 6.0
    assert payload["largest_delta_field"] == {"field": "ALB", "max_abs_diff": 6.0}
    assert any("not strict-gate evidence" in reason for reason in payload["reasons"])


def test_missing_required_provider_source_metadata_fails_closed(tmp_path: Path) -> None:
    candidate = tmp_path / "missing_provider.nc"
    _write_netcdf(candidate, _staging_attrs() | _delta_attrs())

    payload = json.loads(
        report_to_json(audit_diagnostic_adapter_provider_source(candidate))
    )

    assert payload["status"] == "failed"
    assert payload["gate_evidence"] is False
    assert payload["provider_source"]["present"] is False
    assert "VERSION" in payload["provider_source"]["missing_required_attrs"]
    assert any("provider_source missing required attrs" in reason for reason in payload["reasons"])


def test_unsafe_guard_values_fail_closed(tmp_path: Path) -> None:
    candidate = tmp_path / "unsafe.nc"
    attrs = _valid_attrs()
    attrs[f"{PROVIDER}GATE_CANDIDATE"] = _bool(True)
    attrs[f"{STAGING}USES_DIRECT_P_SHORTCUT"] = _bool(True)
    attrs[f"{DELTA}WRITES_CANDIDATE"] = _bool(True)
    _write_netcdf(candidate, attrs)

    payload = json.loads(
        report_to_json(audit_diagnostic_adapter_provider_source(candidate))
    )

    assert payload["status"] == "failed"
    assert payload["provider_source"]["gate_candidate"] is True
    assert payload["source_staging"]["uses_direct_p_shortcut"] is True
    assert payload["source_child_delta"]["writes_candidate"] is True
    assert "provider_source.gate_candidate is True, expected False" in payload["reasons"]
    assert (
        "source_staging.uses_direct_p_shortcut is True, expected False"
        in payload["reasons"]
    )
    assert (
        "source_child_delta.writes_candidate is True, expected False"
        in payload["reasons"]
    )


def test_largest_delta_field_uses_per_field_max_abs_diff(tmp_path: Path) -> None:
    candidate = tmp_path / "largest_delta.nc"
    attrs = _provider_attrs() | _staging_attrs() | _delta_attrs(
        max_by_field={"PHB": 3.0, "MUB": 4.0, "HT": 1.0, "PB": 9.5, "T_INIT": 2.0, "ALB": 8.0}
    )
    _write_netcdf(candidate, attrs)

    payload = json.loads(
        report_to_json(audit_diagnostic_adapter_provider_source(candidate))
    )

    assert payload["status"] == "passed"
    assert payload["largest_delta_field"] == {"field": "PB", "max_abs_diff": 9.5}
    assert payload["source_child_delta"]["fields"]["PB"]["max_abs_diff"] == 9.5
