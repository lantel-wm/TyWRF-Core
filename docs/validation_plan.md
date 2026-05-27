# Validation Plan

Validation is based on 6 h cycle tests rather than requiring a full 168 h free
forecast to remain field-close.

## Cycle Procedure

1. Start from the matching WRF reference initial state.
2. Integrate a 6 h TyWRF-Core segment.
3. Write WRF-compatible core `wrfout` files.
4. Compare against the WRF reference output for the same valid time.
5. Repeat over multiple 6 h windows.

## Default Field Thresholds

- `U`, `V`, `T`, `PH`, `MU`, `P`, `QVAPOR`: normalized RMSE <= `5%`
- `W` and hydrometeors: normalized RMSE <= `10%`
- `PSFC`, `U10`, `V10`, `T2`, `Q2`: normalized RMSE <= `10%`

## Typhoon Diagnostics

- d02 storm center error <= `20 km`
- MSLP error <= `5 hPa`
- Vmax error <= `5 m s-1`

## Python Validation Tools

`tools/extract_reference.py` inspects WRF NetCDF files and reports dimensions,
present core variables, missing core variables, and metadata for each present
core variable. It exits nonzero when required core variables are missing unless
the expected variable list is narrowed with `--core-variables`.

`tools/compare_wrfout.py` compares selected numeric fields between a WRF
reference file and a TyWRF-Core candidate file. The JSON report contains:

- per-variable status: `ok`, `threshold_exceeded`, `missing_reference`,
  `missing_candidate`, `shape_mismatch`, `non_numeric`, or `no_valid_samples`;
- RMSE, normalized RMSE, max absolute error, threshold, sample counts, and
  reference/candidate shapes when metrics are available;
- summary counts and top-level pass/fail status.

Default normalized RMSE thresholds follow the field thresholds above. Use
`--no-thresholds` to report metrics without threshold failures, or repeat
`--threshold VARIABLE=VALUE` to override individual thresholds.

`tools/run_6h_cycle_test.py --dry-run` resolves the expected
`wrfout_d01_YYYY-MM-DD_HH:MM:SS` and `wrfout_d02_YYYY-MM-DD_HH:MM:SS` files for
the cycle start and end times without invoking the integrator. The dry-run plan
reports expected WRF reference and TyWRF-Core candidate files for both cycle
endpoints and both domains by default.

Example:

```bash
UV_CACHE_DIR=.uv-cache UV_PYTHON_INSTALL_DIR=.uv-python uv run python tools/run_6h_cycle_test.py \
  --reference-dir /path/to/reference/WRF \
  --candidate-dir /path/to/tywrf/output \
  --start 2025-07-26_00:00:00 \
  --dry-run \
  --pretty
```

TC diagnostics remain explicitly pending in JSON reports until a tested center,
MSLP, Vmax, and rainfall diagnostic hook is added.
