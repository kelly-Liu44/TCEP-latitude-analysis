"""Verify completeness and consistency of the full 1980--2023 TCPF archive."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


REQUIRED_COLUMNS = {
    "SID",
    "TIME",
    "YEAR",
    "BASIN",
    "TC_LAT",
    "TC_LON_360",
    "GRID_LAT",
    "GRID_LON",
    "PRECIP_3H_MM",
    "PRECIP_MMHR",
    "ROW",
    "IS_INTERPOLATED_TIMESTEP",
}


def main() -> None:
    """Cross-check annual files against their completed streaming audits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--annual-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    args = parser.parse_args()

    expected_years = list(range(args.start_year, args.end_year + 1))
    expected_names = {f"tc_gridded_{year}.parquet" for year in expected_years}
    actual_names = {path.name for path in args.input_dir.glob("tc_gridded_*.parquet")}
    missing_files = sorted(expected_names - actual_names)
    extra_files = sorted(actual_names - expected_names)
    if missing_files or extra_files:
        raise ValueError(
            f"TCPF file set mismatch; missing={missing_files}, extra={extra_files}"
        )
    if not args.annual_report.exists():
        raise FileNotFoundError(args.annual_report)

    annual = pd.read_csv(args.annual_report)
    if annual["year"].duplicated().any():
        raise ValueError("Annual audit report contains duplicate years")
    if sorted(annual["year"].astype(int).tolist()) != expected_years:
        raise ValueError("Annual audit report does not contain exactly 1980--2023")
    annual = annual.set_index(annual["year"].astype(int), drop=False)

    rows = []
    for year in expected_years:
        path = args.input_dir / f"tc_gridded_{year}.parquet"
        parquet = pq.ParquetFile(path)
        missing_columns = REQUIRED_COLUMNS.difference(parquet.schema.names)
        if missing_columns:
            raise ValueError(f"{year}: missing columns {sorted(missing_columns)}")
        metadata_rows = int(parquet.metadata.num_rows)
        audited_rows = int(annual.loc[year, "rows"])
        if metadata_rows != audited_rows:
            raise ValueError(
                f"{year}: Parquet rows {metadata_rows} != audited rows {audited_rows}"
            )
        interpolation_count = int(annual.loc[year, "interpolated_timestep_count"])
        expected_interpolation_count = {2011: 1, 2020: 7, 2023: 1}.get(year, 0)
        if interpolation_count != expected_interpolation_count:
            raise ValueError(
                f"{year}: expected {expected_interpolation_count} represented interpolated timestamps, "
                f"got {interpolation_count}"
            )
        rows.append(
            {
                "YEAR": year,
                "FILE": path.name,
                "FILE_SIZE_BYTES": path.stat().st_size,
                "PARQUET_ROWS": metadata_rows,
                "ROW_GROUPS": parquet.num_row_groups,
                "MIN_PRECIP_3H_MM": float(
                    annual.loc[year, "minimum_precipitation_mm_3h"]
                ),
                "INTERPOLATED_TIMESTEP_COUNT": interpolation_count,
                "STATUS": "PASS",
            }
        )

    result = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    result.to_csv(temporary, index=False)
    temporary.replace(args.output)
    print(
        f"PASS: {len(result)} yearly TCPF files, "
        f"{result['PARQUET_ROWS'].sum():,} total attributed cells; report={args.output}"
    )


if __name__ == "__main__":
    main()
