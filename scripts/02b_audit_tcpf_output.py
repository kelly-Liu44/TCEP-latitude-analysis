"""Audit one yearly TCPF Parquet file before downstream calculations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
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
ALLOWED_INTERPOLATED_TIMES = {
    pd.Timestamp(year=2020, month=12, day=31, hour=hour, tz="UTC")
    for hour in (3, 6, 9, 12, 15, 18, 21)
}
ALLOWED_INTERPOLATED_TIMES.add(pd.Timestamp("2011-09-30T03:00:00Z"))
ALLOWED_INTERPOLATED_TIMES.add(pd.Timestamp("2023-05-12T21:00:00Z"))


def audit(path: Path, expected_year: int) -> dict[str, float | int]:
    """Check schema and stream row groups so the audit does not load a full year."""
    parquet = pq.ParquetFile(path)
    missing = REQUIRED_COLUMNS.difference(parquet.schema.names)
    if missing:
        raise ValueError(f"Missing TCPF columns: {sorted(missing)}")
    row_count = 0
    minimum_precipitation = np.inf
    maximum_unit_error = 0.0
    maximum_latitude_error = 0.0
    maximum_longitude_error = 0.0
    years_seen: set[int] = set()
    basins_seen: set[str] = set()
    interpolated_times: set[pd.Timestamp] = set()

    # Audit each physical Parquet row group independently to bound memory use.
    for row_group in range(parquet.num_row_groups):
        frame = parquet.read_row_group(
            row_group, columns=sorted(REQUIRED_COLUMNS)
        ).to_pandas()
        row_count += len(frame)
        if frame.empty:
            continue
        if (
            frame[
                [
                    "SID",
                    "TIME",
                    "BASIN",
                    "TC_LAT",
                    "TC_LON_360",
                    "GRID_LAT",
                    "GRID_LON",
                    "PRECIP_3H_MM",
                    "PRECIP_MMHR",
                    "ROW",
                ]
            ]
            .isna()
            .any()
            .any()
        ):
            raise ValueError(f"Null required values in row group {row_group}")
        years_seen.update(frame["YEAR"].astype(int).unique().tolist())
        basins_seen.update(frame["BASIN"].astype(str).unique().tolist())
        flags = frame["IS_INTERPOLATED_TIMESTEP"]
        if not pd.api.types.is_bool_dtype(flags.dtype):
            raise ValueError(
                f"Interpolation provenance flag is not Boolean in row group {row_group}"
            )
        timestamps = pd.to_datetime(frame["TIME"], utc=True)
        expected_flags = timestamps.isin(ALLOWED_INTERPOLATED_TIMES)
        if not np.array_equal(flags.to_numpy(dtype=bool), expected_flags):
            raise ValueError(
                f"Interpolation provenance does not match source timestamps in row group {row_group}"
            )
        interpolated_times.update(timestamps.loc[flags].tolist())
        precipitation = frame["PRECIP_3H_MM"].to_numpy(dtype=np.float64)
        rate = frame["PRECIP_MMHR"].to_numpy(dtype=np.float64)
        if np.any(precipitation <= np.float64(np.float32(0.3))):
            raise ValueError(
                f"TCPF contains a value that does not satisfy strict >0.3 mm/3h in row group {row_group}"
            )
        minimum_precipitation = min(minimum_precipitation, float(precipitation.min()))
        maximum_unit_error = max(
            maximum_unit_error, float(np.max(np.abs(precipitation / rate - 3.0)))
        )
        rows = frame["ROW"].to_numpy(dtype=np.int64)
        expected_latitude = 89.95 - (rows // 3600) * 0.1
        expected_longitude = (rows % 3600) * 0.1 + 0.05
        maximum_latitude_error = max(
            maximum_latitude_error,
            float(
                np.max(
                    np.abs(frame["GRID_LAT"].to_numpy(dtype=float) - expected_latitude)
                )
            ),
        )
        maximum_longitude_error = max(
            maximum_longitude_error,
            float(
                np.max(
                    np.abs(frame["GRID_LON"].to_numpy(dtype=float) - expected_longitude)
                )
            ),
        )
        if np.any((frame["TC_LON_360"] < 0) | (frame["TC_LON_360"] >= 360)):
            raise ValueError(f"TC longitude outside [0, 360) in row group {row_group}")
        if np.any((frame["GRID_LON"] < 0) | (frame["GRID_LON"] >= 360)):
            raise ValueError(
                f"Grid longitude outside [0, 360) in row group {row_group}"
            )

    if row_count == 0:
        raise ValueError("TCPF output contains no rows")
    if years_seen != {expected_year}:
        raise ValueError(f"Expected only {expected_year}; found {sorted(years_seen)}")
    if maximum_unit_error > 1e-5:
        raise ValueError(f"Precipitation unit-ratio error is {maximum_unit_error}")
    if maximum_latitude_error > 3e-4 or maximum_longitude_error > 3e-4:
        raise ValueError(
            f"Grid coordinate residuals are too large: lat={maximum_latitude_error}, "
            f"lon={maximum_longitude_error}"
        )
    if "NA" not in basins_seen and expected_year in range(1980, 2024):
        print(f"WARNING: {expected_year} contains no North Atlantic TCPF rows")
    if not interpolated_times.issubset(ALLOWED_INTERPOLATED_TIMES):
        raise ValueError(
            f"Unexpected interpolated timestamps: {sorted(interpolated_times - ALLOWED_INTERPOLATED_TIMES)}"
        )
    return {
        "year": expected_year,
        "rows": row_count,
        "minimum_precipitation_mm_3h": minimum_precipitation,
        "maximum_unit_ratio_error": maximum_unit_error,
        "maximum_latitude_residual_deg": maximum_latitude_error,
        "maximum_longitude_residual_deg": maximum_longitude_error,
        "basin_count": len(basins_seen),
        "interpolated_timestep_count": len(interpolated_times),
    }


def main() -> None:
    """Audit one file and optionally append one concise CSV record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = audit(args.input, args.year)
    print("PASS " + ", ".join(f"{key}={value}" for key, value in result.items()))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        record = pd.DataFrame([result])
        if args.report.exists():
            existing = pd.read_csv(args.report)
            existing = existing.loc[existing["year"] != args.year]
            record = pd.concat([existing, record], ignore_index=True).sort_values(
                "year"
            )
        temporary = args.report.with_suffix(args.report.suffix + ".tmp")
        record.to_csv(temporary, index=False)
        temporary.replace(args.report)


if __name__ == "__main__":
    main()
