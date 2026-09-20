"""Audit one annual TCEP event table and its strict-POT99 extreme cells."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


EVENT_COLUMNS = {
    "SID",
    "TIME",
    "YEAR",
    "HEMISPHERE",
    "BASIN",
    "TC_LAT",
    "TCEP_INTENSITY_MM_3H",
    "N_EXTREME_CELLS",
    "IS_INTERPOLATED_TIMESTEP",
}
CELL_COLUMNS = {
    "SID",
    "TIME",
    "YEAR",
    "BASIN",
    "TC_LAT",
    "PRECIP_3H_MM",
    "POT99_MM_3H",
    "ROW",
    "IS_INTERPOLATED_TIMESTEP",
}


def audit(
    events_path: Path, cells_path: Path, expected_year: int
) -> dict[str, int | float]:
    """Stream extreme cells and prove that their aggregation matches events."""
    events = pd.read_parquet(events_path)
    missing_events = EVENT_COLUMNS.difference(events.columns)
    if missing_events:
        raise ValueError(f"Missing event columns: {sorted(missing_events)}")
    if events.empty:
        raise ValueError(f"No TCEP events in {expected_year}")
    if set(events["YEAR"].astype(int).unique()) != {expected_year}:
        raise ValueError(f"Event YEAR values do not match {expected_year}")
    if events.duplicated(["SID", "TIME"]).any():
        raise ValueError("Event table contains duplicate SID + TIME keys")
    if not set(events["HEMISPHERE"].astype(str)).issubset({"NH", "SH"}):
        raise ValueError("Unexpected hemisphere code")
    expected_hemisphere = np.where(
        events["TC_LAT"].to_numpy(dtype=float) >= 0, "NH", "SH"
    )
    if not np.array_equal(
        events["HEMISPHERE"].astype(str).to_numpy(), expected_hemisphere
    ):
        raise ValueError("Hemisphere is inconsistent with TC latitude")
    if not pd.api.types.is_bool_dtype(events["IS_INTERPOLATED_TIMESTEP"].dtype):
        raise ValueError("Event interpolation flag is not Boolean")

    source = pq.ParquetFile(cells_path)
    missing_cells = CELL_COLUMNS.difference(source.schema.names)
    if missing_cells:
        raise ValueError(f"Missing extreme-cell columns: {sorted(missing_cells)}")
    partials = []
    cell_count = 0
    minimum_margin = np.inf
    for row_group in range(source.num_row_groups):
        frame = source.read_row_group(
            row_group, columns=sorted(CELL_COLUMNS)
        ).to_pandas()
        if frame.empty:
            continue
        if set(frame["YEAR"].astype(int).unique()) != {expected_year}:
            raise ValueError(f"Extreme-cell YEAR values do not match {expected_year}")
        precipitation = frame["PRECIP_3H_MM"].to_numpy(dtype=np.float64)
        thresholds = frame["POT99_MM_3H"].to_numpy(dtype=np.float64)
        if np.any(~np.isfinite(precipitation)) or np.any(~np.isfinite(thresholds)):
            raise ValueError(
                f"Non-finite extreme value or threshold in row group {row_group}"
            )
        margins = precipitation - thresholds
        if np.any(margins <= 0):
            raise ValueError(f"Non-strict POT99 exceedance in row group {row_group}")
        if not pd.api.types.is_bool_dtype(frame["IS_INTERPOLATED_TIMESTEP"].dtype):
            raise ValueError(
                f"Cell interpolation flag is not Boolean in row group {row_group}"
            )
        minimum_margin = min(minimum_margin, float(margins.min()))
        cell_count += len(frame)
        work = frame.assign(INTENSITY_SUM_MM_3H=precipitation)
        partials.append(
            work.groupby(["SID", "TIME"], as_index=False).agg(
                INTENSITY_SUM_MM_3H=("INTENSITY_SUM_MM_3H", "sum"),
                N_EXTREME_CELLS=("PRECIP_3H_MM", "size"),
                IS_INTERPOLATED_TIMESTEP=("IS_INTERPOLATED_TIMESTEP", "max"),
            )
        )
    if not partials:
        raise ValueError(f"No extreme cells in {expected_year}")

    reconstructed = (
        pd.concat(partials, ignore_index=True)
        .groupby(["SID", "TIME"], as_index=False)
        .agg(
            INTENSITY_SUM_MM_3H=("INTENSITY_SUM_MM_3H", "sum"),
            N_EXTREME_CELLS=("N_EXTREME_CELLS", "sum"),
            IS_INTERPOLATED_TIMESTEP=("IS_INTERPOLATED_TIMESTEP", "max"),
        )
    )
    reconstructed["TCEP_INTENSITY_MM_3H"] = (
        reconstructed["INTENSITY_SUM_MM_3H"] / reconstructed["N_EXTREME_CELLS"]
    )
    comparison = events.merge(
        reconstructed,
        on=["SID", "TIME"],
        how="outer",
        validate="one_to_one",
        suffixes=("_EVENT", "_CELL"),
        indicator=True,
    )
    if not comparison["_merge"].eq("both").all():
        raise ValueError("Event and extreme-cell SID + TIME keys differ")
    if not np.array_equal(
        comparison["N_EXTREME_CELLS_EVENT"].to_numpy(dtype=np.int64),
        comparison["N_EXTREME_CELLS_CELL"].to_numpy(dtype=np.int64),
    ):
        raise ValueError("Event cell counts do not match extreme-cell aggregation")
    if not np.allclose(
        comparison["TCEP_INTENSITY_MM_3H_EVENT"],
        comparison["TCEP_INTENSITY_MM_3H_CELL"],
        rtol=1e-10,
        atol=1e-10,
    ):
        raise ValueError("Event mean intensities do not match extreme-cell aggregation")
    if not np.array_equal(
        comparison["IS_INTERPOLATED_TIMESTEP_EVENT"].to_numpy(dtype=bool),
        comparison["IS_INTERPOLATED_TIMESTEP_CELL"].to_numpy(dtype=bool),
    ):
        raise ValueError("Interpolation provenance differs between cells and events")
    return {
        "year": expected_year,
        "events": len(events),
        "extreme_cells": cell_count,
        "minimum_exceedance_margin_mm_3h": minimum_margin,
        "interpolated_events": int(events["IS_INTERPOLATED_TIMESTEP"].sum()),
    }


def main() -> None:
    """Audit one year and update an atomic cross-year CSV report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.events, args.cells, args.year)
    record = pd.DataFrame([result])
    args.report.parent.mkdir(parents=True, exist_ok=True)
    if args.report.exists():
        existing = pd.read_csv(args.report)
        existing = existing.loc[existing["year"] != args.year]
        record = pd.concat([existing, record], ignore_index=True).sort_values("year")
    temporary = args.report.with_suffix(args.report.suffix + ".tmp")
    record.to_csv(temporary, index=False)
    temporary.replace(args.report)
    print("PASS " + ", ".join(f"{key}={value}" for key, value in result.items()))


if __name__ == "__main__":
    main()
