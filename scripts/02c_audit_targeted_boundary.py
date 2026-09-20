"""Audit targeted TCPF output for date-line and Southern Hemisphere failures."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tcep.grid import grid_row, haversine_km, wrapped_delta


def audit(
    path: Path,
    storm_id: str,
    hemisphere: str,
    require_dateline_crossing: bool,
    require_zero_seam_crossing: bool,
) -> pd.DataFrame:
    """Check one targeted output and return one concise row per storm time step."""
    frame = pd.read_parquet(path)
    if frame.empty:
        raise ValueError(f"Targeted output is empty: {path}")
    if set(frame["SID"].astype(str)) != {storm_id}:
        raise ValueError(
            f"Expected only {storm_id}; found {sorted(frame['SID'].astype(str).unique())}"
        )
    expected_sign = 1 if hemisphere == "NH" else -1
    if np.any(np.sign(frame["TC_LAT"].to_numpy(dtype=float)) != expected_sign):
        raise ValueError(f"{storm_id} contains a TC latitude outside {hemisphere}")

    # Recompute the canonical row key, including negative Southern Hemisphere latitudes.
    reconstructed_rows = grid_row(frame["GRID_LAT"], frame["GRID_LON"])
    if not np.array_equal(
        reconstructed_rows.astype(np.int64), frame["ROW"].to_numpy(dtype=np.int64)
    ):
        raise ValueError(f"Grid-row round trip failed for {storm_id}")

    # Every attributed cell must remain inside the exact 1000-km storm-centred circle.
    distances = haversine_km(
        frame["TC_LAT"], frame["TC_LON_360"], frame["GRID_LAT"], frame["GRID_LON"]
    )
    if np.nanmax(distances) > 1000.01:
        raise ValueError(f"Attributed grid cell exceeds 1000 km for {storm_id}")

    # Express grid longitude relative to the storm so 0/360 and 180-degree seams
    # cannot create an artificial 300-degree jump in the diagnostic.
    relative_longitude = wrapped_delta(frame["GRID_LON"], frame["TC_LON_360"])
    if np.nanmax(np.abs(relative_longitude)) > 170.0:
        raise ValueError(f"Longitude convention mismatch detected for {storm_id}")

    summary = (
        frame.assign(
            DISTANCE_KM=np.asarray(distances, dtype=float),
            RELATIVE_LON_DEG=np.asarray(relative_longitude, dtype=float),
        )
        .groupby("TIME", as_index=False)
        .agg(
            TC_LAT=("TC_LAT", "first"),
            TC_LON_360=("TC_LON_360", "first"),
            N_ATTRIBUTED_CELLS=("ROW", "size"),
            GRID_LON_MIN=("GRID_LON", "min"),
            GRID_LON_MAX=("GRID_LON", "max"),
            RELATIVE_LON_MIN=("RELATIVE_LON_DEG", "min"),
            RELATIVE_LON_MAX=("RELATIVE_LON_DEG", "max"),
            MAX_DISTANCE_KM=("DISTANCE_KM", "max"),
        )
        .sort_values("TIME")
    )
    if require_dateline_crossing:
        storm_longitudes = summary["TC_LON_360"].to_numpy(dtype=float)
        if not (np.any(storm_longitudes < 180.0) and np.any(storm_longitudes >= 180.0)):
            raise ValueError(
                f"Selected time window does not cross 180 degrees for {storm_id}"
            )
        steps = np.abs(wrapped_delta(storm_longitudes[1:], storm_longitudes[:-1]))
        if np.any(steps > 10.0):
            raise ValueError(
                f"Track contains an artificial longitude jump for {storm_id}: {steps}"
            )
    if require_zero_seam_crossing:
        storm_longitudes = summary["TC_LON_360"].to_numpy(dtype=float)
        if not (np.any(storm_longitudes < 5.0) and np.any(storm_longitudes > 355.0)):
            raise ValueError(
                f"Selected time window does not cross the 0/360 grid seam for {storm_id}"
            )
        steps = np.abs(wrapped_delta(storm_longitudes[1:], storm_longitudes[:-1]))
        if np.any(steps > 10.0):
            raise ValueError(
                f"Track contains an artificial longitude jump for {storm_id}: {steps}"
            )
    return summary


def main() -> None:
    """Parse paths, run the targeted audit, and save human-readable diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--storm-id", required=True)
    parser.add_argument("--hemisphere", choices=["NH", "SH"], required=True)
    parser.add_argument("--require-dateline-crossing", action="store_true")
    parser.add_argument("--require-zero-seam-crossing", action="store_true")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    summary = audit(
        args.input,
        args.storm_id,
        args.hemisphere,
        args.require_dateline_crossing,
        args.require_zero_seam_crossing,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.report, index=False)
    print(f"PASS {args.storm_id}: {len(summary)} time steps; report={args.report}")


if __name__ == "__main__":
    main()
