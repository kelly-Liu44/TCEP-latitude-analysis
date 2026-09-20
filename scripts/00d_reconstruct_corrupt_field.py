"""Reconstruct one unreadable three-hour MSWEP field from valid neighbours."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def read_valid_field(path: Path) -> xr.DataArray:
    """Load one source accumulation and validate its units and grid shape."""
    with xr.open_dataset(path, decode_times=True) as dataset:
        if "precipitation" not in dataset:
            raise KeyError(f"Missing precipitation variable: {path}")
        field = (
            dataset["precipitation"].squeeze(drop=True).transpose("lat", "lon").load()
        )
        units = str(field.attrs.get("units", "")).lower().replace(" ", "")
        if units not in {"mm/3h", "mm3h-1", "mm/3hr", "mm/3hour"}:
            raise ValueError(
                f"Unexpected precipitation units {field.attrs.get('units')!r}: {path}"
            )
        if field.shape != (1800, 3600):
            raise ValueError(f"Unexpected MSWEP grid shape {field.shape}: {path}")
        return field


def main() -> None:
    """Linearly interpolate the target time between two readable source fields."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--before-name", required=True)
    parser.add_argument("--target-name", required=True)
    parser.add_argument("--after-name", required=True)
    parser.add_argument("--before-time", required=True)
    parser.add_argument("--target-time", required=True)
    parser.add_argument("--after-time", required=True)
    args = parser.parse_args()

    before_path = args.source_dir / args.before_name
    corrupt_path = args.source_dir / args.target_name
    after_path = args.source_dir / args.after_name
    target_path = args.output_dir / args.target_name
    before_time = pd.to_datetime(args.before_time, utc=True)
    target_time = pd.to_datetime(args.target_time, utc=True)
    after_time = pd.to_datetime(args.after_time, utc=True)
    if not before_time < target_time < after_time:
        raise ValueError("Target time must lie strictly between the two source times")

    before = read_valid_field(before_path)
    after = read_valid_field(after_path)
    if not np.array_equal(before["lat"], after["lat"]) or not np.array_equal(
        before["lon"], after["lon"]
    ):
        raise ValueError("The two bracketing fields have different coordinates")
    after_weight = (target_time - before_time) / (after_time - before_time)
    before_weight = 1.0 - after_weight
    values = before_weight * np.asarray(
        before, dtype=np.float64
    ) + after_weight * np.asarray(after, dtype=np.float64)

    reconstructed = xr.Dataset(
        data_vars={
            "precipitation": (
                ("time", "lat", "lon"),
                values[np.newaxis].astype(np.float32),
                {
                    "units": "mm 3h-1",
                    "interpolation": "linear between nearest valid fields bracketing an unreadable source field",
                    "before_weight": float(before_weight),
                    "after_weight": float(after_weight),
                    "source_before": before_path.name,
                    "source_after": after_path.name,
                    "unreadable_source": corrupt_path.name,
                },
            )
        },
        coords={
            "time": ("time", [target_time.to_datetime64()], {"long_name": "time"}),
            "lat": before["lat"],
            "lon": before["lon"],
        },
        attrs={
            "history": "Reconstructed unreadable MSWEP V2.8 field by linear time interpolation",
            "is_interpolated": "true",
        },
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temporary = target_path.with_suffix(".nc.tmp")
    if temporary.exists():
        temporary.unlink()
    reconstructed.to_netcdf(
        temporary,
        engine="netcdf4",
        encoding={
            "precipitation": {
                "dtype": "float32",
                "zlib": True,
                "complevel": 2,
                "_FillValue": -9999.0,
            }
        },
    )
    temporary.replace(target_path)

    check = read_valid_field(target_path)
    if not np.isfinite(np.asarray(check)).all():
        raise ValueError(
            f"Reconstructed field contains non-finite values: {target_path}"
        )
    print(f"Created and validated reconstructed field: {target_path}")


if __name__ == "__main__":
    main()
