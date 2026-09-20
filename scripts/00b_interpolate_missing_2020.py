"""Reconstruct the seven missing 2020-12-31 MSWEP fields by time interpolation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


MISSING_HOURS = (3, 6, 9, 12, 15, 18, 21)


def source_name(year: int, day_of_year: int, hour: int) -> str:
    """Build one MSWEP filename from its year, day of year, and UTC hour."""
    return f"{year}{day_of_year:03d}.{hour:02d}.nc"


def read_precipitation(path: Path) -> xr.DataArray:
    """Load one source field and retain its coordinates and unit metadata."""
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
        return field


def interpolate_one(source_dir: Path, output_dir: Path, hour: int) -> Path:
    """Interpolate between the nearest valid fields that bracket the data gap."""
    # Seven consecutive fields are absent between these two valid endpoints.
    before_path = source_dir / source_name(2020, 366, 0)
    after_path = source_dir / source_name(2021, 1, 0)
    target_path = output_dir / source_name(2020, 366, hour)
    if not before_path.exists() or not after_path.exists():
        raise FileNotFoundError(
            f"Interpolation sources are missing: {before_path}, {after_path}"
        )

    # Load both complete fields before closing their source files.
    before = read_precipitation(before_path)
    after = read_precipitation(after_path)
    if not np.array_equal(before["lat"], after["lat"]) or not np.array_equal(
        before["lon"], after["lon"]
    ):
        raise ValueError(
            f"Coordinate mismatch between {before_path.name} and {after_path.name}"
        )

    # The target is hour/24 of the way from 31 December to 1 January.
    before_values = np.asarray(before.to_numpy(), dtype=np.float64)
    after_values = np.asarray(after.to_numpy(), dtype=np.float64)
    after_weight = hour / 24.0
    before_weight = 1.0 - after_weight
    values = before_weight * before_values + after_weight * after_values
    target_time = np.datetime64(f"2020-12-31T{hour:02d}:00:00")
    interpolated = xr.Dataset(
        data_vars={
            "precipitation": (
                ("time", "lat", "lon"),
                values[np.newaxis, :, :].astype(np.float32),
                {
                    "units": "mm 3h-1",
                    "interpolation": "linear between nearest valid fields bracketing the gap",
                    "before_weight": before_weight,
                    "after_weight": after_weight,
                    "source_before": before_path.name,
                    "source_after": after_path.name,
                },
            )
        },
        coords={
            "time": ("time", [target_time], {"long_name": "time"}),
            "lat": before["lat"],
            "lon": before["lon"],
        },
        attrs={
            "history": "Reconstructed missing MSWEP V2.8 field by linear time interpolation",
            "is_interpolated": "true",
        },
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    interpolated.to_netcdf(
        target_path,
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
    return target_path


def main() -> None:
    """Create all seven replacement files and print their provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    for hour in MISSING_HOURS:
        path = interpolate_one(args.source_dir, args.output_dir, hour)
        print(f"Created interpolated field: {path}")


if __name__ == "__main__":
    main()
