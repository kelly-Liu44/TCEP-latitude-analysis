"""Reconstruct the unreadable 2011-09-30 03 UTC MSWEP field."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


def read_valid_field(path: Path) -> xr.DataArray:
    """Load one complete three-hour accumulation and validate its units."""
    with xr.open_dataset(path, decode_times=True) as dataset:
        field = (
            dataset["precipitation"].squeeze(drop=True).transpose("lat", "lon").load()
        )
        units = str(field.attrs.get("units", "")).lower().replace(" ", "")
        if units not in {"mm/3h", "mm3h-1", "mm/3hr", "mm/3hour"}:
            raise ValueError(
                f"Unexpected precipitation units {field.attrs.get('units')!r}: {path}"
            )
        return field


def main() -> None:
    """Average the valid 00 and 06 UTC fields at every grid cell."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    before_path = args.source_dir / "2011273.00.nc"
    corrupt_path = args.source_dir / "2011273.03.nc"
    after_path = args.source_dir / "2011273.06.nc"
    target_path = args.output_dir / corrupt_path.name
    before = read_valid_field(before_path)
    after = read_valid_field(after_path)
    if not np.array_equal(before["lat"], after["lat"]) or not np.array_equal(
        before["lon"], after["lon"]
    ):
        raise ValueError("The two bracketing fields have different coordinates")

    # 03 UTC lies exactly halfway between the valid 00 and 06 UTC fields.
    values = 0.5 * np.asarray(before, dtype=np.float64) + 0.5 * np.asarray(
        after, dtype=np.float64
    )
    reconstructed = xr.Dataset(
        data_vars={
            "precipitation": (
                ("time", "lat", "lon"),
                values[np.newaxis].astype(np.float32),
                {
                    "units": "mm 3h-1",
                    "interpolation": "linear between nearest valid fields bracketing an unreadable source field",
                    "before_weight": 0.5,
                    "after_weight": 0.5,
                    "source_before": before_path.name,
                    "source_after": after_path.name,
                    "unreadable_source": corrupt_path.name,
                },
            )
        },
        coords={
            "time": (
                "time",
                [np.datetime64("2011-09-30T03:00:00")],
                {"long_name": "time"},
            ),
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

    # Reopen and read every value so a damaged replacement cannot enter the run.
    check = read_valid_field(target_path)
    if check.shape != (1800, 3600) or not np.isfinite(np.asarray(check)).all():
        raise ValueError(f"Reconstructed field failed validation: {target_path}")
    print(f"Created and validated reconstructed field: {target_path}")


if __name__ == "__main__":
    main()
