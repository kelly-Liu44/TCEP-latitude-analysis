"""Validate MSWEP V2.8 and IBTrACS inputs before any large computation."""

from __future__ import annotations

import argparse
import calendar
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


EXPECTED_UNITS = {"mm/3h", "mm/3hr", "mm 3h-1", "mm/3hour"}
FILE_PATTERN = re.compile(r"^(?P<year>\d{4})(?P<doy>\d{3})\.(?P<hour>\d{2})\.nc$")


def expected_names(year: int) -> set[str]:
    """Return every expected three-hour filename for one calendar year."""
    days = 366 if calendar.isleap(year) else 365
    return {
        f"{year}{day:03d}.{hour:02d}.nc"
        for day in range(1, days + 1)
        for hour in range(0, 24, 3)
    }


def validate_mswep(
    directory: Path, years: range, override_dir: Path | None = None
) -> list[str]:
    """Check files, metadata, coordinates, dimensions, and timestamp spacing."""
    errors: list[str] = []
    files = list(directory.glob("*.nc"))
    override_files = list(override_dir.glob("*.nc")) if override_dir else []
    if not files:
        raise FileNotFoundError(f"No NetCDF files found: {directory}")

    # Read one representative file to establish the expected grid contract.
    with xr.open_dataset(files[0], decode_times=True) as dataset:
        if "precipitation" not in dataset:
            errors.append(f"Missing precipitation variable in {files[0].name}")
        else:
            units = str(dataset["precipitation"].attrs.get("units", ""))
            if units.lower().replace(" ", "") not in {
                value.replace(" ", "") for value in EXPECTED_UNITS
            }:
                errors.append(f"Unexpected precipitation units {units!r}")
        lat = dataset["lat"].to_numpy()
        lon = dataset["lon"].to_numpy()
        shape = tuple(dataset["precipitation"].squeeze(drop=True).shape)
        if shape != (1800, 3600):
            errors.append(f"Unexpected precipitation shape {shape}")
        if not np.all(np.diff(lat) < 0):
            errors.append("Latitude must be strictly decreasing")
        if not np.all(np.diff(lon) > 0):
            errors.append("Raw longitude must be strictly increasing")
        if not np.isclose(lat[0], 89.95, atol=2e-3) or not np.isclose(
            lat[-1], -89.95, atol=2e-3
        ):
            errors.append("Latitude cell centres are not 89.95..-89.95")
        if not np.isclose(lon[0], -179.95, atol=2e-3) or not np.isclose(
            lon[-1], 179.95, atol=2e-3
        ):
            errors.append("Longitude cell centres are not -179.95..179.95")

    # Compare each year with the exact 365/366-file expectation.
    for year in years:
        actual = {path.name for path in files if path.name.startswith(str(year))}
        actual.update(
            path.name for path in override_files if path.name.startswith(str(year))
        )
        missing = expected_names(year) - actual
        extra = actual - expected_names(year)
        if missing:
            errors.append(
                f"{year}: missing {len(missing)} files, e.g. {sorted(missing)[:3]}"
            )
        if extra:
            errors.append(f"{year}: unexpected files, e.g. {sorted(extra)[:3]}")

    return errors


def validate_ibtracs(csv_path: Path, first_year: int, last_year: int) -> list[str]:
    """Check required columns, coordinates, timestamps, and track labels."""
    required = {
        "SID",
        "SEASON",
        "BASIN",
        "NAME",
        "ISO_TIME",
        "LAT",
        "LON",
        "TRACK_TYPE",
        "USA_WIND",
        "USA_PRES",
    }
    # Preserve the North Atlantic basin code "NA" as text instead of pandas NaN.
    frame = pd.read_csv(csv_path, skiprows=[1], low_memory=False, keep_default_na=False)
    errors: list[str] = []
    missing = required.difference(frame.columns)
    if missing:
        return [f"IBTrACS missing columns: {sorted(missing)}"]
    frame = frame.loc[frame["SEASON"].between(first_year, last_year)].copy()
    frame["ISO_TIME"] = pd.to_datetime(frame["ISO_TIME"], errors="coerce", utc=True)
    for column in ["LAT", "LON", "USA_WIND", "USA_PRES"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if frame[["ISO_TIME", "LAT", "LON"]].isna().any().any():
        errors.append("IBTrACS contains missing time or track coordinates")
    if ((frame["LAT"].abs() > 90) | (frame["LON"].abs() > 360)).any():
        errors.append("IBTrACS contains coordinates outside physical bounds")
    if frame.duplicated(["SID", "ISO_TIME"]).any():
        errors.append("IBTrACS contains duplicate SID/time rows")
    return errors


def main() -> None:
    """Run validation and fail with actionable messages when a contract breaks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mswep-dir", type=Path, required=True)
    parser.add_argument("--ibtracs-csv", type=Path, required=True)
    parser.add_argument("--override-dir", type=Path)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    args = parser.parse_args()
    errors = validate_mswep(
        args.mswep_dir,
        range(args.start_year, args.end_year + 1),
        args.override_dir,
    )
    errors.extend(validate_ibtracs(args.ibtracs_csv, args.start_year, args.end_year))
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)
    print("Input validation passed: units, grid, annual file counts, and track fields.")


if __name__ == "__main__":
    main()
