"""Build one yearly TCPF table from MSWEP V2.8 and prepared IBTrACS tracks."""

from __future__ import annotations

import argparse
import calendar
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tcep.grid import grid_row, lon360, wrapped_delta
from tcep.tcpf import Parameters, attributed_mask
from tcep.units import accumulation_to_rate, validate_units


def parse_args() -> argparse.Namespace:
    """Parse explicit input and output paths so no workstation path is hidden."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--mswep-dir", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--override-dir", type=Path)
    parser.add_argument("--variable", default="precipitation")
    parser.add_argument(
        "--max-times",
        type=int,
        help="Process only the first N timestamps for a smoke test",
    )
    parser.add_argument(
        "--storm-id",
        action="append",
        help="Restrict extraction to one SID; repeat for more storms",
    )
    parser.add_argument(
        "--start-time", help="Optional inclusive UTC start time for a targeted test"
    )
    parser.add_argument(
        "--end-time", help="Optional inclusive UTC end time for a targeted test"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of an existing final output",
    )
    return parser.parse_args()


def expected_files(
    mswep_dir: Path, year: int, override_dir: Path | None = None
) -> dict[str, Path]:
    """Index exactly the expected 3-hour files and fail on missing timestamps."""
    days = 366 if calendar.isleap(year) else 365
    expected = {
        f"{year}{day:03d}.{hour:02d}.nc"
        for day in range(1, days + 1)
        for hour in range(0, 24, 3)
    }
    found = {path.name: path for path in mswep_dir.glob(f"{year}*.nc")}
    # Documented replacement files take precedence without modifying raw data.
    if override_dir:
        found.update({path.name: path for path in override_dir.glob(f"{year}*.nc")})
    missing = sorted(expected - found.keys())
    if missing:
        raise FileNotFoundError(
            f"{year}: missing {len(missing)} MSWEP files, e.g. {missing[:5]}"
        )
    return {name: found[name] for name in expected}


def read_field(path: Path, variable: str):
    """Read one field and return coordinates, data, and interpolation provenance."""
    with xr.open_dataset(path, decode_times=True) as dataset:
        if variable not in dataset:
            raise KeyError(f"{variable!r} is absent from {path.name}")
        field = dataset[variable].squeeze(drop=True)
        validate_units(field.attrs.get("units"))
        if tuple(field.dims) != ("lat", "lon"):
            field = field.transpose("lat", "lon")
        lat = np.asarray(dataset["lat"].to_numpy(), dtype=float)
        raw_lon = np.asarray(dataset["lon"].to_numpy(), dtype=float)
        order = np.argsort(lon360(raw_lon))
        lon = lon360(raw_lon)[order]
        if not np.all(np.diff(lat) < 0) or not np.all(np.diff(lon) > 0):
            raise ValueError(f"Non-monotonic coordinates in {path.name}")
        values = np.asarray(field.to_numpy(), dtype=np.float32)[:, order]
        # Replacement files explicitly carry this dataset-level provenance flag.
        is_interpolated = (
            str(dataset.attrs.get("is_interpolated", "false")).strip().lower() == "true"
        )
    return lat, lon, values, is_interpolated


def search_indices(lat, lon, tc_lat, tc_lon, radius_km=1000.0):
    """Select a rectangular candidate window; the exact circle is applied by TCPF."""
    lat_half_width = radius_km / 111.0
    cos_lat = max(float(np.cos(np.radians(abs(tc_lat)))), 0.08)
    lon_half_width = min(radius_km / (111.0 * cos_lat), 170.0)
    lat_index = np.flatnonzero(np.abs(lat - tc_lat) <= lat_half_width)
    lon_index = np.flatnonzero(np.abs(wrapped_delta(lon, tc_lon)) <= lon_half_width)
    # Put longitudes in TC-centred order so date-line neighbours are adjacent for CCL.
    lon_index = lon_index[np.argsort(wrapped_delta(lon[lon_index], tc_lon))]
    return lat_index, lon_index


def extract(storm, accumulation, lat, lon, timestamp, is_interpolated_timestep=False):
    """Extract attributed cells and preserve both native and rate units."""
    if pd.isna(storm.LAT) or pd.isna(storm.LON):
        return None
    lat_index, lon_index = search_indices(lat, lon, float(storm.LAT), float(storm.LON))
    if not len(lat_index) or not len(lon_index):
        return None
    local_accum = accumulation[np.ix_(lat_index, lon_index)]
    local_accum = np.where(
        np.isfinite(local_accum) & (local_accum >= 0), local_accum, 0.0
    )
    local_lon, local_lat = np.meshgrid(lon[lon_index], lat[lat_index])
    mask = attributed_mask(
        local_accum,
        local_lat,
        local_lon,
        float(storm.LAT),
        float(storm.LON),
        Parameters(),
    )
    if not mask.any():
        return None
    rate = accumulation_to_rate(local_accum[mask])
    frame = pd.DataFrame(
        {
            "SID": str(storm.SID),
            "TIME": timestamp,
            "YEAR": np.int16(timestamp.year),
            "BASIN": str(storm.BASIN),
            "TC_LAT": np.float32(storm.LAT),
            "TC_LON_360": np.float32(storm.LON),
            "GRID_LAT": local_lat[mask].astype(np.float32),
            "GRID_LON": local_lon[mask].astype(np.float32),
            "PRECIP_3H_MM": local_accum[mask].astype(np.float32),
            "PRECIP_MMHR": rate.astype(np.float32),
            # Repeat the time-step provenance on every cell so streamed downstream
            # files remain independently auditable without reopening the NetCDF file.
            "IS_INTERPOLATED_TIMESTEP": bool(is_interpolated_timestep),
        }
    )
    # Store the canonical row key so downstream joins never depend on float text.
    frame["ROW"] = grid_row(frame["GRID_LAT"], frame["GRID_LON"]).astype(np.int32)
    return frame


def main() -> None:
    """Stream one year and atomically replace the requested Parquet output."""
    args = parse_args()
    files = expected_files(args.mswep_dir, args.year, args.override_dir)
    tracks = pd.read_parquet(args.tracks)
    tracks["TIME"] = pd.to_datetime(tracks["TIME"], utc=True)
    tracks = tracks.loc[tracks["TIME"].dt.year == args.year].copy()
    if args.storm_id:
        tracks = tracks.loc[tracks["SID"].isin(args.storm_id)].copy()
    if args.start_time:
        tracks = tracks.loc[
            tracks["TIME"] >= pd.to_datetime(args.start_time, utc=True)
        ].copy()
    if args.end_time:
        tracks = tracks.loc[
            tracks["TIME"] <= pd.to_datetime(args.end_time, utc=True)
        ].copy()
    if tracks.empty:
        raise ValueError(f"No prepared tracks for {args.year}")
    if args.max_times is not None:
        if args.max_times < 1:
            raise ValueError("--max-times must be positive")
        selected_times = (
            tracks["TIME"].drop_duplicates().sort_values().head(args.max_times)
        )
        tracks = tracks.loc[tracks["TIME"].isin(selected_times)].copy()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output already exists; use --overwrite explicitly: {args.output}"
        )
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.parent.mkdir(parents=True, exist_ok=True)
    # An interrupted temporary file is safe to replace; a final file is not.
    if temporary.exists():
        temporary.unlink()
    writer = None
    cached_name, cached_field = None, None
    cached_is_interpolated = False
    lat = lon = None
    try:
        for timestamp, positions in tracks.groupby("TIME", sort=True):
            filename = (
                f"{timestamp.year}{timestamp.dayofyear:03d}.{timestamp.hour:02d}.nc"
            )
            if filename != cached_name:
                lat, lon, cached_field, cached_is_interpolated = read_field(
                    files[filename], args.variable
                )
                cached_name = filename
            batches = [
                extract(
                    storm, cached_field, lat, lon, timestamp, cached_is_interpolated
                )
                for _, storm in positions.iterrows()
            ]
            batches = [batch for batch in batches if batch is not None]
            if not batches:
                continue
            table = pa.Table.from_pandas(
                pd.concat(batches, ignore_index=True), preserve_index=False
            )
            if writer is None:
                writer = pq.ParquetWriter(temporary, table.schema, compression="snappy")
            writer.write_table(table.cast(writer.schema, safe=False))
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        raise RuntimeError(f"No TCPF cells produced for {args.year}")
    temporary.replace(args.output)
    print(f"{args.year}: wrote {args.output}")


if __name__ == "__main__":
    main()
