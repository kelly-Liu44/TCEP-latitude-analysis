"""Create a clean three-hour IBTrACS track table for TCPF extraction."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def interpolate_track(group: pd.DataFrame) -> pd.DataFrame:
    """Interpolate one storm from six-hour synoptic anchors to three hours."""
    group = group.sort_values("TIME").drop_duplicates("TIME").set_index("TIME")
    # Unwrap longitude before interpolation so a 359->1 degree crossing stays local.
    group["LON"] = np.degrees(np.unwrap(np.radians(group["LON"].to_numpy(dtype=float))))
    # Use only complete UTC three-hour ticks inside the observed track lifetime.
    grid = pd.date_range(
        group.index.min().ceil("3h"), group.index.max().floor("3h"), freq="3h", tz="UTC"
    )
    result = group.reindex(grid)
    # Track coordinates must exist at both ends of every interpolated interval.
    result[["LAT", "LON"]] = result[["LAT", "LON"]].interpolate(
        method="time", limit_area="inside"
    )
    # Wind and pressure are optional attributes; interpolate only across valid endpoints.
    result[["USA_WIND", "USA_PRES"]] = result[["USA_WIND", "USA_PRES"]].interpolate(
        method="time", limit_area="inside"
    )
    categorical = ["SID", "NAME", "BASIN", "USA_STATUS", "TRACK_TYPE"]
    for column in categorical:
        if column in result:
            result[column] = result[column].ffill().bfill()
    # Restore the single internal convention after interpolation.
    result["LON"] = result["LON"] % 360.0
    return result.reset_index(names="TIME")


def preserve_native_three_hour_track(group: pd.DataFrame) -> pd.DataFrame:
    """Preserve a valid native 3-hour record when no 6-hour anchor exists."""
    # These rare records cannot be reconstructed from six-hour anchors.  Keeping
    # their observed three-hour positions avoids an implicit storm-strength filter.
    result = group.sort_values("TIME").drop_duplicates("TIME").copy()
    result["LON"] = result["LON"] % 360.0
    return result


def prepare_tracks(
    source: Path, destination: Path, start_year: int, end_year: int
) -> None:
    """Filter, normalize, interpolate, and write the reproducible track table."""
    columns = [
        "SID",
        "SEASON",
        "BASIN",
        "NAME",
        "ISO_TIME",
        "LAT",
        "LON",
        "TRACK_TYPE",
        "USA_STATUS",
        "USA_WIND",
        "USA_PRES",
        "IFLAG",
    ]
    # Preserve the North Atlantic basin code "NA" as text instead of pandas NaN.
    frame = pd.read_csv(
        source,
        skiprows=[1],
        usecols=lambda name: name in columns,
        low_memory=False,
        keep_default_na=False,
    )
    frame["TIME"] = pd.to_datetime(frame.pop("ISO_TIME"), errors="coerce", utc=True)
    for column in ["LAT", "LON", "USA_WIND", "USA_PRES"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame.dropna(subset=["SID", "TIME", "LAT", "LON"])
    # IBTrACS uses both main and spur-* labels; exclude every spur variant.
    frame = frame.loc[~frame["TRACK_TYPE"].fillna("").str.lower().str.contains("spur")]
    frame["LON"] = frame["LON"] % 360.0
    # Keep storms with at least one observation in the analysis period.
    # Adjacent out-of-period points remain temporarily available so interpolation
    # at 1980-01-01 and 2023-12-31 is not biased by truncating a crossing track.
    in_period = frame["TIME"].dt.year.between(start_year, end_year)
    selected_storms = frame.loc[in_period, "SID"].unique()
    frame = frame.loc[frame["SID"].isin(selected_storms)].copy()
    # No wind threshold is applied: weak and wind-missing systems are retained.
    # Reconstruct the manuscript's 6-hour-to-3-hour step from synoptic anchors.
    is_anchor = (
        frame["TIME"].dt.hour.isin([0, 6, 12, 18])
        & frame["TIME"].dt.minute.eq(0)
        & frame["TIME"].dt.second.eq(0)
    )
    anchors = frame.loc[is_anchor].copy()
    pieces = [
        interpolate_track(group) for _, group in anchors.groupby("SID", sort=False)
    ]
    # A handful of non-spur systems in IBTrACS contain only a single 03/09/15/21
    # UTC record and therefore have no six-hour anchor.  Retain those native
    # three-hour observations instead of silently dropping the storms.
    anchored_sids = set(anchors["SID"].unique())
    fallback = frame.loc[~frame["SID"].isin(anchored_sids)].copy()
    fallback = fallback.loc[
        fallback["TIME"].dt.hour.isin([0, 3, 6, 9, 12, 15, 18, 21])
        & fallback["TIME"].dt.minute.eq(0)
        & fallback["TIME"].dt.second.eq(0)
    ]
    pieces.extend(
        preserve_native_three_hour_track(group)
        for _, group in fallback.groupby("SID", sort=False)
    )
    if not pieces:
        raise ValueError("No six-hour IBTrACS anchors were available for interpolation")
    tracks = pd.concat(pieces, ignore_index=True)
    # Apply the analysis period only after each selected storm has been interpolated.
    tracks = tracks.loc[tracks["TIME"].dt.year.between(start_year, end_year)].copy()
    tracks["YEAR"] = tracks["TIME"].dt.year.astype("int16")
    tracks["MONTH"] = tracks["TIME"].dt.month.astype("int8")
    tracks["INTENSITY_CAT"] = pd.cut(
        tracks["USA_WIND"],
        bins=[-np.inf, 33, 63, 82, 95, 112, 136, np.inf],
        labels=["TD", "TS", "C1", "C2", "C3", "C4", "C5"],
    ).astype("string")
    destination.parent.mkdir(parents=True, exist_ok=True)
    tracks.to_parquet(destination, index=False)
    print(
        f"Wrote {len(tracks):,} rows for {tracks['SID'].nunique():,} storms: {destination}"
    )


def main() -> None:
    """Parse command-line paths and prepare tracks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    args = parser.parse_args()
    prepare_tracks(args.input, args.output, args.start_year, args.end_year)


if __name__ == "__main__":
    main()
