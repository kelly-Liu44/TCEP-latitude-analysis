"""Combine yearly TCEP events and compute hemispheric latitude decomposition."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tcep.bootstrap import bootstrap_component_trends
from tcep.decomposition import annual_band_table, decompose_hemisphere


def load_events(directory: Path, start_year: int, end_year: int) -> pd.DataFrame:
    """Read every required year and reject a partially computed event archive."""
    frames = []
    for year in range(start_year, end_year + 1):
        path = directory / f"tcep_events_{year}_pot99.parquet"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_parquet(path)
        if not frame.empty and set(frame["YEAR"].astype(int).unique()) != {year}:
            raise ValueError(f"Wrong YEAR values in {path}")
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    """Write annual band tables, exact components, trends, and 4-year bootstrap CIs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    parser.add_argument("--bootstrap-repetitions", type=int, default=1000)
    parser.add_argument("--bootstrap-block-years", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    events = load_events(args.events_dir, args.start_year, args.end_year)
    years = np.arange(args.start_year, args.end_year + 1)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for hemisphere in ("NH", "SH"):
        annual = annual_band_table(events, hemisphere, years)
        components, trends = decompose_hemisphere(annual)
        intervals, confidence_band = bootstrap_component_trends(
            components,
            repetitions=args.bootstrap_repetitions,
            block_length=args.bootstrap_block_years,
            seed=args.seed,
        )
        intervals["HEMISPHERE"] = hemisphere
        confidence_band["HEMISPHERE"] = hemisphere
        annual.to_csv(args.output_dir / f"annual_bands_{hemisphere}.csv", index=False)
        components.to_csv(
            args.output_dir / f"annual_components_{hemisphere}.csv", index=False
        )
        trends.to_csv(
            args.output_dir / f"component_trends_{hemisphere}.csv", index=False
        )
        intervals.to_csv(
            args.output_dir / f"bootstrap_intervals_{hemisphere}.csv", index=False
        )
        confidence_band.to_csv(
            args.output_dir / f"trend_confidence_band_{hemisphere}.csv", index=False
        )
        print(f"{hemisphere}: decomposition and bootstrap outputs written")


if __name__ == "__main__":
    main()
