"""Audit the completed local wet-period POT99 threshold table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TOTAL_CELLS = 1800 * 3600


def main() -> None:
    """Validate row keys, units, wet counts, and exact tail refinement."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected-exact-cells", type=int, default=92)
    args = parser.parse_args()
    table = pd.read_parquet(args.input)
    required = {
        "ROW",
        "GRID_LAT",
        "GRID_LON",
        "N_WET",
        "POT99_MM_3H",
        "THRESHOLD_METHOD",
        "IS_EXACT_TAIL_REFINEMENT",
        "WET_SAMPLE_RULE",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Threshold table is missing {sorted(missing)}")
    rows = table["ROW"].to_numpy(dtype=np.int64)
    if len(np.unique(rows)) != len(rows):
        raise ValueError("Threshold table contains duplicate grid rows")
    if np.any((rows < 0) | (rows >= TOTAL_CELLS)):
        raise ValueError("Threshold table contains an invalid grid row")
    values = table["POT99_MM_3H"].to_numpy(dtype=np.float64)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.3):
        raise ValueError(
            "POT99 values must be finite and greater than the wet threshold"
        )
    if np.any(table["N_WET"].to_numpy(dtype=np.int64) <= 0):
        raise ValueError("Every stored threshold must have a positive wet-sample count")
    exact = table["IS_EXACT_TAIL_REFINEMENT"].to_numpy(dtype=bool)
    if int(exact.sum()) != args.expected_exact_cells:
        raise ValueError(
            f"Expected {args.expected_exact_cells} exact tail cells; found {int(exact.sum())}"
        )
    if np.any(values[exact] <= 120.0):
        raise ValueError(
            "Every exact tail refinement must resolve a POT99 above 120 mm/3h"
        )
    if np.any(values[~exact] > 120.0):
        raise ValueError("A histogram-derived threshold exceeds its 120 mm/3h ceiling")
    if set(table["WET_SAMPLE_RULE"].astype(str).unique()) != {">0.3 mm/3h"}:
        raise ValueError("Unexpected wet-sample rule")
    print(
        f"PASS: {len(table):,} local POT99 thresholds; "
        f"{int(exact.sum())} exact unbounded tail refinements; "
        f"range={values.min():.6g}..{values.max():.6g} mm/3h"
    )


if __name__ == "__main__":
    main()
