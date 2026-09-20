"""Local-POT extreme-cell selection and storm/time-step aggregation."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .grid import NLAT, NLON, grid_row


def threshold_lookup(table: pd.DataFrame, column: str = "POT99_MM_3H") -> np.ndarray:
    """Expand one threshold per global grid row into a fast array lookup."""
    required = {"ROW", column}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Threshold table is missing {sorted(missing)}")
    rows = table["ROW"].to_numpy(dtype=np.int64)
    if len(np.unique(rows)) != len(rows) or np.any((rows < 0) | (rows >= NLAT * NLON)):
        raise ValueError("Threshold rows must be unique and within the global grid")
    lookup = np.full(NLAT * NLON, np.nan, dtype=np.float32)
    lookup[rows] = table[column].to_numpy(dtype=np.float32)
    return lookup


def select_extreme_cells(cells: pd.DataFrame, thresholds: np.ndarray) -> pd.DataFrame:
    """Return TCPF cells strictly exceeding their local wet-period POT99."""
    required = {
        "SID",
        "TIME",
        "BASIN",
        "TC_LAT",
        "GRID_LAT",
        "GRID_LON",
        "PRECIP_3H_MM",
    }
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"TCPF cells are missing {sorted(missing)}")
    work = cells.copy()
    work["TIME"] = pd.to_datetime(work["TIME"], utc=True)
    # Prefer the stored canonical row key; recompute only for compatible inputs.
    rows = (
        work["ROW"].to_numpy(dtype=np.int64)
        if "ROW" in work
        else grid_row(work["GRID_LAT"].to_numpy(), work["GRID_LON"].to_numpy())
    )
    local_threshold = thresholds[rows]
    selected = work.loc[
        np.isfinite(local_threshold)
        & (work["PRECIP_3H_MM"].to_numpy() > local_threshold)
    ].copy()
    selected["ROW"] = rows[
        np.isfinite(local_threshold)
        & (work["PRECIP_3H_MM"].to_numpy() > local_threshold)
    ]
    selected["POT99_MM_3H"] = local_threshold[
        np.isfinite(local_threshold)
        & (work["PRECIP_3H_MM"].to_numpy() > local_threshold)
    ]
    return selected


def build_events(cells: pd.DataFrame, thresholds: np.ndarray) -> pd.DataFrame:
    """Average strict local-POT exceedances by storm and three-hour time step."""
    selected = select_extreme_cells(cells, thresholds)
    if selected.empty:
        return pd.DataFrame(
            columns=[
                "SID",
                "TIME",
                "YEAR",
                "HEMISPHERE",
                "BASIN",
                "TC_LAT",
                "TCEP_INTENSITY_MM_3H",
                "N_EXTREME_CELLS",
                "IS_INTERPOLATED_TIMESTEP",
            ]
        )
    if "IS_INTERPOLATED_TIMESTEP" not in selected:
        selected["IS_INTERPOLATED_TIMESTEP"] = False
    events = selected.groupby(["SID", "TIME"], as_index=False).agg(
        BASIN=("BASIN", "first"),
        TC_LAT=("TC_LAT", "first"),
        TCEP_INTENSITY_MM_3H=("PRECIP_3H_MM", "mean"),
        N_EXTREME_CELLS=("PRECIP_3H_MM", "size"),
        IS_INTERPOLATED_TIMESTEP=("IS_INTERPOLATED_TIMESTEP", "max"),
    )
    events["YEAR"] = events["TIME"].dt.year.astype("int16")
    events["HEMISPHERE"] = np.where(events["TC_LAT"] >= 0, "NH", "SH")
    return events[
        [
            "SID",
            "TIME",
            "YEAR",
            "HEMISPHERE",
            "BASIN",
            "TC_LAT",
            "TCEP_INTENSITY_MM_3H",
            "N_EXTREME_CELLS",
            "IS_INTERPOLATED_TIMESTEP",
        ]
    ]
