"""Verify strict POT99 selection and storm/time-step event aggregation."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tcep.events import build_events
from tcep.grid import NLAT, NLON, grid_row


def test_strict_threshold_and_event_mean():
    """Equality is excluded; two strict exceedances form one storm-time event."""
    time = pd.Timestamp("2000-01-01 00:00", tz="UTC")
    cells = pd.DataFrame(
        {
            "SID": ["A", "A", "A"],
            "TIME": [time, time, time],
            "BASIN": ["WP", "WP", "WP"],
            "TC_LAT": [10.0, 10.0, 10.0],
            "GRID_LAT": [0.05, 0.05, 0.05],
            "GRID_LON": [0.05, 0.15, 0.25],
            "PRECIP_3H_MM": [10.0, 12.0, 14.0],
        }
    )
    lookup = np.full(NLAT * NLON, np.nan, dtype=np.float32)
    lookup[grid_row(cells.GRID_LAT, cells.GRID_LON)] = [10.0, 11.0, 13.0]
    events = build_events(cells, lookup)
    assert len(events) == 1
    assert events.loc[0, "N_EXTREME_CELLS"] == 2
    assert np.isclose(events.loc[0, "TCEP_INTENSITY_MM_3H"], 13.0)
