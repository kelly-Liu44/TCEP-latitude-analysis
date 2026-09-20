"""Verify that chunked event aggregation matches one-pass aggregation."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / "scripts" / "03_build_tcep_events.py"
spec = importlib.util.spec_from_file_location("event_script", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_event_crossing_two_chunks():
    """Two chunks from the same storm/time become one correctly weighted event."""
    time = pd.Timestamp("2000-01-01", tz="UTC")
    first = pd.DataFrame(
        {
            "SID": ["A"],
            "TIME": [time],
            "BASIN": ["WP"],
            "TC_LAT": [10.0],
            "PRECIP_3H_MM": [10.0],
            "IS_INTERPOLATED_TIMESTEP": [False],
        }
    )
    second = pd.DataFrame(
        {
            "SID": ["A", "A"],
            "TIME": [time, time],
            "BASIN": ["WP", "WP"],
            "TC_LAT": [10.0, 10.0],
            "PRECIP_3H_MM": [20.0, 30.0],
            "IS_INTERPOLATED_TIMESTEP": [True, True],
        }
    )
    events = module.finalize_events(
        [
            module.aggregate_partial_events(first),
            module.aggregate_partial_events(second),
        ]
    )
    assert len(events) == 1
    assert events.loc[0, "N_EXTREME_CELLS"] == 3
    assert np.isclose(events.loc[0, "TCEP_INTENSITY_MM_3H"], 20.0)
    assert bool(events.loc[0, "IS_INTERPOLATED_TIMESTEP"])
