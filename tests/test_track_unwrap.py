"""Verify that dateline-crossing tracks interpolate along the short arc."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / "scripts" / "01_prepare_ibtracs.py"
spec = importlib.util.spec_from_file_location("prepare", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

times = pd.date_range("2000-01-01", periods=2, freq="6h", tz="UTC")
frame = pd.DataFrame(
    {
        "TIME": times,
        "LON": [359.0, 1.0],
        "LAT": [10.0, 10.0],
        "USA_WIND": [40.0, 40.0],
        "USA_PRES": [990.0, 990.0],
        "SID": ["x", "x"],
        "NAME": ["x", "x"],
        "BASIN": ["WP", "WP"],
        "USA_STATUS": ["TS", "TS"],
        "TRACK_TYPE": ["main", "main"],
    }
)
result = module.interpolate_track(frame)
assert np.isclose(
    result.loc[result["TIME"] == times[0] + pd.Timedelta(hours=3), "LON"].iloc[0], 0.0
)
