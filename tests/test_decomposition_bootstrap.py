"""Synthetic closure and residual moving-block bootstrap tests."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tcep.bootstrap import bootstrap_component_trends, moving_block_indices
from tcep.decomposition import annual_band_table, assign_band, decompose_hemisphere


def test_asymmetric_latitude_band_boundaries():
    """Exact NH and SH boundaries enter the documented higher-latitude band."""
    assert assign_band([0.0, 14.999, 15.0, 25.0, 35.0], "NH").tolist() == [
        0,
        0,
        1,
        2,
        3,
    ]
    assert assign_band([0.0, -9.999, -10.0, -15.0, -20.0], "SH").tolist() == [
        0,
        0,
        1,
        2,
        3,
    ]


def synthetic_events() -> pd.DataFrame:
    """Build ten years with two populated NH bands and changing intensities."""
    records = []
    for year in range(2000, 2010):
        for event in range(10):
            records.append(
                {
                    "YEAR": year,
                    "HEMISPHERE": "NH",
                    "TC_LAT": 10.0 if event < 6 else 20.0,
                    "TCEP_INTENSITY_MM_3H": 20.0 + 0.5 * (year - 2000) + event,
                }
            )
        # Keep the two higher bands populated so the fixed four-band contract is valid.
        records.append(
            {
                "YEAR": year,
                "HEMISPHERE": "NH",
                "TC_LAT": 30.0,
                "TCEP_INTENSITY_MM_3H": 15.0 + year - 2000,
            }
        )
        records.append(
            {
                "YEAR": year,
                "HEMISPHERE": "NH",
                "TC_LAT": 40.0,
                "TCEP_INTENSITY_MM_3H": 12.0 + year - 2000,
            }
        )
    return pd.DataFrame(records)


def test_exact_decomposition_closure():
    """Annual reference plus three anomalies exactly reconstructs the mean."""
    annual = annual_band_table(synthetic_events(), "NH", np.arange(2000, 2010))
    components, _ = decompose_hemisphere(annual)
    assert components["CLOSURE_ERROR_MM_3H"].abs().max() < 1e-10


def test_four_year_residual_bootstrap_shapes():
    """Four-year moving blocks produce finite intervals and fitted-line bands."""
    annual = annual_band_table(synthetic_events(), "NH", np.arange(2000, 2010))
    components, _ = decompose_hemisphere(annual)
    intervals, band = bootstrap_component_trends(
        components, repetitions=20, block_length=4, seed=42
    )
    assert len(intervals) == 4
    assert len(band) == 10
    assert np.isfinite(intervals.filter(like="CI_").to_numpy()).all()
    indices = moving_block_indices(10, 4, np.random.default_rng(42))
    assert len(indices) == 10
