"""Check scientific distinctions that can otherwise silently alter conclusions."""

import numpy as np
import pandas as pd
from tcep.decomposition import (
    annual_band_table,
    decompose_hemisphere,
    theil_sen_per_decade,
)
from tcep.figure_inputs import stats_from


def test_sen_does_not_distribute_over_sum():
    """A median slope is not a linear operator, even for complete annual series."""
    years = np.arange(8)
    first = np.array([2.0, 0.0, 4.0, 3.0, 0.0, 7.0, 1.0, 5.0])
    second = np.array([0.0, 6.0, 1.0, 2.0, 8.0, 0.0, 9.0, 3.0])
    separate = theil_sen_per_decade(first, years) + theil_sen_per_decade(second, years)
    together = theil_sen_per_decade(first + second, years)
    assert not np.isclose(separate, together)


def test_event_average_does_not_weight_by_cell_count():
    """A large low-intensity storm must not dominate the main event-mean estimand."""
    rows = []
    for year in range(2000, 2005):
        for lat, intensity, cells in [
            (5.0, 10.0, 1000),
            (5.0, 30.0, 1),
            (20.0, 20.0, 1),
            (30.0, 20.0, 1),
            (40.0, 20.0, 1),
        ]:
            rows.append(
                {
                    "YEAR": year,
                    "HEMISPHERE": "NH",
                    "TC_LAT": lat,
                    "TCEP_INTENSITY_MM_3H": intensity,
                    "N_EXTREME_CELLS": cells,
                }
            )
    annual = annual_band_table(pd.DataFrame(rows), "NH", np.arange(2000, 2005))
    comp, _ = decompose_hemisphere(annual)
    assert np.allclose(comp.ANNUAL_MEAN_MM_3H, 20.0)


def test_missing_band_year_preserves_annual_identity():
    """Zero frequency and a missing band mean do not remove that year's events."""
    rows = []
    for year in range(2000, 2006):
        for lat in [5.0, 20.0, 30.0, 40.0]:
            if year == 2002 and lat == 20.0:
                continue
            rows.append(
                {
                    "YEAR": year,
                    "HEMISPHERE": "NH",
                    "TC_LAT": lat,
                    "TCEP_INTENSITY_MM_3H": lat + year - 2000,
                }
            )
    events = pd.DataFrame(rows)
    annual = annual_band_table(events, "NH", np.arange(2000, 2006))
    comp, _ = decompose_hemisphere(annual)
    assert np.max(np.abs(comp.CLOSURE_ERROR_MM_3H)) < 1e-10
    assert np.allclose(
        comp.ANNUAL_MEAN_MM_3H, events.groupby("YEAR").TCEP_INTENSITY_MM_3H.mean()
    )


def test_zero_event_basin_band_has_zero_internal_weight():
    """No-event placeholders remain finite internally and are masked for display."""
    counts = np.ones((44, 4))
    counts[:, 3] = 0
    values = counts * 20.0
    result = stats_from(counts, values)
    assert np.allclose(result["mean"], 20.0)
    assert result["mig"][3] == 0 and result["inten"][3] == 0
