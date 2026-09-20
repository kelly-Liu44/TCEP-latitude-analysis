"""Connected-component TCPF attribution for one storm and one time step."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import label

from .grid import centroid_lon, haversine_km


@dataclass(frozen=True)
class Parameters:
    """TCPF criteria in the units used by the manuscript."""

    wet_threshold_mm_3h: float = 0.3
    search_radius_km: float = 1000.0
    centroid_radius_km: float = 500.0


FOUR_NEIGHBOUR = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8)


def attributed_mask(
    accumulation_mm_3h, lat_grid, lon_grid, tc_lat, tc_lon, parameters=Parameters()
):
    """Return cells belonging to a feature whose centroid is within 500 km."""
    accumulation = np.asarray(accumulation_mm_3h)
    if not np.issubdtype(accumulation.dtype, np.floating):
        accumulation = accumulation.astype(float)
    lat_grid, lon_grid = (
        np.asarray(lat_grid, dtype=float),
        np.asarray(lon_grid, dtype=float),
    )
    if accumulation.shape != lat_grid.shape or accumulation.shape != lon_grid.shape:
        raise ValueError("Rainfall and coordinate grids must have identical shapes")
    distance = haversine_km(tc_lat, tc_lon, lat_grid, lon_grid)
    # Compare native accumulated depth directly to avoid /3 floating-point drift.
    native_threshold = np.asarray(
        parameters.wet_threshold_mm_3h, dtype=accumulation.dtype
    )
    wet = (
        np.isfinite(accumulation)
        & (accumulation > native_threshold)
        & (distance <= parameters.search_radius_km)
    )
    feature_ids, count = label(wet, structure=FOUR_NEIGHBOUR)
    selected = np.zeros(wet.shape, dtype=bool)
    for feature_id in range(1, count + 1):
        feature = feature_ids == feature_id
        centre_lat = float(np.mean(lat_grid[feature]))
        centre_lon = centroid_lon(lon_grid[feature], tc_lon)
        centre_distance = float(haversine_km(tc_lat, tc_lon, centre_lat, centre_lon))
        if centre_distance <= parameters.centroid_radius_km:
            selected[feature] = True
    return selected
