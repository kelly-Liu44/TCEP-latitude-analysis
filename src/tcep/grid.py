"""Dateline-safe spherical geometry and MSWEP grid indexing."""

from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0
GRID_DEG = 0.1
NLAT, NLON = 1800, 3600


def lon360(values):
    """Return longitudes in the single internal [0, 360) convention."""
    return np.mod(np.asarray(values, dtype=float), 360.0)


def wrapped_delta(longitudes, reference):
    """Return the shortest signed longitude difference in degrees."""
    return (np.asarray(longitudes) - reference + 180.0) % 360.0 - 180.0


def haversine_km(lat1, lon1, lat2, lon2):
    """Calculate great-circle distance while correctly crossing the date line."""
    lat1_rad, lat2_rad = np.radians(lat1), np.radians(lat2)
    dlat_rad = lat2_rad - lat1_rad
    dlon_rad = np.radians(wrapped_delta(lon2, lon1))
    a = (
        np.sin(dlat_rad / 2.0) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon_rad / 2.0) ** 2
    )
    return 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def centroid_lon(longitudes, reference):
    """Average a feature longitude after unwrapping it around the TC centre."""
    return float(
        (reference + np.mean(wrapped_delta(lon360(longitudes), reference))) % 360.0
    )


def grid_row(latitudes, longitudes):
    """Map 0.1-degree MSWEP cell centres to flattened row indices."""
    lat_index = np.rint((89.95 - np.asarray(latitudes)) / GRID_DEG).astype(np.int64)
    lon_index = np.rint((lon360(longitudes) - 0.05) / GRID_DEG).astype(np.int64)
    if np.any(
        (lat_index < 0) | (lat_index >= NLAT) | (lon_index < 0) | (lon_index >= NLON)
    ):
        raise ValueError("Coordinates are outside the 0.1-degree MSWEP grid")
    return lat_index * NLON + lon_index
