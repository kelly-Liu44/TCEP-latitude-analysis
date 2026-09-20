"""Small deterministic tests that run without the 44-year data archive."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tcep.grid import centroid_lon, grid_row, haversine_km, lon360
from tcep.tcpf import attributed_mask
from tcep.units import accumulation_to_rate, validate_units


def test_units_and_threshold():
    """A 0.3 mm/3h cell is exactly the 0.1 mm/h wet threshold."""
    validate_units("mm/3h")
    assert np.isclose(accumulation_to_rate(0.3), 0.1)
    lat = np.array([[0.0]])
    lon = np.array([[0.0]])
    assert attributed_mask(np.array([[0.31]]), lat, lon, 0.0, 0.0).all()
    assert not attributed_mask(np.array([[0.30]]), lat, lon, 0.0, 0.0).all()
    assert not attributed_mask(np.array([[np.float32(0.30)]]), lat, lon, 0.0, 0.0).all()


def test_four_neighbour_connectivity():
    """Diagonal wet cells are separate features under the manuscript method."""
    accumulation = np.array([[0.31, 0.0], [0.0, 0.31]])
    lat = np.array([[0.05, 0.05], [-0.05, -0.05]])
    lon = np.array([[0.05, 0.15], [0.05, 0.15]])
    selected = attributed_mask(accumulation, lat, lon, 0.0, 0.1)
    assert selected.sum() == 2


def test_dateline_geometry():
    """Dateline distances and centroid longitudes use the shortest arc."""
    assert np.isclose(lon360(-0.1), 359.9)
    assert np.isclose(centroid_lon([359.8, 0.0, 0.2], 0.0) % 360.0, 0.0)
    assert haversine_km(20.0, 359.0, 20.0, 1.0) < 500.0


def test_grid_round_trip():
    """Cell centres map to valid flattened MSWEP row indices."""
    rows = grid_row(np.array([89.95, 0.05, -89.95]), np.array([0.05, 180.05, 359.95]))
    assert rows.tolist() == [0, 3238200, 6479999]
