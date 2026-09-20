"""Conversions and validation for the native three-hour MSWEP quantity."""

from __future__ import annotations

import numpy as np

ACCUMULATION_HOURS = 3.0
ACCEPTED_UNITS = {"mm/3h", "mm/3hr", "mm/3hour", "mm 3h-1", "mm3h-1"}


def validate_units(units: object) -> None:
    """Reject a source variable unless it declares a three-hour accumulation."""
    normalized = str(units or "").lower().replace(" ", "")
    accepted = {item.replace(" ", "") for item in ACCEPTED_UNITS}
    if normalized not in accepted:
        raise ValueError(f"Expected MSWEP units equivalent to mm/3h, got {units!r}")


def accumulation_to_rate(values_mm_3h) -> np.ndarray:
    """Convert a three-hour accumulation in millimetres to mean mm h-1."""
    return np.asarray(values_mm_3h, dtype=float) / ACCUMULATION_HOURS
