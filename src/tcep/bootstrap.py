"""Residual moving-block bootstrap for annual TCEP decomposition components."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


COMPONENT_COLUMNS = [
    "ANNUAL_MEAN_MM_3H",
    "REDISTRIBUTION_MM_3H",
    "WITHIN_BAND_MM_3H",
    "INTERACTION_MM_3H",
]


def moving_block_indices(
    length: int, block_length: int, rng: np.random.Generator
) -> np.ndarray:
    """Sample overlapping residual blocks until the original series length is met."""
    if not 1 <= block_length <= length:
        raise ValueError("block_length must be between 1 and series length")
    number_of_blocks = int(np.ceil(length / block_length))
    starts = rng.integers(0, length - block_length + 1, size=number_of_blocks)
    return np.concatenate([np.arange(start, start + block_length) for start in starts])[
        :length
    ]


def sen_line(values: np.ndarray, years: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the fitted Theil-Sen line and residuals for one annual series."""
    result = stats.theilslopes(values, years)
    fitted = result.intercept + result.slope * years
    return fitted, values - fitted


def bootstrap_component_trends(
    components: pd.DataFrame,
    repetitions: int = 1000,
    block_length: int = 4,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bootstrap component slopes and the annual-mean fitted-line confidence band."""
    ordered = components.sort_values("YEAR").reset_index(drop=True)
    years = ordered["YEAR"].to_numpy(dtype=float)
    values = ordered[COMPONENT_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(
            "Bootstrap component series contain missing or infinite values"
        )
    fitted = np.empty_like(values)
    residuals = np.empty_like(values)
    for index in range(values.shape[1]):
        fitted[:, index], residuals[:, index] = sen_line(values[:, index], years)

    rng = np.random.default_rng(seed)
    slopes = np.empty((repetitions, len(COMPONENT_COLUMNS)), dtype=float)
    fitted_total = np.empty((repetitions, len(years)), dtype=float)
    for repetition in range(repetitions):
        indices = moving_block_indices(len(years), block_length, rng)
        reconstructed = fitted + residuals[indices, :]
        for column_index in range(reconstructed.shape[1]):
            result = stats.theilslopes(reconstructed[:, column_index], years)
            slopes[repetition, column_index] = 10.0 * float(result.slope)
            if column_index == 0:
                fitted_total[repetition] = result.intercept + result.slope * years

    interval = pd.DataFrame(
        {
            "COMPONENT": [name.replace("_MM_3H", "") for name in COMPONENT_COLUMNS],
            "CI_LOWER_MM_3H_PER_DECADE": np.quantile(slopes, 0.025, axis=0),
            "CI_UPPER_MM_3H_PER_DECADE": np.quantile(slopes, 0.975, axis=0),
        }
    )
    central_total, _ = sen_line(values[:, 0], years)
    confidence_band = pd.DataFrame(
        {
            "YEAR": years.astype(int),
            "FIT_MM_3H": central_total,
            "FIT_LOWER_MM_3H": np.quantile(fitted_total, 0.025, axis=0),
            "FIT_UPPER_MM_3H": np.quantile(fitted_total, 0.975, axis=0),
        }
    )
    return interval, confidence_band
