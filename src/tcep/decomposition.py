"""Frequency-weighted latitude-band decomposition of TCEP event intensity."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


LATITUDE_BANDS = {
    "NH": ((0.0, 15.0), (15.0, 25.0), (25.0, 35.0), (35.0, 90.0001)),
    "SH": ((0.0, 10.0), (10.0, 15.0), (15.0, 20.0), (20.0, 90.0001)),
}


def band_labels(hemisphere: str) -> list[str]:
    """Return stable text labels for the manuscript latitude bands."""
    return [f"{low:g}-{high:g}" for low, high in LATITUDE_BANDS[hemisphere]]


def assign_band(latitude, hemisphere: str) -> np.ndarray:
    """Assign absolute storm-centre latitude to one fixed hemispheric band."""
    if hemisphere not in LATITUDE_BANDS:
        raise ValueError("hemisphere must be NH or SH")
    absolute = np.abs(np.asarray(latitude, dtype=float))
    result = np.full(absolute.shape, -1, dtype=np.int8)
    for index, (low, high) in enumerate(LATITUDE_BANDS[hemisphere]):
        result[(absolute >= low) & (absolute < high)] = index
    return result


def annual_band_table(
    events: pd.DataFrame, hemisphere: str, years: np.ndarray
) -> pd.DataFrame:
    """Aggregate event counts and event-intensity sums by year and latitude band."""
    required = {"YEAR", "HEMISPHERE", "TC_LAT", "TCEP_INTENSITY_MM_3H"}
    missing = required.difference(events.columns)
    if missing:
        raise ValueError(f"Event table is missing {sorted(missing)}")
    selected = events.loc[events["HEMISPHERE"] == hemisphere].copy()
    selected["BAND"] = assign_band(selected["TC_LAT"], hemisphere)
    if (selected["BAND"] < 0).any():
        raise ValueError(
            f"{hemisphere} contains an event outside its fixed latitude bands"
        )
    grouped = selected.groupby(["YEAR", "BAND"], as_index=False).agg(
        N_EVENTS=("TCEP_INTENSITY_MM_3H", "size"),
        INTENSITY_SUM_MM_3H=("TCEP_INTENSITY_MM_3H", "sum"),
    )
    index = pd.MultiIndex.from_product(
        [years.astype(int), range(len(LATITUDE_BANDS[hemisphere]))],
        names=["YEAR", "BAND"],
    )
    annual = (
        grouped.set_index(["YEAR", "BAND"]).reindex(index, fill_value=0).reset_index()
    )
    annual["MEAN_INTENSITY_MM_3H"] = np.divide(
        annual["INTENSITY_SUM_MM_3H"],
        annual["N_EVENTS"],
        out=np.full(len(annual), np.nan, dtype=float),
        where=annual["N_EVENTS"].to_numpy() > 0,
    )
    annual["BAND_LABEL"] = annual["BAND"].map(dict(enumerate(band_labels(hemisphere))))
    total = annual.groupby("YEAR")["N_EVENTS"].transform("sum")
    annual["RELATIVE_FREQUENCY"] = np.divide(
        annual["N_EVENTS"],
        total,
        out=np.full(len(annual), np.nan, dtype=float),
        where=total.to_numpy() > 0,
    )
    annual["HEMISPHERE"] = hemisphere
    return annual


def theil_sen_per_decade(values, years) -> float:
    """Return a Theil-Sen slope per decade after dropping missing years."""
    values = np.asarray(values, dtype=float)
    years = np.asarray(years, dtype=float)
    valid = np.isfinite(values) & np.isfinite(years)
    if valid.sum() < 3:
        return float("nan")
    return 10.0 * float(stats.theilslopes(values[valid], years[valid]).slope)


def mann_kendall_pvalue(values) -> float:
    """Return the two-sided Mann-Kendall p-value with tie correction."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    count = len(values)
    if count < 3:
        return float("nan")
    first, second = np.triu_indices(count, k=1)
    score = float(np.sign(values[second] - values[first]).sum())
    _, tie_counts = np.unique(values, return_counts=True)
    tie_term = np.sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5))
    variance = (count * (count - 1) * (2 * count + 5) - tie_term) / 18.0
    if variance <= 0:
        return 1.0
    z_score = (score - np.sign(score)) / np.sqrt(variance) if score != 0 else 0.0
    return float(2.0 * stats.norm.sf(abs(z_score)))


def decompose_hemisphere(annual: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute exact annual components and their Theil-Sen trend contributions."""
    hemisphere_values = annual["HEMISPHERE"].unique()
    if len(hemisphere_values) != 1:
        raise ValueError("Annual band table must contain one hemisphere")
    years = np.sort(annual["YEAR"].unique())
    alpha = annual.pivot(
        index="YEAR", columns="BAND", values="RELATIVE_FREQUENCY"
    ).reindex(years)
    intensity = annual.pivot(
        index="YEAR", columns="BAND", values="MEAN_INTENSITY_MM_3H"
    ).reindex(years)
    if alpha.isna().all(axis=1).any():
        raise ValueError("At least one year contains no TCEP events in this hemisphere")
    # A band with no events in one year has zero frequency; its intensity does
    # not contribute that year. Fill its intensity with the band climatology only
    # for the algebraic terms, while alpha remains exactly zero.
    alpha = alpha.fillna(0.0)
    intensity_bar = intensity.mean(axis=0, skipna=True)
    if intensity_bar.isna().any():
        raise ValueError("A latitude band has no TCEP events during the full period")
    intensity_filled = intensity.fillna(intensity_bar)
    alpha_bar = alpha.mean(axis=0)
    alpha_anomaly = alpha - alpha_bar
    intensity_anomaly = intensity_filled - intensity_bar

    components = pd.DataFrame(index=years)
    components.index.name = "YEAR"
    components["ANNUAL_MEAN_MM_3H"] = (alpha * intensity_filled).sum(axis=1)
    components["REFERENCE_MM_3H"] = float((alpha_bar * intensity_bar).sum())
    components["REDISTRIBUTION_MM_3H"] = (alpha_anomaly * intensity_bar).sum(axis=1)
    components["WITHIN_BAND_MM_3H"] = (alpha_bar * intensity_anomaly).sum(axis=1)
    components["INTERACTION_MM_3H"] = (alpha_anomaly * intensity_anomaly).sum(axis=1)
    reconstruction = components[
        [
            "REFERENCE_MM_3H",
            "REDISTRIBUTION_MM_3H",
            "WITHIN_BAND_MM_3H",
            "INTERACTION_MM_3H",
        ]
    ].sum(axis=1)
    components["CLOSURE_ERROR_MM_3H"] = components["ANNUAL_MEAN_MM_3H"] - reconstruction
    components["HEMISPHERE"] = hemisphere_values[0]
    if components["CLOSURE_ERROR_MM_3H"].abs().max() > 1e-10:
        raise RuntimeError("Annual decomposition identity does not close")

    trend_rows = []
    for column, label in [
        ("ANNUAL_MEAN_MM_3H", "TOTAL"),
        ("REDISTRIBUTION_MM_3H", "REDISTRIBUTION"),
        ("WITHIN_BAND_MM_3H", "WITHIN_BAND"),
        ("INTERACTION_MM_3H", "INTERACTION"),
    ]:
        trend_rows.append(
            {
                "HEMISPHERE": hemisphere_values[0],
                "COMPONENT": label,
                "THEIL_SEN_MM_3H_PER_DECADE": theil_sen_per_decade(
                    components[column], years
                ),
                "MANN_KENDALL_P": mann_kendall_pvalue(components[column]),
            }
        )
    trends = pd.DataFrame(trend_rows)
    return components.reset_index(), trends
