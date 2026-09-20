"""Verify source hashes and manuscript numerical anchors without raw MSWEP data."""

from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tcep.decomposition import (
    decompose_hemisphere,
    theil_sen_per_decade,
    mann_kendall_pvalue,
)
from tcep.bootstrap import bootstrap_component_trends


def read(path):
    """Preserve North Atlantic's NA code while recognizing missing numeric cells."""
    return pd.read_csv(ROOT / path, keep_default_na=False, na_values=[""])


def verify(bootstrap=False):
    """Raise on any mismatch and return the number of successful checks."""
    checks = []

    def equal(name, actual, expected):
        if not np.allclose(actual, expected, rtol=1e-10, atol=1e-10, equal_nan=True):
            raise AssertionError(name)
        checks.append(name)

    sources = json.loads((ROOT / "metadata/source_data_manifest.json").read_text())
    for item in sources:
        digest = hashlib.sha256((ROOT / item["path"]).read_bytes()).hexdigest()
        if digest != item["sha256"]:
            raise AssertionError("Changed source: " + item["path"])
        checks.append("sha256 " + item["path"])

    expected = json.loads((ROOT / "metadata/expected_results.json").read_text())
    annual = read("figures/Fig01_global_trends/data/Fig01_annual_series.csv")
    for item in expected["global_and_hemisphere_trends"]:
        q = annual[annual.domain == item["domain"]].sort_values("year")
        equal(
            item["domain"] + " slope",
            theil_sen_per_decade(q.TCP_mean, q.year),
            item["trend_per_decade"],
        )
        equal(item["domain"] + " MK", mann_kendall_pvalue(q.TCP_mean), item["p_value"])

    table1 = "tables/TableS01_latitude_decomposition/data/"
    for hemisphere in ["NH", "SH"]:
        bands = read(table1 + f"TableS01_source_annual_bands_{hemisphere}.csv")
        components, trends = decompose_hemisphere(bands)
        saved = read(table1 + f"TableS01_source_component_trends_{hemisphere}.csv")
        equal(hemisphere + " annual closure", components.CLOSURE_ERROR_MM_3H, 0)
        equal(
            hemisphere + " component slopes",
            trends.THEIL_SEN_MM_3H_PER_DECADE,
            saved.THEIL_SEN_MM_3H_PER_DECADE,
        )
        if bootstrap:
            intervals, _ = bootstrap_component_trends(components)
            saved_ci = read(
                table1 + f"TableS01_source_bootstrap_intervals_{hemisphere}.csv"
            )
            equal(
                hemisphere + " 1000-replicate intervals",
                intervals.filter(like="CI_"),
                saved_ci.filter(like="CI_"),
            )
    base = read(table1 + "TableS01C_intensity_baselines.csv")
    equal("event count", int(base.n_events.sum()), expected["events"])
    for row in base.itertuples():
        samples = read(
            f"figures/FigS04_intensity_distributions/data/FigS04_samples_{row.hemisphere}.csv"
        )
        values = samples.loc[samples.band_idx == row.band_idx, "TCEP_INTENSITY_MM_3H"]
        equal(
            f"{row.hemisphere} band {row.band_idx} median",
            values.median(),
            row.pooled_median,
        )
        equal(f"{row.hemisphere} band {row.band_idx} count", len(values), row.n_events)

    spatial = read("figures/Fig01_global_trends/data/Fig01_spatial_trend_1deg.csv")
    equal("map cells", len(spatial), expected["valid_map_cells"])
    equal(
        "significant map cells",
        int((spatial.p_value < 0.05).sum()),
        expected["significant_map_cells"],
    )
    lat = read("figures/FigS03_poleward_shifts/data/FigS03_annual_mean_latitude.csv")
    for item in expected["latitude_trends"]:
        q = lat[
            (lat.record_type == item["record_type"]) & (lat.domain == item["domain"])
        ].sort_values("year")
        equal(
            "latitude slope " + item["record_type"] + item["domain"],
            theil_sen_per_decade(q.mean_lat, q.year),
            item["trend_per_decade"],
        )
        equal(
            "latitude MK " + item["record_type"] + item["domain"],
            mann_kendall_pvalue(q.mean_lat),
            item["p_value"],
        )

    basin = read(
        "tables/TableS02_basin_decomposition/data/TableS02_basin_band_components.csv"
    )
    empty = basin[(basin.domain == "NI") & (basin.band_idx == 3)]
    assert len(empty) == 1 and int(empty.n_events.iloc[0]) == 0
    assert empty.filter(like="_point").isna().all().all()
    checks.append("NI empty band is unavailable")
    surface = read("figures/FigS05_land_ocean/data/FigS05_annual_land_ocean.csv")
    assert (surface.p_value > 0.05).all()
    checks.append("surface total trends are nonsignificant")
    print(f"{len(checks)} checks passed. Bootstrap recalculated: {bootstrap}.")
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bootstrap",
        action="store_true",
        help="Also repeat 1000-replicate aggregate component intervals",
    )
    verify(parser.parse_args().bootstrap)


if __name__ == "__main__":
    main()
