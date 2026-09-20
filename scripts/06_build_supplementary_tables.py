"""Build compact manuscript Tables S1-S3 and their full-precision supporting CSVs."""

from pathlib import Path
import argparse
import shutil
import numpy as np
import pandas as pd
from scipy.stats import theilslopes

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "figures"
BAND_COLUMNS = [
    "0–15°N (0–10°S)",
    "15–25°N (10–15°S)",
    "25–35°N (15–20°S)",
    "35–90°N (20–90°S)",
]


def read(folder, name):
    return pd.read_csv(
        HERE / folder / "data" / name, keep_default_na=False, na_values=[""]
    )


def slope(a, t):
    return float(theilslopes(np.asarray(a, float), np.asarray(t, float)).slope * 10)


def number(x):
    """Display two decimals without negative zero; retain missing values as dashes."""
    if pd.isna(x):
        return "—"
    if abs(x) < 0.005:
        x = 0
    return f"{x:.2f}".replace("-", "−")


def band_cell(row):
    """Format redistribution / within-band / residual, using stored intervals."""
    if pd.isna(row["mig_point"]):
        return "—"
    values = []
    for key in ("mig", "int", "resid"):
        significant = row[key + "_ci_low"] > 0 or row[key + "_ci_high"] < 0
        values.append(number(row[key + "_point"]) + (r"\*" if significant else ""))
    return " / ".join(values)


def compact_band_rows(data, domains):
    """Keep hemisphere/basin rows and the four manuscript latitude columns."""
    rows = []
    for code, label in domains:
        subset = data.loc[data.domain == code].set_index("band_idx")
        rows.append([label] + [band_cell(subset.loc[k]) for k in range(4)])
    return rows


def copy_source(source, destination):
    """Allow in-place table rebuilding without copying a source onto itself."""
    if source.resolve() != destination.resolve():
        shutil.copy2(source, destination)


def mdtable(columns, rows):
    return "\n".join(
        [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join(["---"] * len(columns)) + " |",
        ]
        + ["| " + " | ".join(map(str, r)) + " |" for r in rows]
    )


def savecsv(d, folder, name):
    d.to_csv(folder / "data" / name, index=False, float_format="%.17g")


def main():
    global HERE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-root", type=Path, default=ROOT / "figures")
    parser.add_argument("--table-source-root", type=Path, default=ROOT / "tables")
    parser.add_argument(
        "--component-dir",
        type=Path,
        help="Outputs from 04_analyze_decomposition.py; otherwise use included sources",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "tables")
    args = parser.parse_args()
    HERE = args.figure_root
    names = [
        "TableS01_latitude_decomposition",
        "TableS02_basin_decomposition",
        "TableS03_land_ocean",
    ]
    T1, T2, T3 = [args.output_dir / n for n in names]
    for target in [T1, T2, T3]:
        (target / "data").mkdir(parents=True, exist_ok=True)
    for h in ["NH", "SH"]:
        for stem in ["component_trends", "bootstrap_intervals", "annual_bands"]:
            name = f"TableS01_source_{stem}_{h}.csv"
            source = (
                args.component_dir / f"{stem}_{h}.csv"
                if args.component_dir
                else args.table_source_root / names[0] / "data" / name
            )
            copy_source(source, T1 / "data" / name)
    name = "TableS02_event_counts.csv"
    source = (
        (args.figure_root / names[1] / "data" / name)
        if args.component_dir
        else (args.table_source_root / names[1] / "data" / name)
    )
    copy_source(source, T2 / "data" / name)

    allcomp = []
    baseline = []
    for h in ["NH", "SH"]:
        tr = pd.read_csv(T1 / "data" / f"TableS01_source_component_trends_{h}.csv")
        ci = pd.read_csv(
            T1 / "data" / f"TableS01_source_bootstrap_intervals_{h}.csv"
        ).replace({"COMPONENT": {"ANNUAL_MEAN": "TOTAL"}})
        merged = tr.merge(ci, on=["HEMISPHERE", "COMPONENT"], validate="one_to_one")
        merged["bootstrap_significant"] = (merged.CI_LOWER_MM_3H_PER_DECADE > 0) | (
            merged.CI_UPPER_MM_3H_PER_DECADE < 0
        )
        allcomp.append(merged)
        ab = pd.read_csv(T1 / "data" / f"TableS01_source_annual_bands_{h}.csv")
        samples = read("FigS04_intensity_distributions", f"FigS04_samples_{h}.csv")
        for k, q in ab.groupby("BAND"):
            values = samples[samples.band_idx == k].TCEP_INTENSITY_MM_3H
            baseline.append(
                {
                    "hemisphere": h,
                    "band_idx": k,
                    "band_label": q.BAND_LABEL.iloc[0].replace("90.0001", "90"),
                    "climatological_mean": q.MEAN_INTENSITY_MM_3H.mean(),
                    "pooled_median": values.median(),
                    "n_events": len(values),
                }
            )
    comp = pd.concat(allcomp, ignore_index=True)
    base = pd.DataFrame(baseline)
    bands = read(
        "Fig02_latitude_decomposition", "Fig02_decomposition_by_latitude_band.csv"
    )
    savecsv(comp, T1, "TableS01A_hemisphere_components.csv")
    savecsv(bands, T1, "TableS01B_latitude_band_components.csv")
    savecsv(base, T1, "TableS01C_intensity_baselines.csv")
    text = (
        "# Table S1. Latitude-band contributions to hemispheric TCEP intensity trends\n\n"
        + mdtable(
            ["Hemisphere", *BAND_COLUMNS],
            compact_band_rows(
                bands, [("NH", "Northern Hemisphere"), ("SH", "Southern Hemisphere")]
            ),
        )
    )
    text += (
        "\n\nUnits: mm (3 h)⁻¹ decade⁻¹. Each cell lists latitudinal redistribution / "
        "within-latitude-band intensity / residual contributions. Latitude bands outside "
        "and inside parentheses apply to the Northern and Southern Hemispheres, respectively. "
        "An asterisk indicates a 95% residual moving-block bootstrap interval excluding "
        "zero (1,000 repetitions; four-year blocks). Residual = total band trend − "
        "redistribution trend − within-band trend; it includes interaction effects and "
        "Theil–Sen non-additivity and is not the trend of the annual interaction series alone. "
        "Estimates use two decimals; a displayed 0.00 can represent a small nonzero value. "
        "Supporting CSVs retain full precision.\n"
    )
    (T1 / "TableS01.md").write_text(text, encoding="utf-8")
    basins = read(
        "Fig04_basin_decomposition", "Fig04_decomposition_by_basin_and_band.csv"
    )
    counts = pd.read_csv(T2 / "data/TableS02_event_counts.csv", keep_default_na=False)
    basins = basins.merge(
        counts.rename(columns={"BASIN": "domain"}),
        on=["domain", "band_idx"],
        how="left",
        validate="one_to_one",
    )
    basins["n_events"] = basins.n_events.fillna(0).astype(int)
    num = [c for c in basins if c.endswith(("_point", "_ci_low", "_ci_high"))]
    basins.loc[basins.n_events == 0, num] = np.nan
    savecsv(basins, T2, "TableS02_basin_band_components.csv")
    text = (
        "# Table S2. Latitude-band contributions to basin-mean TCEP intensity trends\n\n"
        + mdtable(
            ["Basin", *BAND_COLUMNS],
            compact_band_rows(
                basins,
                [
                    ("WP", "WNP"),
                    ("EP", "EP"),
                    ("NA", "NA"),
                    ("NI", "NIO"),
                    ("SI", "SIO"),
                    ("SP", "SP"),
                ],
            ),
        )
    )
    text += (
        "\n\nUnits: mm (3 h)⁻¹ decade⁻¹. Each cell lists latitudinal redistribution / "
        "within-latitude-band intensity / residual contributions, with residuals and "
        "asterisks defined as in Table S1. A dash denotes a basin–band combination with "
        "no events. WNP, EP, NA, NIO, SIO and SP denote Western Pacific, Eastern Pacific, "
        "North Atlantic, North Indian Ocean, South Indian Ocean and South Pacific, "
        "respectively. Northern bands apply to WNP/EP/NA/NIO; southern bands apply to "
        "SIO/SP. The archived CSV codes WP, NI and SI correspond to WNP, NIO and SIO. "
        "Supporting CSVs retain full precision.\n"
    )
    (T2 / "TableS02.md").write_text(text, encoding="utf-8")
    ab = read("FigS05_land_ocean", "FigS05_annual_surface_band.csv")
    ann = read("FigS05_land_ocean", "FigS05_annual_land_ocean.csv")
    rows = []
    series = []
    for h, s in [("NH", "ocean"), ("SH", "ocean"), ("NH", "land"), ("SH", "land")]:
        q = ab[(ab.domain == h) & (ab.surface == s)]
        n = q.pivot(index="year", columns="band_idx", values="n_j")
        v = q.pivot(index="year", columns="band_idx", values="vol_j")
        alpha = n.div(n.sum(axis=1), axis=0)
        mu = v / n
        mu = mu.fillna(mu.mean())
        abar = alpha.mean()
        mubar = mu.mean()
        da = alpha - abar
        dm = mu - mubar
        total = (alpha * mu).sum(axis=1)
        r = (da * mubar).sum(axis=1)
        w = (abar * dm).sum(axis=1)
        inter = (da * dm).sum(axis=1)
        reference = float((abar * mubar).sum())
        assert np.max(np.abs(total - reference - r - w - inter)) < 1e-10
        p = float(ann[(ann.domain == h) & (ann.surface == s)].p_value.iloc[0])
        rows.append(
            {
                "hemisphere": h,
                "surface": s,
                "total": slope(total, n.index),
                "mk_p": p,
                "redistribution": slope(r, n.index),
                "within_band": slope(w, n.index),
                "interaction": slope(inter, n.index),
            }
        )
        series.append(
            pd.DataFrame(
                {
                    "year": n.index,
                    "hemisphere": h,
                    "surface": s,
                    "total": total.values,
                    "reference": reference,
                    "redistribution": r.values,
                    "within_band": w.values,
                    "interaction": inter.values,
                }
            )
        )
    agg = pd.DataFrame(rows)
    surface = read("FigS05_land_ocean", "FigS05_land_ocean_decomposition.csv")
    savecsv(agg, T3, "TableS03A_surface_components.csv")
    savecsv(surface, T3, "TableS03B_surface_band_components.csv")
    savecsv(pd.concat(series), T3, "TableS03_annual_component_series.csv")
    text = (
        "# Table S3. Land–ocean TCEP intensity trends and decomposition\n\n"
        + mdtable(
            ["Surface and hemisphere", "Total trend", "Decomposition"],
            [
                [
                    f"{x.surface.capitalize()}, {x.hemisphere}",
                    number(x.total),
                    " / ".join(
                        number(v)
                        for v in (x.redistribution, x.within_band, x.interaction)
                    ),
                ]
                for x in agg.itertuples()
            ],
        )
    )
    text += (
        "\n\nUnits: mm (3 h)⁻¹ decade⁻¹. Each decomposition cell lists latitudinal "
        "redistribution / within-latitude-band intensity / interaction contributions. "
        "These are Theil–Sen trends of annually aggregated components; their slopes "
        "are not additive. Here interaction denotes the annual interaction-series "
        "trend, distinct from the band residual in Tables S1 and S2. This analysis "
        "weights extreme-grid-cell records equally, whereas the main analysis weights "
        "storm–3-hour events equally. Surface type follows precipitation-grid location; "
        "hemisphere and latitude band follow storm-center latitude. Total trends are "
        "not significant at p < 0.05; their Mann–Kendall p values are retained in the "
        "supporting CSV. Component confidence intervals are not estimated, so component "
        "significance is not assigned. Estimates use two decimals.\n"
    )
    (T3 / "TableS03.md").write_text(text, encoding="utf-8")
    print(
        "Built compact Tables S1–S3; supporting CSVs retain computational precision; display estimates use two decimals."
    )


if __name__ == "__main__":
    main()
