"""Prepare manuscript figure inputs from saved storm events, extreme cells and tracks."""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
EVENT_DIR = ROOT / "outputs" / "events"
CELL_DIR = ROOT / "outputs" / "extreme_cells"
TRACK_FILE = ROOT / "data" / "tracks_processed.parquet"
OUT = ROOT / "results" / "figure_inputs"
EVENT_SUFFIX = "pot99"
YEARS = np.arange(1980, 2024)
B = 1000
BLOCK = 4
SEED = 42
BANDS = {
    "NH": [(0, 15), (15, 25), (25, 35), (35, 90.0001)],
    "SH": [(0, 10), (10, 15), (15, 20), (20, 90.0001)],
}
BASINS = [
    ("WP", "NH"),
    ("EP", "NH"),
    ("NA", "NH"),
    ("NI", "NH"),
    ("SI", "SH"),
    ("SP", "SH"),
]


def folder(name):
    """Create one numbered output folder and its data directory."""
    p = OUT / name
    (p / "data").mkdir(parents=True, exist_ok=True)
    return p


def events():
    """Read the complete annual POT99 event archive; preserve the NA basin code."""
    d = pd.concat(
        [
            pd.read_parquet(EVENT_DIR / f"tcep_events_{y}_{EVENT_SUFFIX}.parquet")
            for y in YEARS
        ],
        ignore_index=True,
    )
    d["YEAR"] = d.YEAR.astype(int)
    d["BASIN"] = d.BASIN.replace({"NATL": "NA"})
    return d


def band_index(lat, hemi):
    """Assign absolute storm-center latitude to manuscript latitude bands."""
    a = np.abs(np.asarray(lat, float))
    out = np.full(len(a), -1, dtype=np.int8)
    for k, (lo, hi) in enumerate(BANDS[hemi]):
        out[(a >= lo) & (a < hi)] = k
    return out


def labels(hemi):
    """Return display labels without the numerical inclusive-pole sentinel."""
    return [f"{lo:g}-{hi:g}".replace("90.0001", "90") for lo, hi in BANDS[hemi]]


def sen(y, t=None):
    """Estimate a per-year Theil–Sen slope from finite observations."""
    if t is None:
        t = YEARS
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    m = np.isfinite(y) & np.isfinite(t)
    return np.nan if m.sum() < 3 else float(stats.theilslopes(y[m], t[m]).slope)


def mk(y):
    """Return the continuity- and tie-corrected two-sided Mann–Kendall p value."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    n = len(y)
    if n < 3:
        return np.nan
    i, j = np.triu_indices(n, 1)
    s = np.sign(y[j] - y[i]).sum()
    _, c = np.unique(y, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(c * (c - 1) * (2 * c + 5))) / 18
    z = (s - np.sign(s)) / np.sqrt(var) if s and var > 0 else 0
    return float(2 * stats.norm.sf(abs(z)))


def sen_line(y, t=None):
    """Fit a Sen slope with the joint median-residual intercept."""
    if t is None:
        t = YEARS
    sl = sen(y, t)
    sl = 0 if not np.isfinite(sl) else sl
    y = np.asarray(y, float)
    t = np.asarray(t, float)
    m = np.isfinite(y)
    ic = float(np.median(y[m] - sl * t[m]))
    fit = ic + sl * t
    return fit, y - fit


def block_idx(rng, n=None):
    """Draw overlapping residual blocks using the supplied reproducible generator."""
    if n is None:
        n = len(YEARS)
    starts = rng.integers(0, n - BLOCK + 1, size=int(np.ceil(n / BLOCK)))
    return np.concatenate([np.arange(s, s + BLOCK) for s in starts])[:n]


def annual_band(frame, hemi):
    """Tabulate event counts and intensity sums, retaining empty year-band cells."""
    f = frame.copy()
    f["band_idx"] = band_index(f.TC_LAT, hemi)
    g = (
        f.groupby(["YEAR", "band_idx"])
        .TCEP_INTENSITY_MM_3H.agg(["size", "sum"])
        .rename(columns={"size": "n_j", "sum": "vol_j"})
    )
    idx = pd.MultiIndex.from_product([YEARS, range(4)], names=["YEAR", "band_idx"])
    return g.reindex(idx, fill_value=0).reset_index()


def matrices(ab):
    """Convert a complete annual-band table to aligned year-by-band arrays."""
    n = (
        ab.pivot(index="YEAR", columns="band_idx", values="n_j")
        .reindex(YEARS, fill_value=0)
        .to_numpy(float)
    )
    v = (
        ab.pivot(index="YEAR", columns="band_idx", values="vol_j")
        .reindex(YEARS, fill_value=0)
        .to_numpy(float)
    )
    return n, v


def stats_from(n, v):
    """Calculate band slopes; residual includes interaction and Sen non-additivity."""
    nt = n.sum(1)
    alpha = np.divide(n, nt[:, None], out=np.zeros_like(n), where=nt[:, None] > 0)
    p = np.divide(v, n, out=np.full_like(v, np.nan), where=n > 0)
    q = np.divide(v, nt[:, None], out=np.zeros_like(v), where=nt[:, None] > 0)
    abar = alpha.mean(0)
    valid_p = np.isfinite(p)
    pbar = np.divide(
        np.nansum(p, axis=0),
        valid_p.sum(axis=0),
        out=np.zeros(p.shape[1], dtype=float),
        where=valid_p.sum(axis=0) > 0,
    )
    mig = np.array(
        [
            (sen(alpha[:, k]) * 10 if np.isfinite(sen(alpha[:, k])) else 0) * pbar[k]
            for k in range(n.shape[1])
        ]
    )
    inte = np.array(
        [
            abar[k] * (sen(p[:, k]) * 10 if np.isfinite(sen(p[:, k])) else 0)
            for k in range(n.shape[1])
        ]
    )
    total = np.array(
        [
            (sen(q[:, k]) * 10 if np.isfinite(sen(q[:, k])) else 0)
            for k in range(n.shape[1])
        ]
    )
    resid = total - mig - inte
    mean = q.sum(1)
    return dict(
        mig=mig,
        inten=inte,
        total=total,
        resid=resid,
        mean=mean,
        slope=sen(mean) * 10,
        p=mk(mean),
    )


def bootstrap_decomp(ab, hemi, domain):
    """Jointly resample detrended counts and intensity sums for band intervals."""
    n, v = matrices(ab)
    pt = stats_from(n, v)
    rng = np.random.default_rng(SEED)
    nf = np.empty_like(n)
    nr = np.empty_like(n)
    vf = np.empty_like(v)
    vr = np.empty_like(v)
    for k in range(4):
        nf[:, k], nr[:, k] = sen_line(n[:, k])
        vf[:, k], vr[:, k] = sen_line(v[:, k])
    tm = np.empty((B, 4))
    ti = np.empty((B, 4))
    tr = np.empty((B, 4))
    tt = np.empty((B, 4))
    bf = np.empty((B, len(YEARS)))
    meanfit, meanres = sen_line(pt["mean"])
    for b in range(B):
        ix = block_idx(rng)
        st = stats_from(np.clip(nf + nr[ix], 0, None), np.clip(vf + vr[ix], 0, None))
        tm[b] = st["mig"]
        ti[b] = st["inten"]
        tr[b] = st["resid"]
        tt[b] = st["total"]
        star = meanfit + meanres[ix]
        bf[b], _ = sen_line(star)
    rows = []
    for k, label in enumerate(labels(hemi)):
        row = {"domain": domain, "hemisphere": hemi, "band_idx": k, "band_label": label}
        for key, arr, boot in [
            ("mig", pt["mig"], tm),
            ("int", pt["inten"], ti),
            ("resid", pt["resid"], tr),
            ("tot", pt["total"], tt),
        ]:
            row[f"{key}_point"] = arr[k]
            row[f"{key}_ci_low"], row[f"{key}_ci_high"] = np.quantile(
                boot[:, k], [0.025, 0.975]
            )
        rows.append(row)
    annual = pd.DataFrame(
        {
            "domain": domain,
            "year": YEARS,
            "N_t": n.sum(1),
            "TCP_mean": pt["mean"],
            "fit": meanfit,
            "fit_lower": np.quantile(bf, 0.025, 0),
            "fit_upper": np.quantile(bf, 0.975, 0),
            "trend_per_decade": pt["slope"],
            "p_value": pt["p"],
        }
    )
    return pd.DataFrame(rows), annual


def prepare_decomposition(ev):
    """Export hemispheric and basin band contributions and basin mean series."""
    f2 = folder("Fig02_latitude_decomposition")
    all_rows = []
    for h in ["NH", "SH"]:
        rows, annual = bootstrap_decomp(annual_band(ev[ev.HEMISPHERE == h], h), h, h)
        all_rows.append(rows)
        annual.to_csv(f2 / "data" / f"Fig02_annual_{h}.csv", index=False)
    pd.concat(all_rows).to_csv(
        f2 / "data" / "Fig02_decomposition_by_latitude_band.csv", index=False
    )
    f3 = folder("Fig04_basin_decomposition")
    rows = []
    annuals = []
    for basin, h in BASINS:
        sub = ev[ev.BASIN == basin]
        r, a = bootstrap_decomp(annual_band(sub, h), h, basin)
        rows.append(r)
        annuals.append(a)
    pd.concat(rows).to_csv(
        f3 / "data" / "Fig04_decomposition_by_basin_and_band.csv", index=False
    )
    pd.concat(annuals).to_csv(
        f3 / "data" / "Fig04_annual_basin_series.csv", index=False
    )


def prepare_fig1(ev):
    """Export equal-event annual means and separately aggregated fixed-grid trends."""
    f = folder("Fig01_global_trends")
    annual = []
    for domain, sub in [
        ("GLOBAL", ev),
        ("NH", ev[ev.HEMISPHERE == "NH"]),
        ("SH", ev[ev.HEMISPHERE == "SH"]),
    ]:
        a = (
            sub.groupby("YEAR")
            .TCEP_INTENSITY_MM_3H.agg(["size", "mean"])
            .reindex(YEARS)
        )
        fit, res = sen_line(a["mean"])
        rng = np.random.default_rng(SEED)
        bfit = np.array([sen_line(fit + res[block_idx(rng)])[0] for _ in range(B)])
        annual.append(
            pd.DataFrame(
                {
                    "domain": domain,
                    "year": YEARS,
                    "N_t": a["size"].values,
                    "TCP_mean": a["mean"].values,
                    "fit": fit,
                    "fit_lower": np.quantile(bfit, 0.025, 0),
                    "fit_upper": np.quantile(bfit, 0.975, 0),
                    "trend_per_decade": sen(a["mean"]) * 10,
                    "p_value": mk(a["mean"]),
                }
            )
        )
    pd.concat(annual).to_csv(f / "data" / "Fig01_annual_series.csv", index=False)
    # One-degree spatial annual means from extreme cells, streamed year by year.
    pieces = []
    for y in YEARS:
        d = pd.read_parquet(
            CELL_DIR / f"tcep_extreme_cells_{y}_{EVENT_SUFFIX}.parquet",
            columns=["GRID_LAT", "GRID_LON", "PRECIP_3H_MM"],
        )
        d["lat"] = np.floor(d.GRID_LAT).astype(np.int16)
        d["lon"] = np.floor(d.GRID_LON % 360).astype(np.int16)
        g = d.groupby(["lat", "lon"]).PRECIP_3H_MM.agg(["sum", "count"]).reset_index()
        g["year"] = y
        pieces.append(g)
    ag = pd.concat(pieces)
    pivot_sum = ag.pivot_table(index=["lat", "lon"], columns="year", values="sum")
    pivot_n = ag.pivot_table(index=["lat", "lon"], columns="year", values="count")
    vals = pivot_sum / pivot_n
    out = []
    for (lat, lon), row in vals.iterrows():
        if row.notna().sum() >= 10:
            out.append(
                {
                    "lat": lat + 0.5,
                    "lon": lon + 0.5,
                    "trend": sen(row.reindex(YEARS).values) * 10,
                    "p_value": mk(row.reindex(YEARS).values),
                    "valid_years": row.notna().sum(),
                }
            )
    pd.DataFrame(out).to_csv(f / "data" / "Fig01_spatial_trend_1deg.csv", index=False)


def prepare_supp(ev):
    """Export full conditional intensity samples and absolute-latitude trends."""
    f = folder("FigS04_intensity_distributions")
    for h in ["NH", "SH"]:
        sub = ev[ev.HEMISPHERE == h].copy()
        sub["band_idx"] = band_index(sub.TC_LAT, h)
        sample = sub[["band_idx", "TCEP_INTENSITY_MM_3H"]].copy()
        sample["band_label"] = sample.band_idx.map(dict(enumerate(labels(h))))
        sample.to_csv(f / "data" / f"FigS04_samples_{h}.csv", index=False)
        sub.groupby("band_idx").TCEP_INTENSITY_MM_3H.quantile(
            [0.1, 0.5, 0.9, 0.99]
        ).unstack().reset_index().to_csv(
            f / "data" / f"FigS04_quantiles_{h}.csv", index=False
        )
    f = folder("FigS03_poleward_shifts")
    rows = []
    tr = pd.read_parquet(TRACK_FILE, columns=["YEAR", "LAT"])
    tr["HEMISPHERE"] = np.where(tr.LAT >= 0, "NH", "SH")
    tr["ABS_LAT"] = tr.LAT.abs()
    for kind, frame, col in [
        ("track", tr, "ABS_LAT"),
        ("event", ev.assign(ABS_LAT=ev.TC_LAT.abs()), "ABS_LAT"),
    ]:
        for h in ["NH", "SH"]:
            a = (
                frame[frame.HEMISPHERE == h]
                .groupby("YEAR")[col]
                .agg(["mean", "size"])
                .reindex(YEARS)
            )
            fit, res = sen_line(a["mean"])
            rng = np.random.default_rng(SEED)
            bf = np.array([sen_line(fit + res[block_idx(rng)])[0] for _ in range(B)])
            rows.append(
                pd.DataFrame(
                    {
                        "record_type": kind,
                        "domain": h,
                        "year": YEARS,
                        "mean_lat": a["mean"].values,
                        "n_records": a["size"].values,
                        "fit": fit,
                        "fit_lower": np.quantile(bf, 0.025, 0),
                        "fit_upper": np.quantile(bf, 0.975, 0),
                        "trend_per_decade": sen(a["mean"]) * 10,
                        "p_value": mk(a["mean"]),
                    }
                )
            )
    pd.concat(rows).to_csv(f / "data" / "FigS03_annual_mean_latitude.csv", index=False)
    f = folder("Fig03_basin_intensity_trends")
    src = folder("Fig04_basin_decomposition") / "data" / "Fig04_annual_basin_series.csv"
    if src.exists():
        pd.read_csv(src, keep_default_na=False, na_values=[""]).to_csv(
            f / "data" / "Fig03_annual_basin_series.csv", index=False
        )


def prepare_surface_cells():
    """Classify precipitation cells by surface and retain archived CSV precision."""
    f = folder("FigS05_land_ocean")
    from global_land_mask import globe

    rows = []
    for y in YEARS:
        d = pd.read_parquet(
            CELL_DIR / f"tcep_extreme_cells_{y}_{EVENT_SUFFIX}.parquet",
            columns=["TC_LAT", "GRID_LAT", "GRID_LON", "PRECIP_3H_MM"],
        )
        lon = ((d.GRID_LON.to_numpy(float) + 180) % 360) - 180
        land = globe.is_land(
            np.clip(d.GRID_LAT.to_numpy(float), -89.999, 89.999),
            np.clip(lon, -179.999, 179.999),
        )
        d["surface"] = np.where(land, "land", "ocean")
        d["domain"] = np.where(d.TC_LAT >= 0, "NH", "SH")
        for (h, s), q in d.groupby(["domain", "surface"]):
            bi = band_index(q.TC_LAT, h)
            q = q.assign(band_idx=bi)
            g = (
                q.groupby("band_idx")
                .PRECIP_3H_MM.agg(["size", "sum"])
                .reindex(range(4), fill_value=0)
            )
            for k in range(4):
                rows.append(
                    {
                        "year": y,
                        "domain": h,
                        "surface": s,
                        "band_idx": k,
                        "band_label": labels(h)[k],
                        "n_j": g.loc[k, "size"],
                        "vol_j": g.loc[k, "sum"],
                    }
                )
    ab = pd.DataFrame(rows)
    ab.to_csv(f / "data" / "FigS05_annual_surface_band.csv", index=False)
    # Match the archived analysis: serialized float32 sums are read as float64
    # before deriving annual trends and latitude-band contributions.
    ab = pd.read_csv(f / "data" / "FigS05_annual_surface_band.csv")
    annual = []
    decomp = []
    for h in ["NH", "SH"]:
        for s in ["land", "ocean"]:
            x = ab[(ab.domain == h) & (ab.surface == s)]
            n, v = matrices(x.rename(columns={"year": "YEAR"}))
            st = stats_from(n, v)
            fit, res = sen_line(st["mean"])
            rng = np.random.default_rng(SEED)
            bf = np.array([sen_line(fit + res[block_idx(rng)])[0] for _ in range(B)])
            annual.append(
                pd.DataFrame(
                    {
                        "domain": h,
                        "surface": s,
                        "year": YEARS,
                        "N_t": n.sum(1),
                        "TCP_mean": st["mean"],
                        "fit": fit,
                        "fit_lower": np.quantile(bf, 0.025, 0),
                        "fit_upper": np.quantile(bf, 0.975, 0),
                        "trend_per_decade": st["slope"],
                        "p_value": st["p"],
                    }
                )
            )
            for k, label in enumerate(labels(h)):
                decomp.append(
                    {
                        "domain": h,
                        "surface": s,
                        "band_idx": k,
                        "band_label": label,
                        "migration": st["mig"][k],
                        "intensity": st["inten"][k],
                        "residual": st["resid"][k],
                        "total": st["total"][k],
                    }
                )
    pd.concat(annual).to_csv(f / "data" / "FigS05_annual_land_ocean.csv", index=False)
    pd.DataFrame(decomp).to_csv(
        f / "data" / "FigS05_land_ocean_decomposition.csv", index=False
    )


def main():
    """Dispatch reproducible manuscript input preparation from explicit archive paths."""
    global EVENT_DIR, CELL_DIR, TRACK_FILE, OUT, YEARS, B, BLOCK, SEED
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-dir", type=Path, required=True)
    parser.add_argument("--cell-dir", type=Path, required=True)
    parser.add_argument("--track-file", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    parser.add_argument("--bootstrap-repetitions", type=int, default=1000)
    parser.add_argument("--block-years", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--sections",
        nargs="+",
        choices=["global", "decomposition", "distributions", "surface"],
        default=["global", "decomposition", "distributions", "surface"],
    )
    args = parser.parse_args()
    EVENT_DIR = args.event_dir
    CELL_DIR = args.cell_dir
    TRACK_FILE = args.track_file
    OUT = args.out_dir
    YEARS = np.arange(args.start_year, args.end_year + 1)
    B = args.bootstrap_repetitions
    BLOCK = args.block_years
    SEED = args.seed
    if len(YEARS) < 3 or not 1 <= BLOCK <= len(YEARS) or B < 1:
        raise ValueError("Invalid period or bootstrap settings")
    OUT.mkdir(parents=True, exist_ok=True)
    ev = events()
    if "global" in args.sections:
        prepare_fig1(ev)
    if "decomposition" in args.sections:
        prepare_decomposition(ev)
    if "distributions" in args.sections:
        prepare_supp(ev)
    if "surface" in args.sections:
        prepare_surface_cells()
    counts = ev.assign(band_idx=-1)
    for h in ["NH", "SH"]:
        selected = counts.HEMISPHERE == h
        counts.loc[selected, "band_idx"] = band_index(counts.loc[selected, "TC_LAT"], h)
    count_path = (
        folder("TableS02_basin_decomposition") / "data" / "TableS02_event_counts.csv"
    )
    counts.groupby(["BASIN", "band_idx"]).size().rename(
        "n_events"
    ).reset_index().to_csv(count_path, index=False)
    print("Prepared manuscript figure inputs:", OUT)


if __name__ == "__main__":
    main()
