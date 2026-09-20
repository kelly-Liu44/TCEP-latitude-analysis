"""Render Figure 3: six basin annual intensity series and stored trend estimates."""

from pathlib import Path
import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import setup, box, save, ptext, figure_paths

HERE = Path(__file__).resolve().parent
B = [
    ("NI", "a", "North Indian"),
    ("WP", "b", "Western Pacific"),
    ("EP", "c", "Eastern Pacific"),
    ("NA", "d", "North Atlantic"),
    ("SI", "e", "South Indian"),
    ("SP", "f", "South Pacific"),
]
POS = [
    (0.195, 0.52, 0.150, 0.22),
    (0.35, 0.52, 0.150, 0.22),
    (0.505, 0.52, 0.150, 0.22),
    (0.66, 0.52, 0.150, 0.22),
    (0.195, 0.215, 0.150, 0.22),
    (0.505, 0.215, 0.150, 0.22),
]


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup()
    d = pd.read_csv(
        data_dir / "Fig03_annual_basin_series.csv",
        keep_default_na=False,
        na_values=[""],
    )
    lo = min(d.TCP_mean.min(), d.fit_lower.min())
    hi = max(d.TCP_mean.max(), d.fit_upper.max())
    ylim = (np.floor((lo - 1) / 2) * 2, np.ceil((hi + 1) / 2) * 2)
    fig = plt.figure(figsize=(173.57 / 25.4, 92 / 25.4))
    m = fig.add_axes(
        [0.06, 0.035, 0.88, 0.93], projection=ccrs.Robinson(central_longitude=180)
    )
    m.set_global()
    m.add_feature(
        cfeature.LAND.with_scale("110m"), facecolor="#F4F2EE", edgecolor="none"
    )
    m.coastlines("110m", lw=0.45, color="#606060")
    m.gridlines(
        lw=0.28,
        color="#DFE6EE",
        alpha=0.65,
        xlocs=np.arange(-180, 181, 60),
        ylocs=np.arange(-60, 61, 30),
    )
    for pos, (bas, p, title) in zip(POS, B):
        ax = fig.add_axes(pos, facecolor="none")
        q = d[d.domain == bas].sort_values("year")
        ax.fill_between(
            q.year, q.fit_lower, q.fit_upper, color="#D7A092", alpha=0.52, lw=0
        )
        ax.plot(
            q.year,
            q.TCP_mean,
            color="#4A4A4A",
            lw=0.85,
            marker="o",
            ms=1.7,
            mec="white",
            mew=0.18,
        )
        ax.plot(q.year, q.fit, color="#2F2F2F", ls=(0, (4, 2.4)), lw=0.9)
        ax.set_xlim(1979, 2024)
        ax.set_ylim(*ylim)
        box(ax)
        ax.text(0, 1.055, p, transform=ax.transAxes, fontweight="bold", fontsize=9.6)
        ax.text(
            0.5,
            1.055,
            title,
            transform=ax.transAxes,
            ha="center",
            fontsize=7.6,
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.965,
            f"{q.trend_per_decade.iloc[0]:.2f} mm (3 h)$^{{-1}}$ decade$^{{-1}}$",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6.0,
            fontweight="bold",
        )
        ax.text(
            0.985,
            0.025,
            ptext(q.p_value.iloc[0]),
            transform=ax.transAxes,
            ha="right",
            fontsize=6.2,
            fontweight="bold",
        )
        ax.tick_params(labelsize=6.4)
        if p not in {"a", "e"}:
            ax.set_yticks([])
            ax.tick_params(axis="y", left=False, labelleft=False)
        if p not in {"e", "f"}:
            ax.set_xticks([])
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
    save(fig, output_dir / "Fig03_basin_intensity_trends")
    plt.close(fig)


if __name__ == "__main__":
    main()
