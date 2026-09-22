"""Figure 1 in the manuscript five-panel GRL layout."""

from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import setup, heading, box, save, ptext, figure_paths

HERE = Path(__file__).resolve().parent


def grid_points(d):
    lat = np.arange(-89.5, 90, 1)
    lon = np.arange(0.5, 360, 1)
    trend = np.full((180, 360), np.nan)
    pval = np.full((180, 360), np.nan)
    ii = np.rint(d.lat.to_numpy() + 89.5).astype(int)
    jj = np.rint((d.lon.to_numpy() - 0.5) % 360).astype(int)
    trend[ii, jj] = d.trend.to_numpy()
    pval[ii, jj] = d.p_value.to_numpy()
    return lat, lon, trend, pval


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup(8)
    annual = pd.read_csv(data_dir / "Fig01_annual_series.csv")
    spatial = pd.read_csv(data_dir / "Fig01_spatial_trend_1deg.csv")
    lat, lon, trend, pval = grid_points(spatial)
    # GRL extra-wide figure*: 173.57 mm (41 pc); preserve the physical canvas.
    # Preserve the manuscript geometry: wide canvas with a compact top row
    # and a dominant map/profile row below it.
    fig = plt.figure(figsize=(173.57 / 25.4, 123 / 25.4))
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[0.62, 1.0],
        left=0.065,
        right=0.965,
        top=0.95,
        bottom=0.075,
        hspace=0.32,
    )
    top = outer[0].subgridspec(1, 3, wspace=0.18)
    ylo = np.floor((annual.TCP_mean.min() - 0.5) * 2) / 2
    yhi = np.ceil((annual.TCP_mean.max() + 0.5) * 2) / 2
    for i, (dom, panel, title) in enumerate(
        [
            ("GLOBAL", "a", "Global temporal trend"),
            ("NH", "b", "Northern Hemisphere"),
            ("SH", "c", "Southern Hemisphere"),
        ]
    ):
        ax = fig.add_subplot(top[0, i])
        q = annual[annual.domain == dom].sort_values("year")
        ax.fill_between(
            q.year, q.fit_lower, q.fit_upper, color="#E79A62", alpha=0.20, lw=0
        )
        ax.plot(
            q.year,
            q.TCP_mean,
            color="#111111",
            lw=1.05,
            marker="o",
            ms=2.7,
            mec="white",
            mew=0.35,
        )
        ax.plot(q.year, q.fit, color="#111111", lw=1.8, ls=(0, (4, 2.4)))
        ax.set_xlim(1979, 2024)
        ax.set_ylim(ylo, yhi)
        ax.xaxis.set_major_locator(mticker.FixedLocator([1980, 1990, 2000, 2010, 2020]))
        ax.set_xlabel("Year", fontweight="bold")
        heading(ax, panel, title)
        box(ax)
        ax.text(
            0.045,
            0.94,
            f"{q.trend_per_decade.iloc[0]:+.2f} mm (3 h)$^{{-1}}$ decade$^{{-1}}$\n{ptext(q.p_value.iloc[0])}",
            transform=ax.transAxes,
            va="top",
        )
        if i == 0:
            ax.set_ylabel("TCEP intensity (mm (3 h)$^{-1}$)", fontweight="bold")
        else:
            ax.tick_params(axis="y", labelleft=False)
    bottom = outer[1].subgridspec(
        2,
        2,
        height_ratios=[1, 0.10],
        width_ratios=[0.79, 0.21],
        hspace=0.10,
        wspace=0.018,
    )
    pc = ccrs.PlateCarree()
    ax = fig.add_subplot(bottom[0, 0], projection=ccrs.Robinson(central_longitude=0))
    ax.set_global()
    ax.set_extent([-179.999, 179.999, -60.0, 89.999], crs=pc)
    ax.set_anchor("W")
    ax.spines["geo"].set_linewidth(1.15)
    ax.spines["geo"].set_edgecolor("#242424")
    ax.spines["geo"].set_capstyle("round")
    ax.spines["geo"].set_joinstyle("round")
    ax.add_feature(
        cfeature.LAND.with_scale("50m"),
        facecolor="#E9E8E5",
        edgecolor="none",
        zorder=0.5,
    )
    # Reorder the 0.5...359.5 grid to a -180...180 display for a 0-degree-centered map.
    display_lon = ((lon + 180) % 360) - 180
    order = np.argsort(display_lon)
    display_lon = display_lon[order]
    display_trend = trend[:, order]
    cmap = LinearSegmentedColormap.from_list(
        "trend_contrast",
        [
            (0.0, "#466F98"),
            (0.25, "#AFC7D5"),
            (0.50, "#F7F5F0"),
            (0.75, "#E7A37C"),
            (1.0, "#C65D45"),
        ],
    )
    mesh = ax.pcolormesh(
        np.arange(-180, 181),
        np.arange(-90, 91),
        display_trend,
        cmap=cmap,
        vmin=-3,
        vmax=3,
        transform=pc,
        shading="flat",
        rasterized=True,
        zorder=1,
    )
    ax.coastlines("50m", lw=0.55, color="#5A5A5A", zorder=2)
    display_llon, display_llat = np.meshgrid(display_lon, lat)
    display_pval = pval[:, order]
    sig = (
        np.isfinite(display_pval)
        & (display_pval < 0.05)
        & (display_llat >= -60)
        & (display_llat <= 80)
    )
    ax.scatter(
        display_llon[sig],
        display_llat[sig],
        s=0.14,
        c="#3A3A3A",
        alpha=0.27,
        transform=pc,
        zorder=3,
    )
    gl = ax.gridlines(
        lw=0.25,
        color="#B8B8B8",
        alpha=0.35,
        ls="--",
        draw_labels={"bottom": "x"},
        ylocs=[-60, -30, 0, 30, 60],
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.xlabel_style = {"size": 7}
    gl.ylabel_style = {"size": 7}
    # Display annotations only; names and bounds share one validated metadata source.
    regions = json.loads((HERE / "Fig01_boxes.json").read_text(encoding="utf-8"))
    for region in regions:
        number = region["box"]
        lon0 = region["lon_min"]
        lat0 = region["lat_min"]
        width = region["lon_max"] - lon0
        height = region["lat_max"] - lat0
        ax.add_patch(
            Rectangle(
                (lon0, lat0),
                width,
                height,
                transform=pc,
                fill=False,
                edgecolor="#222222",
                linewidth=0.85,
                zorder=5,
            )
        )
        ax.text(
            lon0,
            lat0 + height + 1.2,
            f"{number}",
            transform=pc,
            ha="left",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
            color="#222222",
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.82,
                "pad": 0.45,
            },
            zorder=6,
        )
    prof = fig.add_subplot(bottom[0, 1])
    valid_rows = np.isfinite(trend).sum(axis=1) > 0
    zmean = np.full(trend.shape[0], np.nan)
    zstd = np.full(trend.shape[0], np.nan)
    zmean[valid_rows] = np.nanmean(trend[valid_rows], axis=1)
    zstd[valid_rows] = np.nanstd(trend[valid_rows], axis=1)
    good = np.isfinite(zmean) & (lat >= -60) & (lat <= 80)
    prof.fill_betweenx(
        lat[good],
        (zmean - zstd)[good],
        (zmean + zstd)[good],
        color="#C96B4B",
        alpha=0.08,
        lw=0,
    )
    prof.plot(zmean[good], lat[good], color="#C96B4B", lw=1.45)
    prof.axvline(0, color="#303030", lw=0.7, ls=(0, (4, 2.4)))
    prof.set_xlim(-1.5, 1.5)
    prof.set_xticks([-1.5, 0, 1.5])
    prof.set_ylim(-60, 80)
    prof.yaxis.set_ticks_position("left")
    prof.tick_params(axis="y", left=True, labelleft=True, right=False, labelright=False)
    prof.set_yticks([-60, -30, 0, 30, 60], ["60°S", "30°S", "0°", "30°N", "60°N"])
    prof.set_xlabel("Zonal mean ± SD", fontweight="bold", fontsize=8)
    box(prof)
    fig.canvas.draw()
    map_pos = ax.get_position()
    prof_pos = prof.get_position()
    title_y = max(map_pos.y1, prof_pos.y1) + 0.010
    title_dx = 9 / 72 / fig.get_figwidth()
    for pos, panel, title in (
        (map_pos, "d", "Spatial trend map"),
        (prof_pos, "e", "Zonal profile"),
    ):
        fig.text(
            pos.x0,
            title_y,
            panel,
            ha="left",
            va="bottom",
            fontsize=10.5,
            fontweight="bold",
        )
        fig.text(
            pos.x0 + title_dx,
            title_y,
            title,
            ha="left",
            va="bottom",
            fontsize=10.2,
            fontweight="bold",
        )
    cbar_width = 0.72 * map_pos.width
    cbar_x = map_pos.x0 + (map_pos.width - cbar_width) / 2
    cbar_y = max(0.058, map_pos.y0 - 0.066)
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, 0.025])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", extend="neither")
    cb.set_ticks([-3, -1.5, 0, 1.5, 3])
    cb.ax.tick_params(axis="x", pad=1)
    cb.outline.set_linewidth(0.75)
    cb.outline.set_edgecolor("#242424")
    cb.set_label(
        "TCEP intensity trend (mm (3 h)$^{-1}$ decade$^{-1}$)",
        fontweight="bold",
        labelpad=0,
    )
    save(fig, output_dir / "Fig01_global_trends")
    plt.close(fig)


if __name__ == "__main__":
    main()
