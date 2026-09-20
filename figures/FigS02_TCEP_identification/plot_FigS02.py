"""Reproduce the MSWEP V2.8 Walaka TCPF/TCEP example from archived grid records.

No threshold, precipitation, or storm-attribution field is recalculated here.
The saved records identify storm 2018269N11220 at 2018-10-03 03:00 UTC.
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from pyproj import Geod

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from publication_style import setup, save, figure_paths


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup(8.5)
    wet = pd.read_csv(data_dir / "FigS02_Walaka_TCPF.csv")
    extreme = pd.read_csv(data_dir / "FigS02_Walaka_extreme_cells.csv")
    assert len(wet) == 13162 and len(extreme) == 1393
    assert wet.SID.nunique() == 1 and wet.TIME.nunique() == 1
    assert (extreme.PRECIP_3H_MM > extreme.POT99_MM_3H).all()
    lon0 = float(wet.TC_LON_360.iloc[0])
    lat0 = float(wet.TC_LAT.iloc[0])
    # Convert cell centres to integer 0.1-degree indices before placing pixels.
    ix = np.rint((wet.GRID_LON.to_numpy() - 0.05) * 10).astype(int)
    iy = np.rint((wet.GRID_LAT.to_numpy() - 0.05) * 10).astype(int)
    x = np.arange(ix.min(), ix.max() + 1)
    y = np.arange(iy.min(), iy.max() + 1)
    field = np.full((len(y), len(x)), np.nan)
    field[iy - y.min(), ix - x.min()] = wet.PRECIP_3H_MM
    exceed = np.zeros_like(field)
    ex = np.rint((extreme.GRID_LON.to_numpy() - 0.05) * 10).astype(int)
    ey = np.rint((extreme.GRID_LAT.to_numpy() - 0.05) * 10).astype(int)
    exceed[ey - y.min(), ex - x.min()] = 1
    bounds = [0, 1, 5, 10, 15, 20, 30, 40, 60]
    cmap = ListedColormap(
        [
            "#5E4FA2",
            "#3288BD",
            "#66C2A5",
            "#ABDDA4",
            "#E6F598",
            "#FEE08B",
            "#FDAE61",
            "#F46D43",
        ]
    )
    cmap.set_over("#9E0142")
    norm = BoundaryNorm(bounds, cmap.N)
    fig, ax = plt.subplots(figsize=(136.53 / 25.4, 139 / 25.4))
    fig.subplots_adjust(left=0.11, right=0.84, top=0.96, bottom=0.17)
    im = ax.pcolormesh(
        (np.r_[x, x[-1] + 1]) / 10,
        (np.r_[y, y[-1] + 1]) / 10,
        field,
        cmap=cmap,
        norm=norm,
        shading="flat",
        rasterized=True,
    )
    ax.contour(
        (x + 0.5) / 10,
        (y + 0.5) / 10,
        exceed,
        levels=[0.5],
        colors=["#00B8CF"],
        linewidths=0.45,
    )
    geod = Geod(ellps="WGS84")
    az = np.linspace(0, 360, 721)
    for radius, color in [(500, "#A0A0A0"), (1000, "#333333")]:
        lon, lat, _ = geod.fwd(
            np.full_like(az, lon0),
            np.full_like(az, lat0),
            az,
            np.full_like(az, radius * 1000),
        )
        ax.plot(np.mod(lon, 360), lat, color=color, lw=0.9, ls=(0, (4, 3)))
    ax.plot(lon0, lat0, "+", color="black", ms=9, mew=1.4)
    ax.set_xlim(179.2, 200.3)
    ax.set_ylim(7.0, 26.7)
    ax.set_xticks(
        [180, 185, 190, 195, 200], ["180°", "175°W", "170°W", "165°W", "160°W"]
    )
    ax.set_yticks([10, 15, 20, 25], ["10°N", "15°N", "20°N", "25°N"])
    cax = fig.add_axes([0.885, 0.17, 0.034, 0.79])
    cb = fig.colorbar(im, cax=cax, extend="max", ticks=bounds[1:])
    cb.set_label("TC precipitation (mm (3 h)$^{-1}$)", fontweight="bold")
    handles = [
        Line2D([], [], color="#00B8CF", lw=0.8, label="TCEP grid-cell outlines"),
        Line2D([], [], color="#A0A0A0", lw=0.9, ls=(0, (4, 3)), label="500-km radius"),
        Line2D([], [], color="black", marker="+", ls="none", ms=8, label="TC center"),
        Line2D(
            [], [], color="#333333", lw=0.9, ls=(0, (4, 3)), label="1,000-km radius"
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.51, 0.018),
        fontsize=7.5,
        columnspacing=1.5,
    )
    save(fig, output_dir / "FigS02_TCEP_identification")
    plt.close(fig)
    print(
        f"Walaka: {len(wet)} TCPF cells; {len(extreme)} extreme cells; mean {extreme.PRECIP_3H_MM.mean():.8f} mm/3h"
    )


if __name__ == "__main__":
    main()
