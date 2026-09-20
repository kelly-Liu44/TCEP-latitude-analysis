"""Render Figure 4: six basin decompositions in the manuscript's 2 × 3 layout."""

from pathlib import Path
import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import (
    setup,
    PALETTE,
    heading,
    box,
    save,
    interval_bars,
    figure_paths,
)

HERE = Path(__file__).resolve().parent
B = [
    ("WP", "a", "Western Pacific", "N"),
    ("EP", "b", "Eastern Pacific", "N"),
    ("NA", "c", "North Atlantic", "N"),
    ("NI", "d", "North Indian", "N"),
    ("SI", "e", "South Indian", "S"),
    ("SP", "f", "South Pacific", "S"),
]
C = [
    ("tot", "Total", PALETTE["total"]),
    ("mig", "Redistribution", PALETTE["migration"]),
    ("int", "Within-band", PALETTE["intensity"]),
    ("resid", "Residual", PALETTE["residual"]),
]


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup()
    d = pd.read_csv(
        data_dir / "Fig04_decomposition_by_basin_and_band.csv",
        keep_default_na=False,
        na_values=[""],
    )
    ylim = (-1.5, 1.5)
    fig, axs = plt.subplots(
        2,
        3,
        figsize=(173.57 / 25.4, 112 / 25.4),
        gridspec_kw={"wspace": 0.08, "hspace": 0.30},
    )
    handles = []
    for ix, (ax, (bas, p, title, suf)) in enumerate(zip(axs.ravel(), B)):
        q = d[d.domain == bas].sort_values("band_idx")
        x = np.arange(4) * 1.34
        w = 0.23
        for off, (k, name, color) in zip(np.array([-1.5, -0.5, 0.5, 1.5]) * w, C):
            y = q[f"{k}_point"].to_numpy(float).copy()
            valid = ~((bas == "NI") & (q.band_idx.to_numpy() == 3))
            y[~valid] = np.nan
            bars = ax.bar(
                x + off, y, w, color=color, edgecolor=PALETTE["edge"], lw=1, zorder=3
            )
            interval_bars(
                ax,
                (x + off)[valid],
                q.loc[valid, f"{k}_ci_low"],
                q.loc[valid, f"{k}_ci_high"],
                cap_width=0.065,
            )
            if ix == 0:
                handles.append(bars[0])
            for xx, v, l, h in zip(x + off, y, q[f"{k}_ci_low"], q[f"{k}_ci_high"]):
                if l > 0 or h < 0:
                    ax.text(
                        xx,
                        h + 0.2 if v >= 0 else l - 0.2,
                        "*",
                        ha="center",
                        va="bottom" if v >= 0 else "top",
                        fontsize=9,
                        fontweight="bold",
                    )
        ax.set_xlim(x[0] - 0.62, x[-1] + 0.62)
        if bas == "NI":
            ax.text(x[3], 0.08, "No events", ha="center", fontsize=6.2)
        ax.axhline(0, color=PALETTE["zero"], lw=1.15, ls=(0, (4, 2.4)))
        ax.set_ylim(*ylim)
        ax.set_xticks(
            x,
            [f"{z.split('-')[0]}°-{z.split('-')[1]}°{suf}" for z in q.band_label],
            fontsize=7.6,
        )
        ax.set_xlabel("Latitude band", fontweight="bold")
        heading(ax, p, title)
        box(ax)
        ax.tick_params(top=False, right=False)
        if ix % 3:
            ax.tick_params(axis="y", left=False, labelleft=False)
        if ix < 3:
            ax.set_xlabel("")
    axs.ravel()[5].legend(
        handles,
        [z[1] for z in C],
        loc="lower right",
        ncol=2,
        fontsize=6.1,
        frameon=False,
        bbox_to_anchor=(0.965, 0.025),
        borderaxespad=0,
        columnspacing=0.55,
        handlelength=0.9,
        handletextpad=0.35,
        labelspacing=0.25,
    )
    fig.supylabel(
        "TCEP intensity trend (mm (3 h)$^{-1}$ decade$^{-1}$)",
        x=0.018,
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(
        left=0.105, right=0.975, top=0.94, bottom=0.12, wspace=0.08, hspace=0.30
    )
    save(fig, output_dir / "Fig04_basin_decomposition")
    plt.close(fig)
    None


if __name__ == "__main__":
    main()
