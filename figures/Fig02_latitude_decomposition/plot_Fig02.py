"""Render Figure 2 with the manuscript's latitude-band decomposition layout."""

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
COMP = [
    ("tot", "Total", PALETTE["total"]),
    ("mig", "Redistribution", PALETTE["migration"]),
    ("int", "Within-band", PALETTE["intensity"]),
    ("resid", "Residual", PALETTE["residual"]),
]


def lab(v, h):
    return [
        f"{x.split('-')[0]}°-{x.split('-')[1]}° {'N' if h == 'NH' else 'S'}" for x in v
    ]


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup()
    d = pd.read_csv(data_dir / "Fig02_decomposition_by_latitude_band.csv")
    ylim = (-1.5, 1.5)
    fig, axs = plt.subplots(
        2, 1, figsize=(136.53 / 25.4, 112 / 25.4), gridspec_kw={"hspace": 0.29}
    )
    handles = []
    for ax, (h, p, t) in zip(
        axs, [("NH", "a", "Northern Hemisphere"), ("SH", "b", "Southern Hemisphere")]
    ):
        s = d[d.domain == h].sort_values("band_idx")
        x = np.arange(4) * 1.36
        w = 0.22
        for off, (k, name, color) in zip(np.array([-1.5, -0.5, 0.5, 1.5]) * w, COMP):
            bars = ax.bar(
                x + off,
                s[f"{k}_point"],
                w,
                color=color,
                edgecolor=PALETTE["edge"],
                linewidth=0.9,
                zorder=3,
                label=name,
            )
            interval_bars(
                ax, x + off, s[f"{k}_ci_low"], s[f"{k}_ci_high"], cap_width=0.07
            )
            if h == "NH":
                handles.append(bars[0])
            for xi, v, lo, hi in zip(
                x + off, s[f"{k}_point"], s[f"{k}_ci_low"], s[f"{k}_ci_high"]
            ):
                if lo > 0 or hi < 0:
                    ax.text(
                        xi,
                        hi + 0.25 if v >= 0 else lo - 0.25,
                        "*",
                        ha="center",
                        va="bottom" if v >= 0 else "top",
                        fontsize=10.8,
                        fontweight="bold",
                    )
        ax.axhline(0, color=PALETTE["zero"], lw=1.15, ls=(0, (4, 2.4)))
        ax.set_ylim(*ylim)
        ax.set_xlim(x[0] - 0.62, x[-1] + 0.62)
        ax.set_xticks(x, lab(s.band_label, h))
        ax.tick_params(top=False, right=False, pad=4)
        ax.set_ylabel(
            "Trend of TCEP intensity\n(mm (3 h)$^{-1}$ decade$^{-1}$)",
            fontweight="bold",
        )
        heading(ax, p, t)
        box(ax)
    axs[1].set_xlabel("Latitude band", fontweight="bold")
    axs[1].legend(
        handles,
        [x[1] for x in COMP],
        loc="lower right",
        ncol=2,
        fontsize=7.8,
        columnspacing=1.05,
        handlelength=1.25,
    )
    fig.subplots_adjust(left=0.17, right=0.985, top=0.942, bottom=0.105, hspace=0.29)
    save(fig, output_dir / "Fig02_latitude_decomposition")
    plt.close(fig)
    None


if __name__ == "__main__":
    main()
