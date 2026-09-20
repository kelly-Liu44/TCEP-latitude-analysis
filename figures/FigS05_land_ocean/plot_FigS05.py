"""Render Figure S5 with extreme-cell land/ocean labels in a 2 × 2 layout."""

from pathlib import Path
import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import setup, PALETTE, heading, box, save, ptext, figure_paths

HERE = Path(__file__).resolve().parent


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup()
    a = pd.read_csv(data_dir / "FigS05_annual_land_ocean.csv")
    d = pd.read_csv(data_dir / "FigS05_land_ocean_decomposition.csv")
    fig, axs = plt.subplots(
        2,
        2,
        figsize=(173.57 / 25.4, 118 / 25.4),
        gridspec_kw={"wspace": 0.1, "hspace": 0.31},
    )
    cols = {"land": PALETTE["migration"], "ocean": PALETTE["ocean"]}
    ts = (
        np.floor((a.TCP_mean.min() - 3) / 5) * 5,
        np.ceil((a.TCP_mean.max() + 3) / 5) * 5,
    )
    vals = d[["migration", "intensity"]].to_numpy(float)
    bottom_lo = float(vals.min())
    bottom_hi = float(vals.max())
    bottom_pad = 0.08 * (bottom_hi - bottom_lo)
    for j, h in enumerate(["NH", "SH"]):
        ax = axs[0, j]
        for surf in ["ocean", "land"]:
            q = a[(a.domain == h) & (a.surface == surf)].sort_values("year")
            ax.fill_between(
                q.year,
                q.fit_lower,
                q.fit_upper,
                color=cols[surf],
                alpha=0.18,
                lw=0,
                zorder=1,
            )
            ax.plot(
                q.year,
                q.TCP_mean,
                color=cols[surf],
                lw=1.55,
                marker="o",
                ms=3.3,
                mec="white",
                mew=0.35,
                label=surf.title(),
                zorder=3,
            )
            ax.plot(
                q.year, q.fit, color=cols[surf], lw=1.55, ls=(0, (4, 2.4)), zorder=4
            )
            ax.text(
                0.975,
                0.975 if surf == "ocean" else 0.895,
                f"{surf.title()}: {q.trend_per_decade.iloc[0]:.2f} mm (3 h)$^{{-1}}$ decade$^{{-1}}$ ({ptext(q.p_value.iloc[0])})",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7.6,
                fontweight="bold",
                color=cols[surf],
                zorder=6,
                bbox={
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.82,
                    "pad": 0.4,
                },
            )
        ax.set_ylim(*ts)
        ax.set_xlim(1979, 2024)
        ax.set_xlabel("Year", fontweight="bold")
        heading(ax, "ab"[j], f"{'Northern' if h == 'NH' else 'Southern'} Hemisphere")
        box(ax)
        if j == 0:
            ax.set_ylabel("TCEP mean intensity (mm (3 h)$^{-1}$)", fontweight="bold")
        else:
            ax.tick_params(axis="y", left=False, labelleft=False)
        ax = axs[1, j]
        q = d[d.domain == h]
        x = np.arange(4) * 1.42
        w = 0.23
        spec = [
            ("land", "migration", "Land redistribution", PALETTE["migration"], "//"),
            ("land", "intensity", "Land within-band", "#E97252", ""),
            (
                "ocean",
                "migration",
                "Ocean redistribution",
                PALETTE["ocean_light"],
                "//",
            ),
            ("ocean", "intensity", "Ocean within-band", PALETTE["ocean"], ""),
        ]
        for off, (surf, k, name, c, ha) in zip(
            np.array([-1.5, -0.5, 0.5, 1.5]) * w, spec
        ):
            ax.bar(
                x + off,
                q[q.surface == surf].sort_values("band_idx")[k],
                w,
                color=c,
                edgecolor=PALETTE["edge"],
                lw=1,
                hatch=ha,
                label=name,
                alpha=0.82 if ha else 1,
            )
        ax.set_ylim(bottom_lo - bottom_pad, bottom_hi + bottom_pad)
        ax.axhline(0, color=PALETTE["zero"], lw=1.15, ls=(0, (4, 2.4)))
        ax.set_xticks(
            x,
            [
                f"{z.split('-')[0]}°-{z.split('-')[1]}°{'N' if h == 'NH' else 'S'}"
                for z in q[q.surface == "land"].sort_values("band_idx").band_label
            ],
        )
        ax.set_xlabel("Latitude band", fontweight="bold")
        heading(ax, "cd"[j], f"{'Northern' if h == 'NH' else 'Southern'} Hemisphere")
        box(ax)
        if j == 0:
            ax.set_ylabel(
                "Contribution to trend\n(mm (3 h)$^{-1}$ decade$^{-1}$)",
                fontweight="bold",
                labelpad=1,
            )
        else:
            ax.tick_params(axis="y", left=False, labelleft=False)
    axs[0, 1].legend(
        loc="lower right",
        ncol=2,
        fontsize=6.2,
        frameon=False,
        bbox_to_anchor=(0.97, 0.025),
        borderaxespad=0,
        columnspacing=0.6,
        handlelength=1.0,
        handletextpad=0.35,
    )
    axs[1, 1].legend(
        loc="lower right",
        ncol=2,
        fontsize=6.0,
        frameon=False,
        bbox_to_anchor=(0.97, 0.025),
        borderaxespad=0,
        columnspacing=0.45,
        handlelength=0.8,
        handletextpad=0.3,
        labelspacing=0.2,
    )
    fig.subplots_adjust(
        left=0.105, right=0.992, top=0.95, bottom=0.09, wspace=0.1, hspace=0.31
    )
    save(fig, output_dir / "FigS05_land_ocean")
    plt.close(fig)
    None


if __name__ == "__main__":
    main()
