"""Render Figure S3: annual TC-track and TCEP-event absolute latitudes."""

from pathlib import Path
import sys, pandas as pd, matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import setup, PALETTE, heading, box, save, ptext, figure_paths

HERE = Path(__file__).resolve().parent


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup(9.5)
    d = pd.read_csv(data_dir / "FigS03_annual_mean_latitude.csv")
    fig, axs = plt.subplots(
        2,
        2,
        figsize=(173.57 / 25.4, 118 / 25.4),
        gridspec_kw={"wspace": 0.12, "hspace": 0.34},
    )
    P = [
        ("track", "NH", "a", "Northern Hemisphere TC-track"),
        ("track", "SH", "b", "Southern Hemisphere TC-track"),
        ("event", "NH", "c", "Northern Hemisphere TCEP"),
        ("event", "SH", "d", "Southern Hemisphere TCEP"),
    ]
    for ax, (kind, h, p, title) in zip(axs.ravel(), P):
        q = d[(d.record_type == kind) & (d.domain == h)].sort_values("year")
        c = "#8FB4D0" if kind == "track" else "#E29A90"
        line = "#4F7EA1" if kind == "track" else "#C56C60"
        ax.fill_between(q.year, q.fit_lower, q.fit_upper, color=c, alpha=0.85, lw=0)
        ax.plot(
            q.year,
            q.mean_lat,
            color=PALETTE["series"],
            lw=1.4,
            marker="o",
            ms=3.2,
            mec="white",
            mew=0.35,
        )
        ax.plot(q.year, q.fit, color=line, lw=1.65, ls=(0, (4, 2.4)))
        ax.text(
            0.97,
            0.06,
            f"{q.trend_per_decade.iloc[0]:.2f}° decade$^{{-1}}$ ({ptext(q.p_value.iloc[0])})",
            transform=ax.transAxes,
            ha="right",
            fontweight="bold",
            fontsize=8.5,
        )
        heading(ax, p, title, 1.025)
        box(ax)
    for ax in axs[:, 0]:
        ax.set_ylabel("Mean absolute latitude (°)", fontweight="bold")
    for ax in axs[:, 1]:
        ax.tick_params(axis="y", left=False, labelleft=False)
    for ax in axs[1]:
        ax.set_xlabel("Year", fontweight="bold")
    fig.subplots_adjust(
        left=0.085, right=0.992, top=0.95, bottom=0.09, wspace=0.12, hspace=0.34
    )
    save(fig, output_dir / "FigS03_poleward_shifts")
    plt.close(fig)
    None


if __name__ == "__main__":
    main()
