"""Render Figure S4: latitude-band histograms and kernel density estimates."""

from pathlib import Path
import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from publication_style import setup, heading, box, save, figure_paths

HERE = Path(__file__).resolve().parent
F = ["#E29A90", "#8FB4D0", "#A3C99B", "#E6C888"]
L = ["#C56C60", "#4F7EA1", "#5E8F55", "#C29B4E"]


def histogram_density(values, edges):
    counts, _ = np.histogram(values, bins=edges)
    density = counts / (len(values) * np.diff(edges))
    visible_mass = counts.sum() / len(values)
    if not np.isclose(np.sum(density * np.diff(edges)), visible_mass):
        raise RuntimeError(
            "Fig. S4 histogram density does not conserve visible probability mass"
        )
    return density


def main():
    data_dir, output_dir = figure_paths(__file__)
    setup()
    fig, axs = plt.subplots(
        1, 2, figsize=(173.57 / 25.4, 72 / 25.4), gridspec_kw={"wspace": 0.1}
    )
    pan = []
    yt = 0
    for h in ["NH", "SH"]:
        d = pd.read_csv(data_dir / f"FigS04_samples_{h}.csv")
        xmax = d.TCEP_INTENSITY_MM_3H.quantile(0.99)
        edges = np.linspace(0, xmax, 41)
        cx = np.linspace(0, xmax, 400)
        series = []
        for k, label in enumerate(d.sort_values("band_idx").band_label.unique()):
            v = d[d.band_label == label].TCEP_INTENSITY_MM_3H.to_numpy(float)
            den = histogram_density(v, edges)
            kde = gaussian_kde(v)(cx)
            yt = max(yt, den.max(), kde.max())
            series.append((label, den, kde))
        pan.append((h, xmax, edges, cx, series))
    for ax, (h, xmax, e, cx, series), p in zip(axs, pan, ["a", "b"]):
        centers = (e[:-1] + e[1:]) / 2
        w = e[1] - e[0]
        suffix = "N" if h == "NH" else "S"
        for k, (label, den, kde) in enumerate(series):
            ax.bar(
                centers,
                den,
                w,
                color=F[k],
                alpha=0.42,
                edgecolor="#F2F2F2",
                lw=0.6,
                label=f"{label.replace('-', '°-')}°{suffix}",
            )
            ax.plot(cx, kde, color=L[k], lw=2.6)
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, yt * 1.1)
        ax.set_xlabel("TCEP intensity (mm (3 h)$^{-1}$)", fontweight="bold")
        heading(ax, p, "Northern Hemisphere" if h == "NH" else "Southern Hemisphere")
        box(ax)
        ax.legend(
            title="Latitude band", loc="upper right", fontsize=7.5, title_fontsize=7.8
        )
    axs[0].set_ylabel("Probability density", fontweight="bold")
    axs[1].tick_params(axis="y", left=False, labelleft=False)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.92, bottom=0.16, wspace=0.1)
    save(fig, output_dir / "FigS04_intensity_distributions")
    plt.close(fig)
    None


if __name__ == "__main__":
    main()
