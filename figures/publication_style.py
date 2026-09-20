"""Shared manuscript font, palette, layout helpers and export settings."""

from pathlib import Path
import matplotlib as mpl
import matplotlib.font_manager as fm
from matplotlib.transforms import ScaledTranslation

FONT_DIR = Path(__file__).resolve().parent / "fonts"
PALETTE = {
    "text": "#111111",
    "edge": "#1F1F1F",
    "zero": "#707781",
    "series": "#464646",
    "trend": "#1A1A1A",
    "ci": "#AFC7D8",
    "total": "#C8CDD3",
    "migration": "#C86A56",
    "intensity": "#839FB6",
    "ocean": "#36596E",
    "ocean_light": "#97AAB4",
    "residual": "#FFFFFF",
    "land": "#D8D8D8",
    "threshold_lo": "#FBF6E6",
    "threshold_mid": "#F1C86F",
    "threshold_hi": "#C65D45",
}


def setup(font_size=8.5, top_right=True):
    for name in ("Lato-Regular.ttf", "Lato-Bold.ttf"):
        p = FONT_DIR / name
        if p.exists():
            fm.fontManager.addfont(str(p))
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Lato",
                "Arial",
                "DejaVu Sans",
                "Helvetica",
                "sans-serif",
            ],
            "font.size": font_size,
            "axes.labelsize": font_size + 0.5,
            "axes.titlesize": font_size + 1.7,
            "xtick.labelsize": font_size,
            "ytick.labelsize": font_size,
            "axes.linewidth": 1.15,
            "axes.spines.top": top_right,
            "axes.spines.right": top_right,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 3.6,
            "ytick.major.size": 3.6,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "text.color": PALETTE["text"],
            "axes.labelcolor": PALETTE["text"],
            "xtick.color": PALETTE["text"],
            "ytick.color": PALETTE["text"],
        }
    )


def heading(ax, panel, title, y=1.025):
    ax.text(
        0,
        y,
        panel,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        0,
        y,
        title,
        transform=ax.transAxes
        + ScaledTranslation(9 / 72, 0, ax.figure.dpi_scale_trans),
        ha="left",
        va="bottom",
        fontsize=10.2,
        fontweight="bold",
        clip_on=False,
    )


def box(ax):
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(1.15)
        s.set_color(PALETTE["text"])


def ptext(p):
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


EXPORT_FORMATS = ("png",)


def save(fig, stem):
    """Export one PNG by default; extra formats require an explicit request."""
    for extension in EXPORT_FORMATS:
        options = {"facecolor": "white"}
        if extension in {"png", "tif"}:
            options["dpi"] = 600
        if extension == "tif":
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(stem.with_suffix("." + extension), **options)


def interval_bars(ax, x, lower, upper, cap_width=0.07):
    """Draw the actual percentile interval independently of the point estimate.

    A percentile bootstrap interval need not contain the original estimate.
    Clipping negative yerr values would silently change that interval.
    """
    import numpy as np

    x = np.asarray(x, float)
    lower = np.asarray(lower, float)
    upper = np.asarray(upper, float)
    good = np.isfinite(x) & np.isfinite(lower) & np.isfinite(upper)
    ax.vlines(
        x[good], lower[good], upper[good], color=PALETTE["edge"], lw=1.05, zorder=5
    )
    ax.hlines(
        lower[good],
        x[good] - cap_width,
        x[good] + cap_width,
        color=PALETTE["edge"],
        lw=1.05,
        zorder=5,
    )
    ax.hlines(
        upper[good],
        x[good] - cap_width,
        x[good] + cap_width,
        color=PALETTE["edge"],
        lw=1.05,
        zorder=5,
    )


def figure_paths(script):
    """Use each figure folder as the single default image location."""
    import argparse

    global EXPORT_FORMATS
    folder = Path(script).resolve().parent
    parser = argparse.ArgumentParser(description=folder.name.replace("_", " "))
    parser.add_argument("--data-dir", type=Path, default=folder / "data")
    parser.add_argument("--output-dir", type=Path, default=folder)
    parser.add_argument(
        "--formats", nargs="+", choices=["png", "pdf", "svg", "tif"], default=["png"]
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    EXPORT_FORMATS = tuple(args.formats)
    return args.data_dir, args.output_dir
