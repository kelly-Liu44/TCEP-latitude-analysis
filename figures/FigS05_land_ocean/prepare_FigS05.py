"""Rebuild Figure S5 summaries from the cached annual extreme-cell band table."""

from pathlib import Path
import sys, numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tcep import figure_inputs as prep
import argparse

HERE = Path(__file__).resolve().parent
YEARS = prep.YEARS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=HERE / "data/FigS05_annual_surface_band.csv"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=HERE / "data",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ab = pd.read_csv(args.input)
    annual = []
    decomp = []
    for h in ["NH", "SH"]:
        for surface in ["land", "ocean"]:
            x = ab[(ab.domain == h) & (ab.surface == surface)].rename(
                columns={"year": "YEAR"}
            )
            n, v = prep.matrices(x)
            st = prep.stats_from(n, v)
            fit, res = prep.sen_line(st["mean"])
            rng = np.random.default_rng(prep.SEED)
            boot = np.array(
                [
                    prep.sen_line(fit + res[prep.block_idx(rng)])[0]
                    for _ in range(prep.B)
                ]
            )
            annual.append(
                pd.DataFrame(
                    {
                        "domain": h,
                        "surface": surface,
                        "year": YEARS,
                        "N_t": n.sum(1),
                        "TCP_mean": st["mean"],
                        "fit": fit,
                        "fit_lower": np.quantile(boot, 0.025, 0),
                        "fit_upper": np.quantile(boot, 0.975, 0),
                        "trend_per_decade": st["slope"],
                        "p_value": st["p"],
                    }
                )
            )
            for k, label in enumerate(prep.labels(h)):
                decomp.append(
                    {
                        "domain": h,
                        "surface": surface,
                        "band_idx": k,
                        "band_label": label,
                        "migration": st["mig"][k],
                        "intensity": st["inten"][k],
                        "residual": st["resid"][k],
                        "total": st["total"][k],
                    }
                )
    pd.concat(annual).to_csv(
        args.output_dir / "FigS05_annual_land_ocean.csv", index=False
    )
    pd.DataFrame(decomp).to_csv(
        args.output_dir / "FigS05_land_ocean_decomposition.csv", index=False
    )


if __name__ == "__main__":
    main()
