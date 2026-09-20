"""Render manuscript figures from the supplied or independently prepared CSVs."""

from pathlib import Path
import argparse
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
FIGURES = {
    "Fig01": "Fig01_global_trends",
    "Fig02": "Fig02_latitude_decomposition",
    "Fig03": "Fig03_basin_intensity_trends",
    "Fig04": "Fig04_basin_decomposition",
    "FigS01": "FigS01_TCEP_definition",
    "FigS02": "FigS02_TCEP_identification",
    "FigS03": "FigS03_poleward_shifts",
    "FigS04": "FigS04_intensity_distributions",
    "FigS05": "FigS05_land_ocean",
}


def main():
    """Dispatch one plotting process per figure and preserve the S1 source image."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", choices=FIGURES, default=list(FIGURES))
    parser.add_argument("--data-root", type=Path, default=HERE)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument(
        "--formats", nargs="+", choices=["png", "pdf", "svg", "tif"], default=["png"]
    )
    args = parser.parse_args()
    for label in args.only:
        folder = FIGURES[label]
        destination = args.output_dir / folder
        destination.mkdir(parents=True, exist_ok=True)
        if label == "FigS01":
            name = "FigS01_TCEP_definition_original.jpg"
            if (HERE / folder / name).resolve() != (destination / name).resolve():
                shutil.copy2(HERE / folder / name, destination / name)
        else:
            print("Rendering", label, flush=True)
            subprocess.run(
                [
                    sys.executable,
                    str(HERE / folder / f"plot_{label}.py"),
                    "--data-dir",
                    str(args.data_root / folder / "data"),
                    "--output-dir",
                    str(destination),
                    "--formats",
                    *args.formats,
                ],
                check=True,
            )
        caption = HERE / folder / f"{label}_caption.txt"
        if caption.resolve() != (destination / caption.name).resolve():
            shutil.copy2(caption, destination / caption.name)


if __name__ == "__main__":
    main()
