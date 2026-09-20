"""Run an explicitly selected processing stage with user-supplied data paths."""

from pathlib import Path
import argparse
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
STAGES = ["tracks", "tcpf", "thresholds", "events", "analysis"]


def run(name, *arguments):
    """Stop immediately when a processing command fails."""
    command = [sys.executable, str(ROOT / "scripts" / name), *map(str, arguments)]
    print("Running", name, flush=True)
    subprocess.run(command, check=True)


def main():
    """Run all years only when the selected stage explicitly requests them."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=STAGES + ["all"])
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mswep-dir", type=Path)
    parser.add_argument("--ibtracs-csv", type=Path)
    parser.add_argument("--override-dir", type=Path)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    args = parser.parse_args()
    stages = STAGES if args.stage == "all" else [args.stage]
    if "tracks" in stages and args.ibtracs_csv is None:
        parser.error("--ibtracs-csv is required for tracks")
    if any(x in stages for x in ["tcpf", "thresholds"]) and args.mswep_dir is None:
        parser.error("--mswep-dir is required for TCPF and thresholds")
    if args.start_year > args.end_year:
        parser.error("start-year exceeds end-year")
    archive = args.work_dir
    tracks = archive / "data/tracks_processed.parquet"
    outputs = archive / "outputs"
    period = ["--start-year", args.start_year, "--end-year", args.end_year]
    override = ["--override-dir", args.override_dir] if args.override_dir else []
    threshold = outputs / "thresholds/pot99_thresholds.parquet"
    for stage in stages:
        if stage == "tracks":
            run(
                "01_prepare_ibtracs.py",
                "--input",
                args.ibtracs_csv,
                "--output",
                tracks,
                *period,
            )
        elif stage == "tcpf":
            for year in range(args.start_year, args.end_year + 1):
                run(
                    "02_build_tcpf_database.py",
                    "--year",
                    year,
                    "--mswep-dir",
                    args.mswep_dir,
                    "--tracks",
                    tracks,
                    "--output",
                    outputs / f"tcpf/tc_gridded_{year}.parquet",
                    *override,
                )
        elif stage == "thresholds":
            run(
                "03_compute_pot99_thresholds.py",
                "--source-dir",
                args.mswep_dir,
                "--work-dir",
                archive / "threshold_work",
                "--output",
                threshold,
                *override,
                *period,
            )
        elif stage == "events":
            for year in range(args.start_year, args.end_year + 1):
                run(
                    "03_build_tcep_events.py",
                    "--tcpf",
                    outputs / f"tcpf/tc_gridded_{year}.parquet",
                    "--thresholds",
                    threshold,
                    "--output",
                    outputs / f"events/tcep_events_{year}_pot99.parquet",
                    "--cells-output",
                    outputs / f"extreme_cells/tcep_extreme_cells_{year}_pot99.parquet",
                )
        elif stage == "analysis":
            components = outputs / "components"
            inputs = outputs / "figure_inputs"
            run(
                "04_analyze_decomposition.py",
                "--events-dir",
                outputs / "events",
                "--output-dir",
                components,
                *period,
            )
            run(
                "05_prepare_figure_inputs.py",
                "--event-dir",
                outputs / "events",
                "--cell-dir",
                outputs / "extreme_cells",
                "--track-file",
                tracks,
                "--out-dir",
                inputs,
                *period,
            )
            run(
                "06_build_supplementary_tables.py",
                "--figure-root",
                inputs,
                "--component-dir",
                components,
                "--output-dir",
                outputs / "tables",
            )


if __name__ == "__main__":
    main()
