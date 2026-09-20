"""Extract the existing Walaka example from the project's archived V2.8 records."""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
SID = "2018269N11220"
TIME = pd.Timestamp("2018-10-03 03:00:00", tz="UTC")


def main():
    global PROJECT, HERE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-dir",
        type=Path,
        required=True,
        help="Processed archive with outputs/tcpf, outputs/extreme_cells and outputs/events",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="S2 folder; CSVs are written to its data/ subfolder",
    )
    args = parser.parse_args()
    PROJECT = args.project_dir
    HERE = args.output_dir
    (HERE / "data").mkdir(parents=True, exist_ok=True)
    filters = [("SID", "==", SID), ("TIME", "==", TIME)]
    wet = pd.read_parquet(
        PROJECT / "outputs/tcpf/tc_gridded_2018.parquet", filters=filters
    )
    extreme = pd.read_parquet(
        PROJECT / "outputs/extreme_cells/tcep_extreme_cells_2018_pot99.parquet",
        filters=filters,
    )
    event = pd.read_parquet(
        PROJECT / "outputs/events/tcep_events_2018_pot99.parquet", filters=filters
    )
    assert len(event) == 1 and len(wet) == 13162 and len(extreme) == 1393
    assert event.N_EXTREME_CELLS.iloc[0] == len(extreme)
    assert np.isclose(
        extreme.PRECIP_3H_MM.astype("float64").mean(),
        event.TCEP_INTENSITY_MM_3H.iloc[0],
        atol=1e-10,
    )
    wet.to_csv(HERE / "data/FigS02_Walaka_TCPF.csv", index=False)
    extreme.to_csv(HERE / "data/FigS02_Walaka_extreme_cells.csv", index=False)
    print(
        "Extracted the stored Walaka records; no TCPF masks or POT99 thresholds recalculated."
    )


if __name__ == "__main__":
    main()
