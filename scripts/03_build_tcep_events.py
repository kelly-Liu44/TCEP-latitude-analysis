"""Build storm/time-step TCEP events from TCPF cell tables and local POT99."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from tcep.events import select_extreme_cells, threshold_lookup


CELL_COLUMNS = [
    "SID",
    "TIME",
    "YEAR",
    "BASIN",
    "TC_LAT",
    "TC_LON_360",
    "GRID_LAT",
    "GRID_LON",
    "PRECIP_3H_MM",
    "PRECIP_MMHR",
    "ROW",
    "IS_INTERPOLATED_TIMESTEP",
]


def aggregate_partial_events(extreme_cells: pd.DataFrame) -> pd.DataFrame:
    """Aggregate one physical chunk into additive storm-time statistics."""
    if extreme_cells.empty:
        return pd.DataFrame()
    work = extreme_cells.copy()
    work["INTENSITY_SUM_MM_3H"] = work["PRECIP_3H_MM"].astype(np.float64)
    return work.groupby(["SID", "TIME"], as_index=False).agg(
        BASIN=("BASIN", "first"),
        TC_LAT=("TC_LAT", "first"),
        INTENSITY_SUM_MM_3H=("INTENSITY_SUM_MM_3H", "sum"),
        N_EXTREME_CELLS=("PRECIP_3H_MM", "size"),
        IS_INTERPOLATED_TIMESTEP=("IS_INTERPOLATED_TIMESTEP", "max"),
    )


def finalize_events(partials: list[pd.DataFrame]) -> pd.DataFrame:
    """Combine additive chunk statistics into one row per storm and time step."""
    if not partials:
        return pd.DataFrame(
            columns=[
                "SID",
                "TIME",
                "YEAR",
                "HEMISPHERE",
                "BASIN",
                "TC_LAT",
                "TCEP_INTENSITY_MM_3H",
                "N_EXTREME_CELLS",
                "IS_INTERPOLATED_TIMESTEP",
            ]
        )
    combined = pd.concat(partials, ignore_index=True)
    events = combined.groupby(["SID", "TIME"], as_index=False).agg(
        BASIN=("BASIN", "first"),
        TC_LAT=("TC_LAT", "first"),
        INTENSITY_SUM_MM_3H=("INTENSITY_SUM_MM_3H", "sum"),
        N_EXTREME_CELLS=("N_EXTREME_CELLS", "sum"),
        IS_INTERPOLATED_TIMESTEP=("IS_INTERPOLATED_TIMESTEP", "max"),
    )
    events["TCEP_INTENSITY_MM_3H"] = (
        events["INTENSITY_SUM_MM_3H"] / events["N_EXTREME_CELLS"]
    )
    events["YEAR"] = pd.to_datetime(events["TIME"], utc=True).dt.year.astype(np.int16)
    events["HEMISPHERE"] = np.where(events["TC_LAT"] >= 0, "NH", "SH")
    return events[
        [
            "SID",
            "TIME",
            "YEAR",
            "HEMISPHERE",
            "BASIN",
            "TC_LAT",
            "TCEP_INTENSITY_MM_3H",
            "N_EXTREME_CELLS",
            "IS_INTERPOLATED_TIMESTEP",
        ]
    ]


def main() -> None:
    """Process one year without loading the full 44-year archive."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcpf", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cells-output",
        type=Path,
        help="Optional extreme-cell table for maps and surface analyses",
    )
    args = parser.parse_args()
    thresholds = threshold_lookup(pd.read_parquet(args.thresholds))
    source = pq.ParquetFile(args.tcpf)
    missing = set(CELL_COLUMNS).difference(source.schema.names)
    if missing:
        raise ValueError(f"TCPF file is missing {sorted(missing)}")
    partial_events: list[pd.DataFrame] = []
    cell_writer = None
    temporary_cells = None
    if args.cells_output:
        args.cells_output.parent.mkdir(parents=True, exist_ok=True)
        temporary_cells = args.cells_output.with_suffix(
            args.cells_output.suffix + ".tmp"
        )
        if temporary_cells.exists():
            temporary_cells.unlink()

    # Stream physical row groups to keep peak memory independent of yearly size.
    try:
        for row_group in range(source.num_row_groups):
            cells = source.read_row_group(row_group, columns=CELL_COLUMNS).to_pandas()
            extreme_cells = select_extreme_cells(cells, thresholds)
            partial = aggregate_partial_events(extreme_cells)
            if not partial.empty:
                partial_events.append(partial)
            if temporary_cells is not None and not extreme_cells.empty:
                table = pa.Table.from_pandas(extreme_cells, preserve_index=False)
                if cell_writer is None:
                    cell_writer = pq.ParquetWriter(
                        temporary_cells, table.schema, compression="snappy"
                    )
                cell_writer.write_table(table.cast(cell_writer.schema, safe=False))
            if (row_group + 1) % 100 == 0 or row_group + 1 == source.num_row_groups:
                print(
                    f"Processed row group {row_group + 1}/{source.num_row_groups}",
                    flush=True,
                )
    finally:
        if cell_writer is not None:
            cell_writer.close()

    events = finalize_events(partial_events)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_events = args.output.with_suffix(args.output.suffix + ".tmp")
    events.to_parquet(temporary_events, index=False)
    temporary_events.replace(args.output)
    print(f"Wrote {len(events):,} events: {args.output}")
    if args.cells_output:
        if cell_writer is None or temporary_cells is None:
            raise RuntimeError("No extreme cells were produced")
        temporary_cells.replace(args.cells_output)
        print(f"Wrote streamed extreme cells: {args.cells_output}")


if __name__ == "__main__":
    main()
