"""Compute local wet-period POT99 thresholds from native MSWEP V2.8 fields.

The calculation uses a resumable 500-bin global histogram. All values strictly
greater than 0.3 mm/3h enter the wet sample. Values above the histogram ceiling
remain in the sample through a separate overflow count; they are never silently
discarded. Cells whose requested rank lies in the overflow tail are reread and resolved
by exact, unbounded nearest-rank selection. In-range thresholds use bin centres.
"""

from __future__ import annotations

import argparse
import calendar
import os
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset


NLAT, NLON = 1800, 3600
TOTAL_CELLS = NLAT * NLON
WET_MIN_MM_3H = 0.3
HISTOGRAM_MAX_MM_3H = 120.0
N_BINS = 500
BIN_EDGES = np.linspace(
    WET_MIN_MM_3H, HISTOGRAM_MAX_MM_3H, N_BINS + 1, dtype=np.float32
)
BIN_CENTRES = ((BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2.0).astype(np.float32)
LATITUDES_PER_BLOCK = 10
ROWS_PER_BLOCK = LATITUDES_PER_BLOCK * NLON
N_BLOCKS = NLAT // LATITUDES_PER_BLOCK
HEADER_DTYPE = np.dtype(np.int32)


def expected_names(year: int) -> set[str]:
    """Return the exact three-hour filename set for one complete year."""
    days = 366 if calendar.isleap(year) else 365
    return {
        f"{year}{day:03d}.{hour:02d}.nc"
        for day in range(1, days + 1)
        for hour in range(0, 24, 3)
    }


def file_map(source_dir: Path, override_dir: Path | None, year: int) -> dict[str, Path]:
    """Build a complete index with documented override files taking precedence."""
    paths = {path.name: path for path in source_dir.glob(f"{year}*.nc")}
    if override_dir:
        paths.update({path.name: path for path in override_dir.glob(f"{year}*.nc")})
    missing = sorted(expected_names(year) - paths.keys())
    if missing:
        raise FileNotFoundError(
            f"{year}: missing {len(missing)} fields, e.g. {missing[:5]}"
        )
    return {name: paths[name] for name in sorted(expected_names(year))}


def longitude_order(sample: Path) -> np.ndarray:
    """Return the column permutation from raw -180..180 to internal 0..360."""
    with Dataset(sample, "r") as dataset:
        longitude = np.asarray(dataset.variables["lon"][:], dtype=np.float64)
        units = (
            str(dataset.variables["precipitation"].getncattr("units"))
            .lower()
            .replace(" ", "")
        )
    if units not in {"mm/3h", "mm3h-1", "mm/3hr", "mm/3hour"}:
        raise ValueError(f"Unexpected precipitation units in {sample}: {units!r}")
    order = np.argsort(longitude % 360.0)
    converted = (longitude[order] % 360.0).astype(np.float64)
    if not np.allclose(converted, np.arange(0.05, 360.0, 0.1), atol=2e-4):
        raise ValueError(f"Longitude conversion failed for {sample}")
    return order


def annual_histogram_path(work_dir: Path, year: int) -> Path:
    """Return the per-year histogram path used for safe resume."""
    return work_dir / "annual" / f"hist_{year}.uint16.bin"


def annual_overflow_path(work_dir: Path, year: int) -> Path:
    """Return the per-year overflow-count path used for safe resume."""
    return work_dir / "annual" / f"overflow_{year}.uint16.bin"


def annual_done_path(work_dir: Path, year: int) -> Path:
    """Return the marker created only after a complete annual source scan."""
    return work_dir / "annual" / f"scan_{year}.done"


def cumulative_block_path(work_dir: Path, block_index: int) -> Path:
    """Return one independently resumable latitude-block histogram path."""
    return work_dir / "cumulative_blocks" / f"hist_{block_index:03d}.bin"


def cumulative_overflow_block_path(work_dir: Path, block_index: int) -> Path:
    """Return the matching latitude-block overflow-count path."""
    return work_dir / "cumulative_blocks" / f"overflow_{block_index:03d}.bin"


def exact_unbounded_thresholds(
    source_dir: Path,
    override_dir: Path | None,
    rows: np.ndarray,
    expected_wet_counts: np.ndarray,
    start_year: int,
    end_year: int,
) -> pd.DataFrame:
    """Compute exact nearest-rank POT99 for cells whose rank exceeds 120 mm/3h.

    Only the unresolved cells are reread. Their values are never clipped or
    binned, so this refinement has no precipitation ceiling.
    """
    rows = np.asarray(rows, dtype=np.int64)
    expected_wet_counts = np.asarray(expected_wet_counts, dtype=np.uint64)
    if rows.size == 0:
        return pd.DataFrame(columns=["ROW", "POT99_MM_3H"])

    # Convert canonical flattened rows to native MSWEP latitude and longitude indices.
    latitude_indices = rows // NLON
    internal_longitude_indices = rows % NLON
    sample_files = file_map(source_dir, override_dir, start_year)
    order = longitude_order(next(iter(sample_files.values())))
    raw_longitude_indices = order[internal_longitude_indices]

    # The unresolved cells form a compact region; one rectangular read per file
    # is substantially faster than 92 independent NetCDF point reads.
    latitude_start = int(latitude_indices.min())
    latitude_stop = int(latitude_indices.max()) + 1
    longitude_start = int(raw_longitude_indices.min())
    longitude_stop = int(raw_longitude_indices.max()) + 1
    local_latitude_indices = latitude_indices - latitude_start
    local_longitude_indices = raw_longitude_indices - longitude_start
    collected: list[list[float]] = [[] for _ in rows]

    for year in range(start_year, end_year + 1):
        files = file_map(source_dir, override_dir, year)
        for index, path in enumerate(files.values(), start=1):
            with Dataset(path, "r") as dataset:
                variable = dataset.variables["precipitation"]
                units = str(variable.getncattr("units")).lower().replace(" ", "")
                if units not in {"mm/3h", "mm3h-1", "mm/3hr", "mm/3hour"}:
                    raise ValueError(
                        f"Unexpected precipitation units in {path}: {units!r}"
                    )
                rectangle = np.asarray(
                    np.ma.filled(
                        variable[
                            0,
                            latitude_start:latitude_stop,
                            longitude_start:longitude_stop,
                        ],
                        np.nan,
                    ),
                    dtype=np.float32,
                )
            values = rectangle[local_latitude_indices, local_longitude_indices]
            for cell_index, value in enumerate(values):
                if np.isfinite(value) and value > WET_MIN_MM_3H:
                    collected[cell_index].append(float(value))
        print(
            f"[exact refinement] {year}: {len(files)}/{len(files)} fields", flush=True
        )

    exact_values = np.empty(len(rows), dtype=np.float32)
    for cell_index, values in enumerate(collected):
        observed_count = len(values)
        expected_count = int(expected_wet_counts[cell_index])
        if observed_count != expected_count:
            raise RuntimeError(
                f"ROW {rows[cell_index]} exact wet count {observed_count} "
                f"does not match histogram count {expected_count}"
            )
        ordered = np.sort(np.asarray(values, dtype=np.float32))
        rank = int(np.ceil(0.99 * observed_count))
        exact_values[cell_index] = ordered[rank - 1]
    return pd.DataFrame(
        {
            "ROW": rows,
            "POT99_MM_3H": exact_values,
            "THRESHOLD_METHOD": "exact nearest-rank (unbounded tail refinement)",
            "IS_EXACT_TAIL_REFINEMENT": True,
        }
    )


def block_rows(block_index: int) -> tuple[int, int]:
    """Return global flattened-row bounds for one latitude block."""
    start = block_index * ROWS_PER_BLOCK
    return start, min(start + ROWS_PER_BLOCK, TOTAL_CELLS)


def read_block_year(path: Path) -> int | None:
    """Read the last merged year stored in an atomic block header."""
    if not path.exists():
        return None
    with path.open("rb") as stream:
        header = stream.read(HEADER_DTYPE.itemsize)
    if len(header) != HEADER_DTYPE.itemsize:
        raise ValueError(f"Invalid cumulative block header: {path}")
    return int(np.frombuffer(header, dtype=HEADER_DTYPE)[0])


def year_is_fully_merged(work_dir: Path, year: int) -> bool:
    """Return true only when every cumulative block already includes this year."""
    years = [
        read_block_year(cumulative_block_path(work_dir, index))
        for index in range(N_BLOCKS)
    ]
    return all(value is not None and value >= year for value in years)


def scan_year(
    source_dir: Path, override_dir: Path | None, work_dir: Path, year: int
) -> None:
    """Scan one year into uint16 bin and overflow counts for every grid cell."""
    histogram_path = annual_histogram_path(work_dir, year)
    overflow_path = annual_overflow_path(work_dir, year)
    expected_histogram_bytes = TOTAL_CELLS * N_BINS * np.dtype(np.uint16).itemsize
    expected_overflow_bytes = TOTAL_CELLS * np.dtype(np.uint16).itemsize
    if (
        histogram_path.exists()
        and histogram_path.stat().st_size == expected_histogram_bytes
        and overflow_path.exists()
        and overflow_path.stat().st_size == expected_overflow_bytes
        and annual_done_path(work_dir, year).exists()
    ):
        print(f"[{year}] complete annual scan retained", flush=True)
        return

    # Remove only incomplete intermediates for this exact year.
    for path in (histogram_path, overflow_path, annual_done_path(work_dir, year)):
        if path.exists():
            path.unlink()
    histogram_path.parent.mkdir(parents=True, exist_ok=True)
    files = file_map(source_dir, override_dir, year)
    order = longitude_order(next(iter(files.values())))
    counts = np.memmap(
        histogram_path, dtype=np.uint16, mode="w+", shape=(TOTAL_CELLS, N_BINS)
    )
    overflow = np.memmap(overflow_path, dtype=np.uint16, mode="w+", shape=TOTAL_CELLS)
    counts[:] = 0
    overflow[:] = 0

    for index, path in enumerate(files.values(), start=1):
        with Dataset(path, "r") as dataset:
            variable = dataset.variables["precipitation"]
            units = str(variable.getncattr("units")).lower().replace(" ", "")
            if units not in {"mm/3h", "mm3h-1", "mm/3hr", "mm/3hour"}:
                raise ValueError(f"Unexpected precipitation units in {path}: {units!r}")
            raw = np.asarray(np.ma.filled(variable[0], np.nan), dtype=np.float32)[
                :, order
            ]
        finite = np.isfinite(raw)
        in_range = finite & (raw > WET_MIN_MM_3H) & (raw <= HISTOGRAM_MAX_MM_3H)
        flat_rows = np.flatnonzero(in_range)
        if flat_rows.size:
            values = raw.ravel()[flat_rows]
            bin_index = np.searchsorted(BIN_EDGES[1:], values, side="left").clip(
                0, N_BINS - 1
            )
            np.add.at(counts, (flat_rows, bin_index), 1)
        high_rows = np.flatnonzero(finite & (raw > HISTOGRAM_MAX_MM_3H))
        if high_rows.size:
            np.add.at(overflow, high_rows, 1)
        if index % 400 == 0 or index == len(files):
            counts.flush()
            overflow.flush()
            print(f"[{year}] {index}/{len(files)} fields", flush=True)
    counts.flush()
    overflow.flush()
    del counts, overflow
    # A marker distinguishes a fully scanned mmap from a merely preallocated one.
    annual_done_path(work_dir, year).write_text("complete\n", encoding="ascii")


def merge_year(work_dir: Path, year: int) -> None:
    """Atomically merge one annual scan into independently resumable blocks."""
    if not annual_done_path(work_dir, year).exists():
        raise RuntimeError(f"{year}: annual scan has no completion marker")
    annual = np.memmap(
        annual_histogram_path(work_dir, year),
        dtype=np.uint16,
        mode="r",
        shape=(TOTAL_CELLS, N_BINS),
    )
    annual_overflow = np.memmap(
        annual_overflow_path(work_dir, year),
        dtype=np.uint16,
        mode="r",
        shape=TOTAL_CELLS,
    )
    block_dir = cumulative_block_path(work_dir, 0).parent
    block_dir.mkdir(parents=True, exist_ok=True)

    for block_index in range(N_BLOCKS):
        histogram_path = cumulative_block_path(work_dir, block_index)
        overflow_path = cumulative_overflow_block_path(work_dir, block_index)
        histogram_year = read_block_year(histogram_path)
        overflow_year = read_block_year(overflow_path)
        if histogram_year != overflow_year:
            raise RuntimeError(
                f"Block {block_index}: histogram and overflow years disagree"
            )
        if histogram_year is not None and histogram_year >= year:
            continue
        if histogram_year is not None and histogram_year != year - 1:
            raise RuntimeError(
                f"Block {block_index}: cannot merge {year} after {histogram_year}"
            )

        start, end = block_rows(block_index)
        rows = end - start
        temporary_histogram = histogram_path.with_suffix(".tmp")
        temporary_overflow = overflow_path.with_suffix(".tmp")
        for temporary in (temporary_histogram, temporary_overflow):
            if temporary.exists():
                temporary.unlink()

        with temporary_histogram.open("wb") as stream:
            stream.write(np.asarray([year], dtype=HEADER_DTYPE).tobytes())
            stream.truncate(
                HEADER_DTYPE.itemsize + rows * N_BINS * np.dtype(np.uint32).itemsize
            )
        with temporary_overflow.open("wb") as stream:
            stream.write(np.asarray([year], dtype=HEADER_DTYPE).tobytes())
            stream.truncate(HEADER_DTYPE.itemsize + rows * np.dtype(np.uint32).itemsize)

        updated = np.memmap(
            temporary_histogram,
            dtype=np.uint32,
            mode="r+",
            offset=HEADER_DTYPE.itemsize,
            shape=(rows, N_BINS),
        )
        updated_high = np.memmap(
            temporary_overflow,
            dtype=np.uint32,
            mode="r+",
            offset=HEADER_DTYPE.itemsize,
            shape=rows,
        )
        previous = None
        previous_high = None
        if histogram_path.exists():
            previous = np.memmap(
                histogram_path,
                dtype=np.uint32,
                mode="r",
                offset=HEADER_DTYPE.itemsize,
                shape=(rows, N_BINS),
            )
            previous_high = np.memmap(
                overflow_path,
                dtype=np.uint32,
                mode="r",
                offset=HEADER_DTYPE.itemsize,
                shape=rows,
            )
        if previous is None:
            updated[:] = annual[start:end]
            updated_high[:] = annual_overflow[start:end]
        else:
            updated[:] = previous[:] + annual[start:end]
            updated_high[:] = previous_high[:] + annual_overflow[start:end]
        updated.flush()
        updated_high.flush()
        del updated, updated_high, previous, previous_high
        os.replace(temporary_histogram, histogram_path)
        os.replace(temporary_overflow, overflow_path)
        if block_index % 20 == 0 or block_index == N_BLOCKS - 1:
            print(f"[{year}] merged block {block_index + 1}/{N_BLOCKS}", flush=True)

    del annual, annual_overflow
    if any(
        read_block_year(cumulative_block_path(work_dir, index)) != year
        for index in range(N_BLOCKS)
    ):
        raise RuntimeError(f"{year}: not every cumulative block reached this year")
    annual_histogram_path(work_dir, year).unlink()
    annual_overflow_path(work_dir, year).unlink()
    annual_done_path(work_dir, year).unlink()
    print(f"[{year}] merged and annual intermediates removed", flush=True)


def derive_thresholds(
    source_dir: Path,
    override_dir: Path | None,
    work_dir: Path,
    output: Path,
    start_year: int,
    end_year: int,
) -> None:
    """Derive POT99 and exactly refine every rank above the histogram ceiling."""
    if any(
        read_block_year(cumulative_block_path(work_dir, index)) != end_year
        for index in range(N_BLOCKS)
    ):
        raise RuntimeError(
            "Not every cumulative block includes the requested final year"
        )
    frames: list[pd.DataFrame] = []
    unresolved_rows: list[np.ndarray] = []
    unresolved_counts: list[np.ndarray] = []
    for block_index in range(N_BLOCKS):
        start, end = block_rows(block_index)
        rows_in_block = end - start
        counts = np.memmap(
            cumulative_block_path(work_dir, block_index),
            dtype=np.uint32,
            mode="r",
            offset=HEADER_DTYPE.itemsize,
            shape=(rows_in_block, N_BINS),
        )
        overflow = np.memmap(
            cumulative_overflow_block_path(work_dir, block_index),
            dtype=np.uint32,
            mode="r",
            offset=HEADER_DTYPE.itemsize,
            shape=rows_in_block,
        )
        block = np.asarray(counts, dtype=np.uint32)
        in_range_total = block.sum(axis=1, dtype=np.uint64)
        total = in_range_total + np.asarray(overflow, dtype=np.uint64)
        valid = total > 0
        target = np.ceil(0.99 * total).astype(np.uint64)
        unresolved = valid & (target > in_range_total)
        if unresolved.any():
            unresolved_rows.append(start + np.flatnonzero(unresolved))
            unresolved_counts.append(total[unresolved])
        resolved = valid & ~unresolved
        if resolved.any():
            cumulative = np.cumsum(block[resolved], axis=1, dtype=np.uint64)
            indices = np.argmax(cumulative >= target[resolved, None], axis=1)
            rows = start + np.flatnonzero(resolved)
            frames.append(
                pd.DataFrame(
                    {
                        "ROW": rows.astype(np.int64),
                        "GRID_LAT": (89.95 - (rows // NLON) * 0.1).astype(np.float32),
                        "GRID_LON": ((rows % NLON) * 0.1 + 0.05).astype(np.float32),
                        "N_WET": total[resolved].astype(np.uint32),
                        "POT99_MM_3H": BIN_CENTRES[indices],
                        "THRESHOLD_METHOD": "nearest-rank 500-bin histogram",
                        "IS_EXACT_TAIL_REFINEMENT": False,
                    }
                )
            )
        del counts, overflow
    if unresolved_rows:
        rows_to_refine = np.concatenate(unresolved_rows).astype(np.int64)
        counts_to_verify = np.concatenate(unresolved_counts).astype(np.uint64)
        print(
            f"Exactly refining {len(rows_to_refine)} cells without a precipitation ceiling",
            flush=True,
        )
        exact = exact_unbounded_thresholds(
            source_dir,
            override_dir,
            rows_to_refine,
            counts_to_verify,
            start_year,
            end_year,
        )
        exact["GRID_LAT"] = (89.95 - (exact["ROW"] // NLON) * 0.1).astype(np.float32)
        exact["GRID_LON"] = ((exact["ROW"] % NLON) * 0.1 + 0.05).astype(np.float32)
        exact["N_WET"] = counts_to_verify.astype(np.uint32)
        frames.append(exact)
    thresholds = (
        pd.concat(frames, ignore_index=True).sort_values("ROW").reset_index(drop=True)
    )
    if thresholds["ROW"].duplicated().any():
        raise RuntimeError("Threshold output contains duplicate canonical grid rows")
    thresholds["WET_SAMPLE_RULE"] = ">0.3 mm/3h"
    thresholds["QUANTILE_METHOD"] = (
        "nearest-rank; exact unbounded refinement when POT99 >120 mm/3h"
    )
    thresholds["BIN_WIDTH_MM_3H"] = np.float32(
        (HISTOGRAM_MAX_MM_3H - WET_MIN_MM_3H) / N_BINS
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    thresholds.to_parquet(temporary, index=False)
    temporary.replace(output)
    print(f"Wrote {len(thresholds):,} POT99 thresholds: {output}", flush=True)


def main() -> None:
    """Run or resume the annual scans in chronological order."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--override-dir", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2023)
    args = parser.parse_args()
    for year in range(args.start_year, args.end_year + 1):
        if year_is_fully_merged(args.work_dir, year):
            print(f"[{year}] cumulative blocks already include this year", flush=True)
            continue
        scan_year(args.source_dir, args.override_dir, args.work_dir, year)
        merge_year(args.work_dir, year)
    derive_thresholds(
        args.source_dir,
        args.override_dir,
        args.work_dir,
        args.output,
        args.start_year,
        args.end_year,
    )


if __name__ == "__main__":
    main()
