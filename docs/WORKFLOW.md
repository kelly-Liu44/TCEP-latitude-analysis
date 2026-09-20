# Reproduction workflow

Run commands from the repository root after installation. The supplied figure
CSVs support the quick start in README. This document describes regeneration
from externally obtained input datasets and saved analysis archives.

## 1. Supply the original input versions

MSWEP V2.8 files use `YYYYDDD.HH.nc`, where DDD is day of year and HH is a
three-hour UTC tick. Put the files in one directory. Required variable:
`precipitation`; dimensions `lat, lon`; native units `mm/3h`; the global grid has
1,800 latitude rows and 3,600 longitude columns. The raw source longitude axis
is reordered to [0,360) internally.

Use the all-basins IBTrACS V4 `v04r01` CSV with its header and units row. The
method paper is Knapp et al. (2010); it does not replace the dataset citation.

## 2. Handle missing or unreadable precipitation fields

The manuscript's input archive uses nine documented replacements: seven fields
on 31 December 2020 (03 through 21 UTC), 30 September 2011 03 UTC, and 12 May
2023 21 UTC. They are linearly interpolated between their valid bracketing
fields and stored separately. Do not overwrite the provider archive.

```bash
python scripts/00b_interpolate_missing_2020.py --source-dir /path/to/MSWEP --output-dir /path/to/archive/overrides
python scripts/00c_reconstruct_corrupt_2011.py --source-dir /path/to/MSWEP --output-dir /path/to/archive/overrides
python scripts/00d_reconstruct_corrupt_field.py --source-dir /path/to/MSWEP --output-dir /path/to/archive/overrides --before-name 2023132.18.nc --target-name 2023132.21.nc --after-name 2023133.00.nc --before-time 2023-05-12T18:00:00 --target-time 2023-05-12T21:00:00 --after-time 2023-05-13T00:00:00
python scripts/00_validate_inputs.py --mswep-dir /path/to/MSWEP --override-dir /path/to/archive/overrides --ibtracs-csv /path/to/ibtracs.ALL.list.v04r01.csv
```

Replacement scripts verify the bracketing data and mark provenance as
`IS_INTERPOLATED_TIMESTEP` in downstream records. These are archive-specific
preprocessing choices, not a claim that every provider download has these gaps.

## 3. Run the scientific stages

The examples below use `/path/to/archive` as a user-chosen work directory.
On Windows, use quoted native paths in place of the examples. Each command
stops on a failed subprocess. Do not use a threshold work directory from a
different period or dataset version.

```bash
python scripts/run_workflow.py --stage tracks --work-dir /path/to/archive --ibtracs-csv /path/to/ibtracs.ALL.list.v04r01.csv
python scripts/run_workflow.py --stage tcpf --work-dir /path/to/archive --mswep-dir /path/to/MSWEP --override-dir /path/to/archive/overrides
python scripts/run_workflow.py --stage thresholds --work-dir /path/to/archive --mswep-dir /path/to/MSWEP --override-dir /path/to/archive/overrides
python scripts/run_workflow.py --stage events --work-dir /path/to/archive
python scripts/run_workflow.py --stage analysis --work-dir /path/to/archive
```

The default period is 1980–2023. `--stage all` explicitly runs all these stages.
The TCPF step rejects an existing final annual file unless its individual CLI
receives `--overwrite`. The histogram step resumes its chronological blocks.
Other processing stages regenerate their specified outputs; choose a fresh work
directory when preserving a previous run.

Expected archive structure:

```text
data/tracks_processed.parquet
outputs/tcpf/tc_gridded_YYYY.parquet
outputs/thresholds/pot99_thresholds.parquet
outputs/events/tcep_events_YYYY_pot99.parquet
outputs/extreme_cells/tcep_extreme_cells_YYYY_pot99.parquet
outputs/components/annual_bands_NH.csv, annual_bands_SH.csv, ...
outputs/figure_inputs/Fig*/data/*.csv
outputs/tables/TableS*/...
```

The histogram grid contains 6,480,000 cells × 500 bins. One uint32 cumulative
histogram is about 13 GB (decimal), and a uint16 annual histogram about 6.5 GB,
before block temporaries and other files. The raw fields and TCPF archive
require additional, substantially larger storage. Processing streams annual
Parquet records, but memory and runtime depend on input layout and hardware.
No universal runtime or whole-archive storage minimum is asserted.

## 4. Render regenerated figure inputs

The S2 example is extracted separately from the saved 2018 archive:

```bash
python figures/FigS02_TCEP_identification/prepare_FigS02.py --project-dir /path/to/archive --output-dir /path/to/archive/outputs/figure_inputs/FigS02_TCEP_identification
python figures/run_all_figures.py --data-root /path/to/archive/outputs/figure_inputs --output-dir /path/to/archive/outputs/rendered_figures
```

S1 always uses the supplied original illustration. All other figures use the
requested data root. Plotting scripts do not recalculate trends or significance.
Without `--output-dir`, each figure replaces its supplied image in its numbered
`figures/` folder. The default is PNG only; use `--formats pdf svg tif` for
additional formats. The bundled Lato fonts reproduce the manuscript typography.
To regenerate S5 summaries from its compact annual source table alone:

```bash
python figures/FigS05_land_ocean/prepare_FigS05.py
```

`python scripts/06_build_supplementary_tables.py` rebuilds the three compact
Markdown tables and their full-precision supporting CSVs in `tables/`. Table S1
has two hemisphere rows; Table S2 has six basin rows; their cells list
redistribution / within-band / residual contributions. Table S3 has four
surface–hemisphere rows and reports total trends with aggregated
redistribution / within-band / interaction components. Use `--output-dir` to
write a separate table set. See `METHODS.md` for the residual–interaction
distinction and Theil–Sen non-additivity.

## 5. Check the results

Use the `02b`/`02c`/`02d` scripts for TCPF geometry, provenance and archive checks;
`03a`/`03b` for thresholds and event/cell consistency. Every script has `--help`.
The threshold audit's default 92 exact-tail cells applies to the manuscript's
archive; set the expected count explicitly for a different input version.

`scripts/verify_release.py` checks the supplied source-data hashes, core numerical
anchors, latitude decomposition and table consistency. Tests do not fetch raw
MSWEP or start a 44-year calculation. Add `--bootstrap` to repeat the aggregate
component interval calculations. Dependency versions used for local numerical
validation are listed in `requirements-tested.txt`.
