# TCEP latitude analysis

This repository provides the workflow for constructing a global tropical cyclone extreme precipitation (TCEP) event catalogue and analysing its intensity changes during 1980–2023.

It accompanies the manuscript **Poleward migration reinforces global tropical cyclone extreme precipitation intensification**.

Dataset construction combines three-hourly MSWEP V2.8 precipitation at 0.1° resolution with IBTrACS V4 best-track data. Tropical cyclone (TC) precipitation is identified using the tropical cyclone precipitation feature (TCPF) approach ([Jiang et al., 2011](https://doi.org/10.1175/2011JAMC2662.1); [Guzman and Jiang, 2021](https://doi.org/10.1038/s41467-021-25685-2)). Contiguous wet grid cells within 1,000 km of each storm centre form precipitation features. Features with centroids within 500 km are attributed to that storm.

TCEP is identified where precipitation attributed to a TC exceeds the local 99th percentile of all wet-period precipitation during 1980–2023. These extreme cells are grouped by storm and three-hourly timestep into TCEP events. Event intensity is the mean precipitation across the extreme cells in each event.

The analysis adapts the frequency-weighted latitude-band decomposition of [Guo and Tan (2026)](https://doi.org/10.1073/pnas.2535524123) to mean TCEP intensity. It separates annual intensity changes into latitudinal redistribution, within-latitude-band intensity change, and their interaction. The redistribution component captures changes in event shares across latitude bands at fixed climatological band intensities. The within-band component captures intensity changes at fixed climatological event shares. The interaction component captures joint changes in event shares and band intensities. Applied separately to hemispheres and ocean basins, the framework quantifies how latitudinal redistribution contributes to changes in mean TCEP intensity.

The repository includes dataset-construction and analysis code, plotting scripts, and derived CSV inputs supporting Figures 1–4, Figures S1–S5, and Tables S1–S3. The supplied inputs support reproduction of the analyses and plotted results. The documented raw-data workflow provides the processing steps for rebuilding the event catalogue from the original provider datasets.

Version 1.1.0: [GitHub release](https://github.com/kelly-Liu44/TCEP-latitude-analysis/releases/tag/v1.1.0) and [Zenodo archive](https://doi.org/10.5281/zenodo.22857802).

## Start here

Python 3.10 or newer is required; the local validation uses Python 3.11.
Create and activate a virtual environment, then run from the repository root:

```bash
python -m pip install -e ".[figures,test]"
python -m pytest -q
python scripts/verify_release.py
python figures/run_all_figures.py
python scripts/06_build_supplementary_tables.py
```

These commands use the included, full-precision CSV inputs. They do not require
the raw precipitation archive. Figures are regenerated in their numbered
`figures/` folders, and tables in `tables/`. The package contains one image per
manuscript figure; the default commands replace these outputs in place.
Cartopy may download public Natural Earth coastline files on first use.

Figure S1 is an author-supplied JPEG illustration; its original drawing code is
unavailable. The other eight figures are supplied and regenerated as 600-dpi
PNGs using the bundled Lato fonts and the manuscript's plotting parameters.
To keep a separate output set or export other formats, specify these options:

```bash
python figures/run_all_figures.py --output-dir /path/to/figure-exports --formats pdf svg tif
python scripts/06_build_supplementary_tables.py --output-dir /path/to/table-exports
```

To render one figure:

```bash
python figures/run_all_figures.py --only Fig02
```

## Repository structure

```text
src/tcep/          reusable geometry, attribution, events and statistical methods
scripts/           numbered processing stages, audits and release verification
figures/           one image per figure, plotting scripts, captions, CSVs and fonts
tables/            full-precision table sources and readable reference tables
tests/             scientific unit tests and numerical regression checks
docs/              methods, data dictionary and reproduction workflow
metadata/          source hashes, expected numerical anchors and file checksums
.github/workflows/  automated tests for future pushes and pull requests
```

## Figures and tables

| Manuscript item | Folder | Content |
|---|---|---|
| Figure 1 | `figures/Fig01_global_trends` | Global/hemispheric intensity and mapped trends |
| Figure 2 | `figures/Fig02_latitude_decomposition` | Hemispheric latitude-band contributions |
| Figure 3 | `figures/Fig03_basin_intensity_trends` | Six basin intensity trends |
| Figure 4 | `figures/Fig04_basin_decomposition` | Basin latitude-band contributions |
| Figure S1 | `figures/FigS01_TCEP_definition` | Conceptual definition, original image |
| Figure S2 | `figures/FigS02_TCEP_identification` | MSWEP V2.8 Walaka identification example |
| Figure S3 | `figures/FigS03_poleward_shifts` | Track/event mean absolute latitude |
| Figure S4 | `figures/FigS04_intensity_distributions` | Conditional intensity distributions |
| Figure S5 | `figures/FigS05_land_ocean` | Land–ocean trends and decomposition |
| Table S1 | `tables/TableS01_latitude_decomposition` | Hemispheric contributions by latitude band |
| Table S2 | `tables/TableS02_basin_decomposition` | Basin-band components |
| Table S3 | `tables/TableS03_land_ocean` | Surface-specific components |

## Reproduce from the input datasets

Obtain the historical **MSWEP V2.8** three-hourly 0.1° product and **IBTrACS V4**
from their providers. A newer MSWEP version or a revised track archive does not
constitute the same input dataset. The original analysis uses the IBTrACS
`v04r01` CSV. Raw archives, complete TCPF/extreme-cell catalogues and global
threshold fields are not redistributed in this repository.

[Full workflow](docs/WORKFLOW.md) gives the exact stage commands, input filename
conventions, nine documented replacement fields, intermediate filenames,
resource requirements and instructions for rendering independently generated
inputs. The processing entry point is:

```bash
python scripts/run_workflow.py --help
```

The complete raw-field calculation is storage- and computation-intensive.
It does not run as part of the quick start or automated tests. The included tests and `scripts/verify_release.py` check the supplied
inputs and numerical anchors without running the raw-field workflow.

## Data access, licenses and citation

- MSWEP V2.8: <https://www.gloh2o.org/mswep/>; methodology:
  [Beck et al. (2019)](https://doi.org/10.1175/BAMS-D-17-0138.1).
- IBTrACS V4: NOAA/NCEI,
  <https://doi.org/10.25921/82ty-9e16>;
  [Knapp et al. (2010)](https://doi.org/10.1175/2009BAMS2755.1).
- Original source code: [MIT](LICENSE).
- Precipitation-derived CSVs and figure products: see
  [DATA_LICENSE.md](DATA_LICENSE.md). The code license does not override MSWEP's
  non-commercial data terms.

Use [CITATION.cff](CITATION.cff) for software authorship and version information.

To build a ZIP of the public source tree, run `python scripts/build_release.py`.
The archive and its SHA256 digest are written to `dist/`; generated results,
raw archives, caches and local environment files are excluded. Creating this
archive does not publish it.

See [METHODS.md](docs/METHODS.md) and [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md)
for the implemented methods and source-data fields.
