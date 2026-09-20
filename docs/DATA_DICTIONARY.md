# Data dictionary

CSV files are UTF-8, use a decimal point and retain computational precision.
Empty numerical fields denote unavailable values; they are not zero.
Read basin tables with `keep_default_na=False, na_values=['']` so that `NA`
remains **North Atlantic**. Time columns in the S2 example are UTC.

| Field or family | Definition / unit |
|---|---|
| `year`, `YEAR` | Calendar year, 1980–2023 |
| `domain`, `HEMISPHERE` | GLOBAL, NH, SH or basin code as applicable |
| `BASIN` | WP, EP, NA, NI, SI, SP; SA occurs in full global counts |
| `surface` | land or ocean at the precipitation grid location |
| `band_idx`, `BAND` | Zero-based band index, 0–3 |
| `band_label`, `BAND_LABEL` | Absolute storm-center latitude interval in degrees |
| `SID`, `TIME` | Storm identifier and three-hour UTC timestamp |
| `TC_LAT`, `TC_LON_360` | Storm-center coordinates in degrees; longitude [0,360) |
| `GRID_LAT`, `GRID_LON` | Precipitation grid-center coordinates in degrees |
| `ROW` | Canonical zero-based cell key on the 1,800 × 3,600 grid |
| `PRECIP_3H_MM` | Native three-hour precipitation accumulation, mm |
| `PRECIP_MMHR` | Three-hour accumulation / 3, mm h⁻¹ |
| `POT99_MM_3H` | Local wet-period 99th-percentile threshold, mm per 3 h |
| `IS_INTERPOLATED_TIMESTEP` | Whether the input field is a documented replacement |
| `TCEP_INTENSITY_MM_3H` | Mean precipitation of extreme cells in a storm–3-hour event |
| `N_EXTREME_CELLS` | Extreme-grid-record count in the event |
| `N_t`, `N_EVENTS`, `n_j`, `n_records` | Relevant event/record count; see the file's estimand |
| `vol_j`, `INTENSITY_SUM_MM_3H` | Sum of intensities, **not a physical rainwater volume** |
| `TCP_mean`, `ANNUAL_MEAN_MM_3H` | Annual mean TCEP intensity, mm per 3 h |
| `mean_lat` | Mean absolute latitude in degrees |
| `trend_per_decade`, `trend` | Trend in the corresponding variable per decade |
| `p_value`, `MANN_KENDALL_P`, `mk_p` | Two-sided Mann–Kendall p value |
| `fit`, `fit_lower`, `fit_upper` | Fitted value and bootstrap limits at a given year |
| `tot_point`, `total` | Latitude-band contribution to the total trend |
| `mig_point`, `migration`, `redistribution` | Redistribution contribution, not TC speed |
| `int_point`, `intensity`, `within_band` | Within-band intensity contribution |
| `resid_point`, `residual` | Total band slope minus redistribution and within-band slopes |
| `*_ci_low`, `*_ci_high` | 95% bootstrap interval limits in trend units |
| `REFERENCE_MM_3H` | Constant reference intensity C |
| `REDISTRIBUTION_MM_3H` | Annually aggregated R(t) |
| `WITHIN_BAND_MM_3H` | Annually aggregated W(t) |
| `INTERACTION_MM_3H` | Annually aggregated X(t) |
| `bootstrap_significant` | True only when the stored bootstrap interval excludes zero |
| `climatological_mean` | Mean of annual band means, mm per 3 h |
| `pooled_median` | Median across all events in the band, mm per 3 h |
| `valid_years` | Number of finite annual values contributing to a map-cell trend |

Precipitation trend units are mm (3 h)⁻¹ decade⁻¹. Latitude trends are degrees
decade⁻¹. Do not compare their numeric columns without checking the units.

## Supplied input groups

- **Fig01**: 132 global/hemispheric annual rows and 16,536 mapped 1° cells.
- **Fig02**: eight band rows plus hemispheric annual series.
- **Fig03/Fig04**: six basin annual series and 24 basin-band rows. NI's outer
  band has no events; the plotting code masks its stored legacy zero placeholders,
  and Table S2 records missing estimates explicitly.
- **S2**: 13,162 Walaka TCPF wet cells and 1,393 extreme cells. The event's
  mean intensity is 13.53 mm per 3 h when displayed to two decimals.
- **S3**: annual track/event absolute-latitude means in each hemisphere.
- **S4**: full band-labelled intensity samples, with 173,966 Northern Hemisphere
  and 81,805 Southern Hemisphere events, plus selected quantiles.
- **S5**: annual surface-band counts/intensity sums, surface means and band slopes.
- **Table S1 sources**: annual band values, component trend estimates and
  bootstrap intervals from `04_analyze_decomposition.py`.
- **Table S2 count source**: all basin-band event counts, including SA records
  retained for reconciliation with the global total.

The generic names `TCP_mean`, `mig`, `int`, and `vol_j` are retained in archived
CSV schemas for reproducibility; their meanings are defined above. Display
estimates use two decimal places and p values three; p < 0.001 is never written
as p = 0.000. A displayed 0.00 can represent a small nonzero estimate.

The compact Markdown tables follow the manuscript layout. The `A`, `B` and `C`
suffixes on supporting CSV filenames identify data groups, not extra manuscript
tables. Table S2 displays WNP, NIO and SIO for the archived codes WP, NI and SI.
