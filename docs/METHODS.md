# Implemented methods

## Input and event definition

The analysis covers 1980–2023. MSWEP V2.8 provides three-hour precipitation
accumulations on a 0.1° grid. Stored native values are **mm per 3 h**; division
by three gives mm h⁻¹. The wet threshold is strictly **>0.30 mm per 3 h**.

IBTrACS spur tracks are excluded. There is no wind-speed filter and no
storm-latitude cutoff. Coordinates at 00/06/12/18 UTC are linearly interpolated
to three-hour timestamps. Longitudes are unwrapped before interpolation and
then converted to [0, 360). Rare systems with no six-hour anchors retain their
valid native three-hour records. Interpolation precedes analysis-period clipping.

TCPF searches within 1,000 km of the storm center, labels wet cells with
four-neighbor connectivity, and retains features whose geometric centroids are
within 500 km. Each qualifying feature is attributed independently to a storm.
Coincident storms can therefore associate records with the same precipitation
grid cell; the data are storm-associated records, not a deduplicated global
rainfall-volume budget.

Local POT99 is estimated from **all wet MSWEP periods**, not only TC-associated
periods. The production algorithm uses 500 equal-width bins from 0.30 to
120 mm per 3 h (width 0.2394 mm). It locates the nearest-rank 99th percentile
and uses the bin center for in-range ranks. If the required rank exceeds the
histogram range, that cell is reread and its unbounded nearest-rank percentile
is evaluated exactly. Overflow values are retained in the sample count. Thus
the in-range calculation is histogram-based; it is not an exact order-statistic
calculation at every grid cell.

TCEP records strictly exceed their local threshold. All qualifying cells for
one storm and one three-hour timestamp form a TCEP event. Event intensity is
the arithmetic mean of their precipitation values. The annual main-analysis
mean weights events equally, regardless of their area or extreme-cell count.

## Latitude-band decomposition

Let n(j,t) denote the event count, alpha(j,t) its share of the hemispheric count,
and mu(j,t) the event-mean intensity for band j in year t. Bars denote means of
annual values over the analysis period. The bands are:

| Hemisphere | Absolute storm-center latitude bands |
|---|---|
| Northern | [0,15), [15,25), [25,35), [35,90] degrees |
| Southern | [0,10), [10,15), [15,20), [20,90] degrees |

The exact annual identity is:

```text
I(t) = C + R(t) + W(t) + X(t)
C    = sum_j alpha_bar(j) * mu_bar(j)
R(t) = sum_j (alpha(j,t) - alpha_bar(j)) * mu_bar(j)
W(t) = sum_j alpha_bar(j) * (mu(j,t) - mu_bar(j))
X(t) = sum_j (alpha(j,t) - alpha_bar(j)) * (mu(j,t) - mu_bar(j))
```

The annual component series are aggregated across bands **before** estimating
their Theil–Sen slopes. Multiplication by 10 converts annual slopes to decadal
slopes. The estimator is a median of pairwise slopes, so separate component
slopes do not necessarily sum to the total slope. Similarly, the sum of band
slopes is not a substitute for a hemispheric component slope.

Figures 2 and 4 show each band's total, redistribution and within-band slopes.
Their **residual** equals the total band slope minus its redistribution and
within-band slopes. This includes annual interaction and Sen non-additivity.
Only the underlying annual component X(t) is called **interaction**.

A band with no events in a particular year has zero relative frequency. For
annual identities, missing intensity is filled with the climatological band
mean. Band intensity slope estimation omits empty years. North Indian Ocean
35–90°N has no events throughout the record and is displayed as unavailable.
The six basin analyses contain 255,663 events; 108 South Atlantic events are
also present in the global and hemispheric analyses (255,771 in total).

## Trend tests and uncertainty

Total intensity and latitude trends use the two-sided Mann–Kendall test with
continuity and tie correction. There is no autocorrelation correction in this
MK calculation and no spatial multiple-testing correction for map dots.

Uncertainty uses 1,000 residual moving-block bootstrap replicates, four-year
overlapping blocks and seed 42. Hemispheric component intervals jointly
resample the detrended annual component series. Band intervals instead jointly
resample detrended event counts and intensity sums, clip reconstructed negative
counts/sums at zero, and then repeat the decomposition. Percentile intervals
use the 2.5th and 97.5th percentiles. A component is classified as significant
when that interval excludes zero. MK and bootstrap significance need not agree.

The aggregate component implementation uses SciPy's default Sen intercept;
publication fitted lines use median(y - slope*t). These intercept conventions
produce the same point slopes; preserve their respective implementations when
reproducing the archived confidence bands. All bootstrap error bars use their
actual endpoints, including percentile intervals that exclude the point estimate.

## Spatial and land–ocean analyses

Figure 1d estimates trends in annual extreme-record means within fixed 1° cells
with at least 10 valid annual values. Figure 1e shows the longitude-wise spatial
mean and standard deviation of those mapped slopes, not uncertainty of the
global event-mean trend. Hotspot rectangles are descriptive selections; the
box-average slopes are not independent regional significance tests.

For S5, surface type follows precipitation-grid location using `global-land-mask`;
hemisphere and latitude band still follow the storm center. All extreme-grid
records have equal weight. This analysis is not a partition of the event-weighted
main trend. Component confidence intervals are not estimated for this surface
analysis. The intermediate float32 intensity sums are serialized to CSV and
read as float64 before downstream calculations, matching the archived workflow.

Figure S4 uses every event in each band for the histogram/KDE normalization.
The display ends at the hemispheric 99th percentile; the visible portion is
not renormalized. Climatological means are means of annual band means, whereas
pooled medians summarize the complete event sample. Neither statistic measures
event occurrence probability or, by itself, identifies a physical mechanism.
