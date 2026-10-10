# Method

**Images.** One 16-bit TIFF per well, site and channel (default naming: ImageXpress / MetaMorph). Channels are mapped to
markers in the config (per plate if the order differs; `diffquant inspect` lists filter metadata).

**Segmentation.** Cellpose (default `cpsam`, Cellpose ≥ 4) on the cell-type marker and the nucleus channel
(H × W × 3 input: marker, nucleus, empty; automatic diameter; flow threshold 0.4, cell-probability threshold 0).
Any Cellpose model can be used; choose by visual inspection.

**Per-cell measurement.** Mean raw intensity of every channel within each mask; area, centroid, solidity, border contact.

**Object filter (optional).** Before gating, objects can be removed by area (`filters.min_area_px`, `max_area_px`) or
because they touch the image border (`exclude_border`). Off by default: Cellpose already drops objects < 15 px, and
marker-negative debris falls below the cell-type gate. Filtered objects never enter the per-plate rules
(`filtered = True` in `cells_called.csv.gz`).

**Background.** Per site and channel, the median of cell-free pixels (Euclidean distance > 15 px from any mask). A site
uses its own value only if cell-free pixels cover ≥ 20 % of the field; otherwise the plate value (median over qualifying
sites; if fewer than 10 sites qualify, the 10 sites with the largest cell-free fraction). Dense fields are excluded from
their own estimate because their "cell-free" pixels are dominated by halo and unsegmented cells. Cell mean =
max(raw mean − background, 0).

**Cell-type gate (per plate).** On v = log10(mean + 1): the main population peak = mode of a Gaussian-smoothed histogram
(470 bins over −0.1…4.6, σ = 3 bins). Cells with v > peak − Δ are of the cell type; of these, cells above
median + 4 × 1.4826 × MAD are dropped as outliers. Alternatives: Otsu or a fixed value.

**Readout (per plate, among gated cells).** Otsu threshold t on log10(mean + 1) (256 bins); a cell is positive if its mean
> factor × (10^t − 1). Alternative: a fixed value.

**Well and plate summary.** % positive = positive / gated cells per well; wells with < 200 gated cells are flagged.
With a plate map, the per-plate condition mean (valid wells) minus the matched vehicle of the same plate gives Δ;
across plates (plate = replicate): mean, SD, 95 % CI, two-sided one-sample t-test, Benjamini–Hochberg correction.

**Focus QC.** Variance of the Laplacian of the Gaussian-smoothed (σ = 1 px), 2 × down-sampled image divided by the squared
median intensity; sites < 0.5 × plate median are flagged.

**Calibration.** See [recalibration.md](recalibration.md).

**Provenance.** diffquant 0.1.0 is a re-implementation of the analysis scripts of the keratinocyte differentiation screen
(Watt lab, EMBL, 2024–2026). Validation (2026-10-10): run from raw images and the original masks on three full plates
(351,315 cells, 184 wells), diffquant reproduced the original per-cell background-subtracted means, gate and readout calls,
plate thresholds and well percentages with zero differences; `diffquant segment` (Cellpose 4.2.1.1, cpsam, GPU) reproduced
the original mask of the example field pixel for pixel (521 cells).
