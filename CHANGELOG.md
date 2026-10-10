# Changelog

## 0.1.0 (2026-10-10)
First public version, ported from the analysis scripts of the keratinocyte differentiation screen
(validated: on 3 full plates, 351,315 cells / 184 wells, per-cell values, calls, thresholds and well % are identical to the original; `segment` reproduces the original cpsam mask of the example field pixel for pixel).
- `run` / `segment` / `measure` / `threshold` / `focus` / `summarize`
- blinded scoring sets (`score-sets`) + Fiji scorer + `calibrate`
- `overlay` QC images, `inspect` channel metadata, `plates` for SLURM arrays
