# Output reference

`<output>/<plate>/`
| file | content |
|---|---|
| `masks/<well>_s<site>_masks.tif` | label image per site (0 = background), from `segment` |
| `site_background.csv` | per site: cell-free fraction, cell-free median per channel, background used, source (site / plate_median) |
| `cells.csv.gz` | per cell: `well, site, label, area_px, centroid_y/x, solidity, touches_border`, per channel `<M>_mean_raw`, `<M>_background`, `<M>_mean` (= max(raw − background, 0)) |
| `cells_called.csv.gz` | `cells.csv.gz` + `filtered` (removed by `filters:`), `gated` (cell type), `gate_low` (below the gate), `<readout>_pos` |
| `thresholds.csv` | per plate: `gate_peak`, `gate_cut`, `gate_upper`, `frac_gated`; per readout `<M>_otsu`, `<M>_cut`, `<M>_frac_pos`, `<M>_check_histogram` (positive fraction < 2 % or > 95 %: inspect) (intensities background-subtracted, linear) |
| `thresholds.png`, `histograms.npz` | log10 histograms with the cuts |
| `wells.csv` | per well: `n_cells, n_sites, n_filtered, n_gated, n_gate_low, frac_gate_low, n_<M>_pos, pct_<M>_pos, valid` |
| `focus.csv` | per site: `focus_<M>`, relative to the plate median, `focus_flag` (< 0.5 × median) |
| `overlays/overlay_<well>_s<site>.tif` | ImageJ TIFF, outlines: green readout+, magenta readout−, cyan excluded (below gate or filtered) |

`<output>/`
| file | content |
|---|---|
| `wells_all.csv`, `thresholds_all.csv` | all plates (+ plate map columns, + number of focus-flagged sites per well) |
| `conditions_per_plate.csv` | plate map given: mean % positive of valid wells per plate and condition |
| `delta_per_plate.csv`, `delta_summary.csv` | plate map with `vehicle` column: Δ vs vehicle per plate; across plates mean, SD, 95 % CI, t-test p, BH q |
| `scoring/<plate>_<marker>/` | `crops.tif`, `index.csv` (hidden key), `display.txt`, `calls_<scorer>.csv` |
| `calibration/` | `calls.csv`, `<M>_factor_curve.csv`, `<M>_per_plate.csv`, `gate_agreement_curve.csv`, `gate_per_plate.csv`, `calibration_summary.json`, `calibration.png` |

Plate map (optional, `plate_map:` in the config): CSV with `plate, well, condition` and optionally `vehicle` (the
condition on the same plate to compare with; empty for vehicle wells).
