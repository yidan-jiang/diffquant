# Instructions for AI coding assistants (Codex, Claude Code, ...)

diffquant quantifies the fraction of marker-positive cells among cells of one type (gate) in high-content plate images.
Pipeline: segment (Cellpose) → measure (per-cell means, background) → threshold (gate + readout cut per plate) → focus QC →
summarize. Calibration: score-sets → manual scoring in Fiji (by a person) → calibrate.

## Layout
- `src/diffquant/` — package; one module per step (`segment`, `measure`, `threshold`, `focus`, `summarize`, `scoring`,
  `calibrate`, `overlay`), configuration in `config.py`, CLI in `cli.py`.
- `examples/config.yaml` + `examples/example_field/` — one bundled field with precomputed masks.
- `fiji/score_cells.py` — Jython script for blinded manual scoring (runs inside Fiji, not in Python).
- `docs/` — method, calibration guide, outputs, SLURM; `tests/` — pytest suite with expected example outputs.

## Commands
```bash
pip install -e ".[dev]"                     # add [seg] for Cellpose
pytest -q                                   # must pass before any commit
diffquant run examples/config.yaml          # expected: 23.9 % IVL+ in results/example/example/wells.csv
diffquant init screen.yaml && diffquant inspect screen.yaml
diffquant run screen.yaml --plate <name>    # one plate first
diffquant overlay screen.yaml --plate <name> --n-sites 3
diffquant score-sets screen.yaml            # then a person scores in Fiji
diffquant calibrate screen.yaml
```

## Setting up a new dataset
1. `diffquant init`, then `diffquant inspect` to map channel tokens (w1, w2, ...) to markers via their filter metadata;
   set `channels:` per plate if the order differs.
2. Run one plate; report `thresholds.csv`, `thresholds.png` and `wells.csv`, and write overlays for visual inspection.
3. Generate scoring sets and ask the user to score them in Fiji (`fiji/score_cells.py`). Do not create, simulate or edit
   manual calls (`calls_*.csv`).
4. Run `diffquant calibrate`; report the recommended gate delta and readout factor with their curves, and propose the config
   change. Apply it only after the user agrees, then re-run `threshold` and `summarize`.

## Rules
- The default cut-offs (gate delta 0.8, readout factor 1.3) are calibrated for one set-up only. Never present results on a
  new dataset as final without a calibration, and never tune cut-offs to obtain an expected biological result.
- Keep numerical behaviour unchanged unless asked: the example outputs in `tests/expected/` are regression references;
  update them only for an intended change and state why.
- Do not commit `results/`, images other than the bundled example, or data from users' plates.
- Keep the code compatible with Python ≥ 3.9 and add a test for every new rule or option.
