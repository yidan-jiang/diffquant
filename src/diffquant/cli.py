"""Command-line interface:  diffquant <command> config.yaml [options]

  init          write a commented config template
  plates        list plate names (one per line; handy for SLURM array jobs)
  inspect       show the channel tokens of each plate with their wavelength/filter metadata
  segment       Cellpose segmentation (skipped if segmentation.masks is set)
  measure       per-cell means + background subtraction
  threshold     cell-type gate + readout cut per plate, % positive per well
  focus         focus QC per site
  summarize     collect all plates; condition summary and delta vs vehicle if a plate map is given
  run           segment -> measure -> threshold -> focus -> summarize
  score-sets    blinded single-cell crops for manual scoring in Fiji (fiji/score_cells.py)
  calibrate     choose gate delta and readout factor from the manual calls
  overlay       ImageJ TIFFs with outlines coloured by call, for a visual check
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from . import config as C


def _plates(cfg, names):
    if not names:
        return cfg["plates"]
    sel = [p for p in cfg["plates"] if p["name"] in names]
    missing = set(names) - {p["name"] for p in sel}
    if missing:
        raise SystemExit(f"unknown plate(s): {sorted(missing)}")
    return sel


def main(argv=None):
    ap = argparse.ArgumentParser(prog="diffquant", description="Cell-type-gated quantification of marker-positive cells in plate images.")
    ap.add_argument("--version", action="version", version=f"diffquant {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init", help="write a config template")
    p.add_argument("path", nargs="?", default="config.yaml")
    for name in ("plates", "inspect", "segment", "measure", "threshold", "focus", "summarize", "run", "score-sets", "calibrate", "overlay"):
        p = sub.add_parser(name)
        p.add_argument("config")
        p.add_argument("--plate", action="append", help="restrict to this plate (repeatable)")
        p.add_argument("--jobs", type=int, help="parallel processes (overrides config)")
        if name in ("segment", "run"):
            p.add_argument("--overwrite", action="store_true", help="re-segment sites that already have masks")
        if name == "score-sets":
            p.add_argument("--per-decile", type=int, default=6)
            p.add_argument("--n-random", type=int, default=40)
            p.add_argument("--per-band", type=int, default=12)
        if name == "calibrate":
            p.add_argument("--scorer", help="use only calls_<scorer>.csv")
            p.add_argument("--exclude", nargs="*", default=[], help="plates to leave out")
        if name == "overlay":
            p.add_argument("--n-sites", type=int, default=3)
            p.add_argument("--site", nargs=2, action="append", metavar=("WELL", "SITE"))
            p.add_argument("--marker")
    a = ap.parse_args(argv)
    if a.cmd == "init":
        if os.path.exists(a.path):
            raise SystemExit(f"{a.path} exists")
        open(a.path, "w").write(C.TEMPLATE)
        print(f"wrote {a.path}")
        return 0
    cfg = C.load(a.config)
    if getattr(a, "jobs", None):
        cfg["jobs"] = a.jobs
    plates = _plates(cfg, a.plate)
    if a.cmd == "plates":
        print("\n".join(p["name"] for p in plates))
        return 0
    if a.cmd == "inspect":
        from . import io
        for p in plates:
            print(f"{p['name']}  ({p['path']})  configured: {io.channels(cfg, p)}")
            print("\n".join("  " + l for l in io.inspect(p, cfg)))
        return 0
    os.makedirs(cfg["output"], exist_ok=True)
    from . import calibrate, focus, measure, overlay, scoring, segment, summarize, threshold
    if a.cmd == "calibrate":
        calibrate.run(cfg, scorer=a.scorer, exclude=a.exclude)
        return 0
    if a.cmd == "summarize":
        summarize.run(cfg)
        return 0
    for pl in plates:
        C.plate_out(cfg, pl["name"])
        if a.cmd in ("segment", "run"):
            segment.run(cfg, pl, overwrite=a.overwrite)
        if a.cmd in ("measure", "run"):
            measure.run(cfg, pl)
        if a.cmd in ("threshold", "run"):
            threshold.run(cfg, pl)
        if a.cmd == "focus" or (a.cmd == "run" and cfg["focus"].get("enabled", True)):
            focus.run(cfg, pl)
        if a.cmd == "score-sets":
            scoring.run(cfg, pl, per_decile=a.per_decile, n_random=a.n_random, per_band=a.per_band)
        if a.cmd == "overlay":
            overlay.run(cfg, pl, n_sites=a.n_sites, sites=a.site, marker=a.marker)
    if a.cmd == "run":
        summarize.run(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
