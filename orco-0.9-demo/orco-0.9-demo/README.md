# ORCO 0.9 Demo

This is a functional demonstration of ORCO (Optimized Residue Classification Oracle), a lightweight residue assignment cascade system for NMR data. This demo uses manually curated residue statistics and cascaded logic based on side-chain chemical shifts (CA, CB, CG1, CG2, CE).

## Features

- Cascaded classification based on shift availability
- Probabilistic residue assignment with uncertainty thresholding
- Flexible output format (ASAP-compatible)
- Expandable scaffold for future versions

## Usage

Place your `orco_input.csv`, `stats.csv`, and `res_seq.txt` in the `data/` folder and run `main.py`.

## Version

0.9 (Demo)

## Notes

This is a scaffold release. Expect placeholder code and hardcoded thresholds.
