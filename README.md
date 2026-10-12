# Quantum Doppler Covert Inter-Satellite Link Simulator

MATLAB simulator for the Master's thesis *Feasibility Study of a Quantum Doppler Effect-Based Channel for Covert Inter-Satellite Communications* (Computer Science, Università di Padova, A.Y. 2025–2026).

It models a covert optical link between two LEO satellites. Each bit is written as a spectral fringe on the quantum Doppler signature of an Al⁺ emitter, read by a matched-filter receiver, and evaluated against three eavesdroppers over a full orbital pass.

## What it does

1. Propagates a LEO–LEO cross-plane pass (800 km, planes 60° apart) with SGP4.
2. Computes the optical link budget, including pointing loss.
3. Builds the symbol spectrum and picks the operating point that maximises the covert volume.
4. Runs Monte Carlo campaigns for Bob's bit-error rate, both against theory and along the pass.
5. Computes the covert volume against Eve-A (energy detector), Eve-B (key-aware) and Eve-C (scheme-aware, key-unaware).
6. Prints a summary and draws four figures.

## Requirements

- MATLAB R2021a or later
- Satellite Communications Toolbox
- Phased Array System Toolbox
- Signal Processing Toolbox

## Usage

Run the script from MATLAB:

```matlab
CovertSimFix
```

All results are printed to the console. The final summary block holds the headline figures.

Two switches in Section 4 change the run:

| Setting | Default | Effect |
|---|---|---|
| `cfg.hopping` | `true` | `false` gives the fixed-pitch baseline |
| `cfg.fastMC` | `true` | `false` runs the full spectral Monte Carlo (slower, used as a check) |

Results are reproducible: every random number comes from a single seed (`SEED_MASTER = 42`).

## Physical model

The physical model is based on:

P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Dębski, "Quantum time dilation in atomic spectra", *Phys. Rev. Research* 3, 023053 (2021).

## Author

Giulio Morandini — Università di Padova
