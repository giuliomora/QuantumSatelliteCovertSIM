# Quantum Doppler Covert Inter-Satellite Link Simulator

MATLAB simulator for the Master's thesis *Feasibility Study of a Quantum Doppler Effect-Based Channel for Covert Inter-Satellite Communications* (Computer Science, Università di Padova, A.Y. 2025–2026).

The simulator estimates how much information two satellites in low Earth orbit can exchange over one pass before an eavesdropper notices that a transmission is taking place (the *covert volume*).

Each bit is the sign of a spectral fringe written on the quantum Doppler signature of a ²⁷Al⁺ emitter: a superposition of two momentum packets whose relative phase (0 or π) carries the bit, and whose spatial displacement Δx turns the interference term into a regular fringe. The fringe is imprinted on a laser probe and read by a matched-filter receiver. Its pitch changes at every symbol under a shared key (*fringe hopping*), so an eavesdropper who knows the scheme but not the key cannot use a single filter.

## What it simulates

The script runs nine blocks in sequence:

| Block | Content |
|---|---|
| 1 | Physical constants |
| 2 | Orbital geometry: a cross-plane pass (800 km, planes 60° apart, closest approach about 300 km) propagated with SGP4; distance, classical Doppler shift, line-of-sight and solar-exclusion constraints |
| 3 | Optical link budget, with pointing loss at both terminals |
| 4 | Symbol spectrum, choice of the operating point that maximises the covert volume, hopping alphabet, adversary models and run-time checks |
| 5 | Campaign A: simulated bit-error rate against the theoretical curve |
| 6 | Campaign B: bit-error rate along the pass, with pointing fades, Doppler residual and calibration errors |
| 7 | Sensitivity to the declared inputs: calibration errors, optical depth, branch displacement |
| 8 | Capacity and covert volume against three eavesdroppers, on a grid of detection thresholds, and a final summary |
| 9 | The seven figures of the results chapter of the thesis |

The three eavesdroppers are:

- **Eve-A**: energy detector, knows neither the scheme nor the key;
- **Eve-B**: knows the key (key-compromise case);
- **Eve-C**: knows the scheme but not the key; the reference adversary when hopping is on.

## What it does not simulate

- **The atomic source:** described by three declared inputs (packet width Δ, displacement Δx, optical depth OD); its preparation and coherence are assumed.
- **The Doppler compensation:** the classical shift is computed but assumed removed; only the residual offset is simulated.
- **The cipher:** the keyed pitch sequence is not generated; each pitch is drawn uniformly among the admissible ones, as a secure cipher looks to anyone without the key.
- **The channel code:** each symbol carries its capacity.

The full list of assumptions and parameters is in the appendices of the thesis.

## Requirements

- MATLAB R2021a or later
- Satellite Communications Toolbox (`satelliteScenario`, SGP4 propagation)
- Phased Array System Toolbox (`fspl`)
- Signal Processing Toolbox (`findpeaks`)

## Usage

Run the script from MATLAB:

```matlab
CovertSimNew
```

All results are printed to the console, block by block; the summary at the end of block 8 collects the final figures. Block 9 opens the figures but does not save them: the file name for each is printed at the end of the run.

Main switches (block 4):

| Setting | Default | Effect |
|---|---|---|
| `cfg.hopping` | `true` | `false` gives the fixed-pitch scheme, for comparison |
| `cfg.fastMC` | `true` | `false` runs the full spectral Monte Carlo instead of the sufficient-statistic path (much slower; used as a check) |
| `cfg.profileFit` | `true` | `false` removes Bob's per-frame estimation of the classical profile |

## Reproducibility

Every random number comes from a single seed (`SEED_MASTER = 42`), and each Monte Carlo campaign, point and realisation uses its own random substream: the same seed gives the same results, and each run can be repeated on its own.

With the default settings, the summary gives a covert volume against Eve-C of about 4.3 kbit at the detection threshold D = 3, and about 1.4 kbit at D = 1.

## Physical model

The quantum Doppler signature is taken from:

P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Dębski, "Quantum time dilation in atomic spectra", *Phys. Rev. Research* 3, 023053 (2021).

The spatial displacement between the branches, the fringe hopping and the covert-communication analysis are added in this work.

## Author

Giulio Morandini — Università di Padova