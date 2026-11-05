# QuantumSatelliteCovertSIM

Simulator of a covert optical link between two LEO satellites based on the quantum Doppler effect.

Code of the master thesis *Feasibility Study of a Quantum Doppler Effect-Based Channel for Covert Inter-Satellite Communications*, Giulio Morandini, Università degli Studi di Padova, 2026.

## What it does

The bit is written in the sign of a spectral fringe produced by an atomic superposition of two momentum packets and carried by a laser probe across a LEO–LEO fly-by. The program follows the whole chain:

- orbital geometry of a cross-plane pass (SGP4) and classical Doppler shift;
- free-space optical link budget with pointing loss;
- spectrum of the symbol, fringe hopping driven by a key, matched-filter receiver;
- Bob's error rate by Monte Carlo (Campaign A: validation against theory; Campaign B: along the pass);
- three adversaries (Eve-A energy detector, Eve-B scheme-aware without the key, Eve-C with the key) and the covert volume each one allows before detection;
- sensitivity analyses, control cases and the figures of the thesis.


## Requirements

MATLAB with:

- Satellite Communications Toolbox (`satelliteScenario`, `satellite`, `states`)
- Phased Array System Toolbox (`fspl`)
- Signal Processing Toolbox (`findpeaks`)

## Files

| File | Content |
|---|---|
| `CovertSimFix.m` | Main script: the reference run (seed 42) and all the figures of the thesis |
| `runSeed.m` | The part of `CovertSimFix.m` that produces the covert volumes, as a function of the seed |
| `seedStudy.m` | Variability study over seeds 1–30; writes `seedStudy.mat` |

## How to run

**Reference run.** Put all the files in one folder and run

```matlab
CovertSimFix
```

The results are printed block by block in the command window; the figures are drawn but not saved, and the names to export them with are printed at the end of block 9.

**Variability over seeds.** Run

```matlab
seedStudy
```

It first checks that `runSeed(42)` reproduces the reference run and stops if it does not; it then runs seeds 1–30, prints mean, standard deviation and range of every result, and saves `seedStudy.mat`. If `seedStudy.mat` is in the same folder, a later run of `CovertSimFix` adds the variability bands to Fig 2b.

## Switches

Set in block 4 of `CovertSimFix.m`:

- `cfg.hopping = false`: fixed-pitch baseline, without fringe hopping;
- `cfg.fastMC = false`: full spectral Monte Carlo instead of the sufficient-statistic shortcut (much slower; used in the thesis as a check);
- `cfg.profileFit = false`: receiver without the per-frame estimation of the classical profile.

## Reproducibility

Every random number comes from the master seed `SEED_MASTER = 42`, through independent streams for each campaign, point and trial.

## Physical basis

P. T. Grochowski, A. R. H. Smith, A. Dragan, K. Dębski, "Quantum time dilation in atomic spectra", *Physical Review Research* 3, 023053 (2021). Equation numbers quoted in the code refer to that work.

## Citation

If you use this code, please cite the thesis and the archived version:

> G. Morandini, *QuantumSatelliteCovertSIM*, version 1.0, 2026.