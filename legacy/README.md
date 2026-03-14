# Legacy: MATLAB implementation

Original MATLAB implementation and experiment scripts from the **NeurIPS 2020** paper. Kept for reference and reproducibility; the main, optimized implementation is at the **repo root** (`surf/`, `visualize.py`, etc.).

## Requirements

- MATLAB (R2016a or later)

## Usage

From the project root:

```matlab
addpath('legacy');
[I, koi] = surf(sort(rand(1, 2^8 - 1)), 0.25, 1);
```

## Experiment scripts

| File | Description |
|------|-------------|
| `synthetic_experiments.m` | Synthetic data (e.g. beta mixtures). |
| `salary_experiments.m`    | Salary data. Requires `../data/salary.dat`. |
| `sensor_experiments.m`   | Sensor data. Requires `../data/sensor1.dat`. |
| `cover_experiments.m`    | Cover-type data. Requires `../data/covt1.dat`. |

Run from project root: `addpath('legacy'); run('legacy/salary_experiments.m');` etc.

## Note

This implementation is recursive and not optimized for speed. For faster runs and easier integration, use the Python implementation at the repo root (`pip install -e .` then `python visualize.py` etc.).
