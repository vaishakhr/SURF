# SURF

SURF learns a **piecewise polynomial density** from univariate samples on [0,1]. This repo provides a **Python implementation** (repo root) and the original **MATLAB** code in **`legacy/`** (NeurIPS 2020).

---

## Overview

SURF takes a set of samples and returns a piecewise polynomial density estimator over the unit interval. The algorithm recursively merges intervals according to a criterion controlled by a tuning parameter α.

**How SURF works (high level):**
1. **Input:** Sorted samples in (0, 1), with total count `n - 1` and `n = 2^k` (power of two).
2. **Initial partition:** Boundaries are 0, the sorted sample values, and 1; each interval gets a polynomial fit from the empirical distribution.
3. **Merge step:** Pairs of adjacent intervals are merged when a penalized L1 criterion is ≤ 0; the penalty scales like α √(log n / n). Larger α → more merging → fewer pieces.
4. **Output:** `boundaries` (breakpoints) and `piece_coeffs` (polynomial coefficients per piece). Evaluate with `regpoly(x, piece_coeffs[i])` for x in `[boundaries[i], boundaries[i+1]]`.

**Notation (used in code and docs):** `alpha` = tuning parameter; `degree` = polynomial degree per piece (1–4); `boundaries` = 1D array of breakpoints including 0 and 1; `piece_coeffs` = shape `(num_pieces, degree+1)` (constant term first); `n` = number of quantile points = `len(samples) + 1`, must be a power of 2.

---

## Quick Start (Python)

From the project root:

```bash
pip install -e .
python visualize.py
```

With data in `data/`:

```bash
python visualize_salary.py   # uses data/salary.dat
python visualize_sensor.py   # uses data/sensor1.dat
```

Synthetic experiments (beta, gamma, Gaussian mixtures; L1 error and plots):

```bash
python synthetic_experiments.py
```

Plots are saved in the repo root (`surf_plot*.png`, `synthetic_*.png`).

**Scripts at a glance:**

| Script | Purpose | Run |
|--------|---------|-----|
| `visualize.py` | Demo: random sample, fit SURF, plot density + histogram | `python visualize.py` |
| `visualize_salary.py` | Fit on `data/salary.dat`, save plot | `python visualize_salary.py` |
| `visualize_sensor.py` | Fit on `data/sensor1.dat`, save plot | `python visualize_sensor.py` |
| `synthetic_experiments.py` | Beta/gamma/Gaussian mixtures: L1 error + plots | `python synthetic_experiments.py` |
| `run_and_show_polynomials.py` | Small random sample; print piecewise polynomials | `python run_and_show_polynomials.py` |

---

## Requirements

- **Python 3.8+**, **NumPy**, **matplotlib**; **scipy** for `synthetic_experiments.py`. See `requirements.txt` or `pyproject.toml`.
- **MATLAB** (optional): only if you want to run the legacy implementation in `legacy/`.

---

## Project structure

```
SURF/
├── README.md
├── pyproject.toml
├── requirements.txt
├── surf/                   # Python package: surf(), regpoly()
├── visualize.py            # Plot on random sample
├── visualize_salary.py     # Plot on data/salary.dat
├── visualize_sensor.py     # Plot on data/sensor1.dat
├── synthetic_experiments.py # Beta/gamma/Gaussian mixtures, L1, plots
├── run_and_show_polynomials.py
├── legacy/                 # Original MATLAB (NeurIPS 2020)
│   ├── surf.m, merge.m, coeffint.m, ...
│   ├── salary_experiments.m, sensor_experiments.m, ...
│   └── README.md
└── data/                   # Datasets (salary.dat, sensor1.dat, etc.)
```

---

## Python API

- **`surf(samples, alpha=..., degree=...)`** — Returns `(boundaries, piece_coeffs)`. `samples`: sorted 1D array of size `2^k - 1` (e.g. 31, 255, 1023) with values strictly in (0, 1). `alpha`: tuning (larger → fewer pieces). `degree`: 1–4.
- **`regpoly(x, coeffs)`** — Evaluate polynomial at `x`; `coeffs` is constant term first, e.g. `[a0, a1, a2]` for a0 + a1*x + a2*x².

**Data format:** For your own data, sort and scale to (0, 1), then take exactly `2^k - 1` points (e.g. subsample or pad). See `visualize_salary.py` / `visualize_sensor.py` for load-and-scale examples.

**Tuning:** Larger α → fewer pieces (more merging). Higher degree → more detail per piece. The criterion uses only α (no extra constant).

**Visualization:** `visualize.run_and_plot(samples, alpha, degree, out_path, title=..., xlabel=...)` runs SURF and saves a density + histogram plot. Use it from your own scripts after preprocessing your data to a sorted array of size `2^k - 1` in (0, 1).

---

## Legacy (MATLAB)

The original implementation and experiment scripts are in **`legacy/`**. From the project root:

```matlab
addpath('legacy');
[I, koi] = surf(sort(rand(1, 2^8 - 1)), 0.25, 1);
```

See **`legacy/README.md`** for experiment scripts and data paths.

---

## Citation

If you use this code in your work, please cite:

```bibtex
@inproceedings{surf2020,
  title     = {SURF: Learning Univariate Distributions},
  booktitle = {Advances in Neural Information Processing Systems (NeurIPS)},
  year      = {2020}
}
```

---

## License

See the repository for license information.
