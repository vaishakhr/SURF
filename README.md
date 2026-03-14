# SURF

**SURF** is an algorithm for learning univariate distributions from samples. This repository is built around an **efficient Python implementation** at the repo root; the original **MATLAB** code from our NeurIPS 2020 paper lives in **`legacy/`** for reference.

---

## Overview

SURF takes a set of samples and returns a piecewise polynomial density estimator over the unit interval. The algorithm recursively merges intervals according to a criterion controlled by a tuning parameter α.

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

Plots are saved in the repo root (`surf_plot.png`, `surf_plot_salary_*.png`, etc.).

---

## Requirements

- **Python 3.8+**, **NumPy**, **matplotlib**. See `requirements.txt` or `pyproject.toml`.
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
├── run_and_show_polynomials.py
├── legacy/                 # Original MATLAB (NeurIPS 2020)
│   ├── surf.m, merge.m, coeffint.m, ...
│   ├── salary_experiments.m, sensor_experiments.m, ...
│   └── README.md
└── data/                   # Datasets (salary.dat, sensor1.dat, etc.)
```

---

## Python API

- **`surf(samples, alpha, degree)`** — Returns `(boundaries, piece_coeffs)`. `samples`: sorted 1D array of size `2^k - 1` in (0, 1). `alpha`: tuning (larger → fewer pieces). `degree`: 1–4.
- **`regpoly(x, coeffs)`** — Evaluate polynomial (constant term first).

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
