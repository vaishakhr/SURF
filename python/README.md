# SURF — Python implementation

Efficient Python port of the SURF algorithm (learning univariate distributions), with the same API as the MATLAB version and several performance optimizations.

## Optimizations

- **Binary search** for sample counts (`numpy.searchsorted`) instead of scanning the full array
- **Exact integral** of \|polynomial\| in the criterion (roots + antiderivative) instead of numerical integration
- **Reuse of precomputed coefficients** for single-interval blocks in the recursive criterion
- **Restricted sample slices** passed into `coeffint` so all counting is on the relevant slice
- **Vectorized** coefficient matrix and polynomial evaluation

## Install

From the project root:

```bash
pip install -e ./python
```

Or use the environment of your choice and ensure `numpy` is installed (`python/requirements.txt`).

## Usage

```python
import numpy as np
from surf import surf, regpoly

# n = 2^k - 1 samples, sorted in (0, 1)
n = 2**8 - 1
samp = np.sort(np.random.uniform(0.01, 0.99, n))

I, koi = surf(samp, alp=0.25, deg=1)

# Evaluate density at many points (see visualize.py for a full density_at helper)
x = np.linspace(0.01, 0.99, 200)
idx = np.clip(np.searchsorted(I, x, side="right") - 1, 0, koi.shape[0] - 1)
f = np.array([regpoly(xi, koi[idx[i], :]) for i, xi in enumerate(x)])
```

## Tuning

- **α (alp):** Merge penalty scale. Larger α → fewer pieces (more merging); smaller α → more pieces. The criterion uses only α (no extra constant in the formula).
- **deg:** Piecewise polynomial degree (1–4). Higher degree can fit more detail per piece.

## API

- **`surf(samp, alp, deg)`** — `samp`: sorted 1D array of size `2^k - 1` in (0, 1); `alp`: tuning parameter α (larger → fewer pieces); `deg`: 1–4. Returns `(I, koi)` with interval boundaries and piecewise polynomial coefficients.
- **`regpoly(x, a)`** — Evaluate polynomial with coefficients `a` (constant term first) at `x`.

## Visualization

- **`visualize.py`** — Defines `run_and_plot(samp, alp, deg, out_path, ...)` to run SURF and save a density + histogram plot. Run as a script for a random-sample demo (saves `surf_plot.png`).
- **`visualize_salary.py`** — Preprocesses `data/salary.dat` (same as MATLAB `salary_experiments.m`), then calls `visualize.run_and_plot()`. Saves `surf_plot_salary_{alp}.png`. Edit `alp` and `deg` in the script to try different fits.

For other datasets: preprocess to a sorted array `samp` of size `2^k - 1` in (0, 1), then call `run_and_plot(samp, alp, deg, out_path, title=..., xlabel=...)`.

```bash
cd python && PYTHONPATH=. python visualize.py
# or, with salary data in data/salary.dat:
cd python && PYTHONPATH=. python visualize_salary.py
```

Requires `matplotlib` (in `requirements.txt`). To verify the install, run `visualize.py` or `visualize_salary.py` (with `data/salary.dat` in place).
