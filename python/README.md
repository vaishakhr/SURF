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

# num_samples = 2^k - 1, sorted in (0, 1)
num_samples = 2**8 - 1
samples = np.sort(np.random.uniform(0.01, 0.99, num_samples))

boundaries, piece_coeffs = surf(samples, alpha=0.25, degree=1)

# Evaluate density at many points (see visualize.py for density_at helper)
x = np.linspace(0.01, 0.99, 200)
idx = np.clip(np.searchsorted(boundaries, x, side="right") - 1, 0, piece_coeffs.shape[0] - 1)
f = np.array([regpoly(xi, piece_coeffs[idx[i], :]) for i, xi in enumerate(x)])
```

## Tuning

- **alpha:** Merge penalty scale. Larger α → fewer pieces (more merging); smaller α → more pieces. The criterion uses only α (no extra constant in the formula).
- **degree:** Piecewise polynomial degree (1–4). Higher degree can fit more detail per piece.

## API

- **`surf(samples, alpha, degree)`** — `samples`: sorted 1D array of size `2^k - 1` in (0, 1); `alpha`: tuning parameter (larger → fewer pieces); `degree`: 1–4. Returns `(boundaries, piece_coeffs)` with interval boundaries and polynomial coefficients per piece.
- **`regpoly(x, a)`** — Evaluate polynomial with coefficients `a` (constant term first) at `x`.

## Visualization

- **`visualize.py`** — Defines `run_and_plot(samples, alpha, degree, out_path, ...)` to run SURF and save a density + histogram plot. Run as a script for a random-sample demo (saves `surf_plot.png`).
- **`visualize_salary.py`** — Preprocesses `data/salary.dat` (same as MATLAB `salary_experiments.m`), then calls `visualize.run_and_plot()`. Saves `surf_plot_salary_{alpha}.png`. Edit `alpha` and `degree` in the script to try different fits.
- **`visualize_sensor.py`** — Preprocesses `data/sensor1.dat` (same as MATLAB `sensor_experiments.m`, chunk size 2^14), then calls `visualize.run_and_plot()`. Saves `surf_plot_sensor_{alpha}.png`.

For other datasets: preprocess to a sorted array `samples` of size `2^k - 1` in (0, 1), then call `run_and_plot(samples, alpha, degree, out_path, title=..., xlabel=...)`.

```bash
cd python && PYTHONPATH=. python visualize.py
# or, with data in data/:
cd python && PYTHONPATH=. python visualize_salary.py   # data/salary.dat
cd python && PYTHONPATH=. python visualize_sensor.py   # data/sensor1.dat
```

Requires `matplotlib` (in `requirements.txt`). To verify the install, run `visualize.py` or `visualize_salary.py` (with `data/salary.dat` in place).
