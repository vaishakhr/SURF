# SURF

**SURF** is an algorithm for learning univariate distributions from samples. This repository contains a **MATLAB** implementation and experiment scripts from our **NeurIPS 2020** paper, plus an **efficient Python implementation** in `python/`.

---

## Overview

SURF takes a set of samples and returns a piecewise polynomial density estimator over the unit interval. The algorithm works by recursively merging intervals according to a criterion controlled by a tuning parameter α.

---

## Requirements

- **MATLAB** (R2016a or later recommended) for the original implementation
- **Python 3.8+** and **NumPy** for the optimized implementation in `python/` (see `python/README.md`)

---

## Quick Start

1. **Clone or download** this repository.
2. Add the `src` folder to your MATLAB path:
   ```matlab
   addpath('src');
   ```
3. Call the main function with sorted samples in (0, 1), a tuning parameter, and polynomial degree:
   ```matlab
   samp = sort(rand(1, 2^8 - 1));  % n = 2^8 - 1 samples
   [I, koi] = surf(samp, 0.25, 1);
   ```

---

## API

### `[I, koi] = surf(samp, alp, deg)`

| Argument | Description |
|----------|-------------|
| `samp`   | Sorted samples **strictly in (0, 1)**. The number of samples must be **2^k − 1** for some integer k (e.g. 255, 511, 1023). |
| `alp`    | Tuning parameter α. A good default is **0.25**. |
| `deg`    | Degree of the piecewise polynomial. Supported values: **1 ≤ deg ≤ 4**. |

| Output | Description |
|--------|-------------|
| `I`     | Interval boundaries (partition of [0, 1]). |
| `koi`   | Polynomial coefficients for each piece (one row per interval). |

---

## Experiment Scripts

Scripts reproducing experiments from the NeurIPS 2020 paper are in `src/`:

| File | Description |
|------|-------------|
| `synthetic_experiments.m` | Experiments on synthetic data from known distributions (e.g. beta mixtures). |
| `salary_experiments.m`    | Experiments on salary data. Requires `data/salary.dat`. |
| `sensor_experiments.m`    | Sensor data experiments. |
| `cover_experiments.m`     | Cover-type experiments. |

To run from the project root:

```matlab
addpath('src');
run('src/synthetic_experiments.m');   % or salary_experiments.m, etc.
```

**Note:** `salary_experiments.m` expects the dataset at `data/salary.dat`. For large n (e.g. 2^13) and higher degrees, runtimes can be long (see below).

---

## Project Structure

```
SURF/
├── README.md
├── src/                    % MATLAB implementation
│   ├── surf.m              % Main SURF function
│   ├── merge.m             % Interval merging
│   ├── coeffint.m          % Coefficient computation
│   ├── regpoly.m           % Polynomial evaluation
│   ├── matcoeff.m          % Matrix coefficients
│   ├── maxcrit.m           % Criterion maximization
│   ├── differ.m            % Difference utilities
│   ├── synthetic_experiments.m
│   ├── salary_experiments.m
│   ├── sensor_experiments.m
│   └── cover_experiments.m
├── python/                 % Efficient Python implementation (see python/README.md)
│   ├── surf/               % Package: surf(), regpoly()
│   ├── visualize.py        % Plot SURF on random sample
│   ├── visualize_salary.py % Plot SURF on data/salary.dat (like salary_experiments.m)
│   ├── requirements.txt
│   └── README.md
└── data/                   % Place salary.dat here for salary experiments
```

---

## Performance Notes

- The implementation is **recursive** and **not optimized for speed**.
- Running time depends on **number of samples (n)**, **degree (deg)**, and your machine. For large n or deg &gt; 2, runs can take several minutes or more.
- For real-data experiments (e.g. salary), the code is tested up to about **n = 2^13** with **deg ≤ 2**.

---

## Citation

If you use this code in your work, please cite our NeurIPS 2020 paper:

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
