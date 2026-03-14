"""
SURF on salary data (same setup as MATLAB salary_experiments.m).

Only preprocessing: loads data/salary.dat, scales to (0, 1), then calls
visualize.run_and_plot() to fit and save the figure. Edit alp and deg below.

Run: cd python && PYTHONPATH=. python visualize_salary.py
"""

import os
import numpy as np
from visualize import run_and_plot


def load_salary_chunk(data_path, n=2**13, offset=0):
    """
    Load salary data and return SURF-ready samples for one chunk.

    Same preprocessing as MATLAB salary_experiments.m: take chunk of size n+1,
    sort, drop min/max, scale to (0, 1) → 2^k - 1 samples. n must be 2^k.

    Returns
    -------
    samp : 1D array of length n - 1, sorted in (0, 1)
    """
    sal = np.loadtxt(data_path, delimiter=",").ravel()
    sen_train = sal[offset : offset + n + 1]
    sen_train = np.sort(sen_train)
    loc1 = sen_train[0]
    loc2 = sen_train[-1]
    scale = loc2 - loc1
    if scale <= 0:
        raise ValueError("Scale is zero (constant chunk).")
    samp = (sen_train[1:-1] - loc1) / scale
    samp = np.clip(samp, 1e-6, 1 - 1e-6)
    return samp


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_path = os.path.join(project_root, "data", "salary.dat")

    if not os.path.isfile(data_path):
        print(f"Data file not found: {data_path}")
        print("Place salary.dat in the project's data/ folder (see README).")
        return

    alp = 1.0
    deg = 2

    samp = load_salary_chunk(data_path)
    print(f"Preprocessed {len(samp)} samples from salary.dat (one chunk, scaled to (0,1))")

    out_path = os.path.join(script_dir, f"surf_plot_salary_{alp}.png")
    I, koi, num_pieces = run_and_plot(
        samp,
        alp=alp,
        deg=deg,
        out_path=out_path,
        title="SURF on salary data — {num_pieces} pieces",
        xlabel="x (salary scaled to [0,1])",
    )
    print(f"SURF fit: {num_pieces} pieces (deg={deg}, alpha={alp})")
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
