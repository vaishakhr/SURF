"""
SURF on salary data (same setup as MATLAB salary_experiments.m).

Only preprocessing: loads data/salary.dat, scales to (0, 1), then calls
visualize.run_and_plot() to fit and save the figure. Edit alpha and degree below.

Run: cd python && PYTHONPATH=. python visualize_salary.py
"""

import os
import numpy as np
from visualize import run_and_plot


def load_salary_chunk(data_path, chunk_size=2**13, offset=0):
    """
    Load salary data and return SURF-ready samples for one chunk.

    Same preprocessing as MATLAB salary_experiments.m: take chunk of size chunk_size+1,
    sort, drop min/max, scale to (0, 1) → chunk_size - 1 samples. chunk_size must be 2^k.

    Returns
    -------
    samples : 1D array of length chunk_size - 1, sorted in (0, 1)
    """
    raw = np.loadtxt(data_path, delimiter=",").ravel()
    chunk = raw[offset : offset + chunk_size + 1]
    chunk = np.sort(chunk)
    min_val = chunk[0]
    max_val = chunk[-1]
    scale = max_val - min_val
    if scale <= 0:
        raise ValueError("Scale is zero (constant chunk).")
    samples = (chunk[1:-1] - min_val) / scale
    samples = np.clip(samples, 1e-6, 1 - 1e-6)
    return samples


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_path = os.path.join(project_root, "data", "salary.dat")

    if not os.path.isfile(data_path):
        print(f"Data file not found: {data_path}")
        print("Place salary.dat in the project's data/ folder (see README).")
        return

    alpha = 1.0
    degree = 2

    samples = load_salary_chunk(data_path)
    print(f"Preprocessed {len(samples)} samples from salary.dat (one chunk, scaled to (0,1))")

    out_path = os.path.join(script_dir, f"surf_plot_salary_{alpha}.png")
    boundaries, piece_coeffs, num_pieces = run_and_plot(
        samples,
        alpha=alpha,
        degree=degree,
        out_path=out_path,
        title="SURF on salary data — {num_pieces} pieces",
        xlabel="x (salary scaled to [0,1])",
    )
    print(f"SURF fit: {num_pieces} pieces (degree={degree}, alpha={alpha})")
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
