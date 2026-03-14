"""
Run SURF and plot the fitted density + sample histogram. Reusable for any dataset.

Other scripts (e.g. visualize_salary.py) preprocess their data, then call
run_and_plot(samples, alpha, degree, out_path, ...) to fit and save the figure.

Requires matplotlib. Run as script for a small random-sample demo.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from surf import surf, regpoly


def density_at(x, boundaries, piece_coeffs):
    """Evaluate the SURF density at x. x can be a scalar or array."""
    x = np.asarray(x, dtype=np.float64)
    idx = np.searchsorted(boundaries, x, side="right") - 1
    idx = np.clip(idx, 0, piece_coeffs.shape[0] - 1)
    out = np.zeros_like(x)
    for i in np.unique(idx):
        mask = idx == i
        out[mask] = regpoly(x[mask], piece_coeffs[i, :])
    return out


def run_and_plot(
    samples,
    alpha,
    degree,
    out_path,
    *,
    title=None,
    xlabel="x",
):
    """
    Run SURF on samples, plot density + histogram, save to out_path.

    samples: sorted 1D array of size 2^k - 1 in (0, 1)
    alpha, degree: SURF parameters
    out_path: path for the saved PNG
    title: plot title. If it contains "{num_pieces}", that is filled in.
           If None, uses "SURF fit — {num_pieces} pieces".
    xlabel: x-axis label

    Returns (boundaries, piece_coeffs, num_pieces).
    """
    import matplotlib.pyplot as plt

    boundaries, piece_coeffs = surf(samples, alpha=alpha, degree=degree)
    num_pieces = len(boundaries) - 1

    x_plot = np.linspace(0.001, 0.999, 500)
    f_plot = density_at(x_plot, boundaries, piece_coeffs)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.fill_between(x_plot, f_plot, alpha=0.3, label="SURF density")
    ax.plot(x_plot, f_plot, color="C0", lw=1.5, label="p(x)")
    ax.hist(
        samples,
        bins=min(60, max(20, len(samples) // 200)),
        density=True,
        alpha=0.5,
        color="gray",
        label="samples (hist)",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("density")
    ax.legend(loc="upper right")
    if title is None:
        title = f"SURF fit — {num_pieces} pieces"
    else:
        title = title.format(num_pieces=num_pieces)
    ax.set_title(title)
    ax.text(
        0.02, 0.98,
        f"α = {alpha}\nd = {degree}\n# pieces = {num_pieces}",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.9),
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    return boundaries, piece_coeffs, num_pieces


def main():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError:
        print("matplotlib is required. Install with: pip install matplotlib")
        return

    np.random.seed(42)
    num_samples = 2**6 - 1
    samples = np.sort(np.random.uniform(0.02, 0.98, num_samples))

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(script_dir, "surf_plot.png")
    boundaries, piece_coeffs, num_pieces = run_and_plot(
        samples, alpha=1.0, degree=1, out_path=out_path,
        title="SURF fit — {num_pieces} pieces",
    )
    print(f"SURF fit: {num_pieces} pieces")
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
