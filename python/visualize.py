"""
Run SURF and plot the fitted density + sample histogram. Reusable for any dataset.

Other scripts (e.g. visualize_salary.py) preprocess their data, then call
run_and_plot(samp, alp, deg, out_path, ...) to fit and save the figure.

Requires matplotlib. Run as script for a small random-sample demo.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
from surf import surf, regpoly


def density_at(x, I, koi):
    """Evaluate the SURF density at x. x can be a scalar or array."""
    x = np.asarray(x, dtype=np.float64)
    idx = np.searchsorted(I, x, side="right") - 1
    idx = np.clip(idx, 0, koi.shape[0] - 1)
    out = np.zeros_like(x)
    for i in np.unique(idx):
        mask = idx == i
        out[mask] = regpoly(x[mask], koi[i, :])
    return out


def run_and_plot(
    samp,
    alp,
    deg,
    out_path,
    *,
    title=None,
    xlabel="x",
):
    """
    Run SURF on samp, plot density + histogram, save to out_path.

    samp: sorted 1D array of size 2^k - 1 in (0, 1)
    alp, deg: SURF parameters
    out_path: path for the saved PNG
    title: plot title. If it contains "{num_pieces}", that is filled in.
           If None, uses "SURF fit — {num_pieces} pieces".
    xlabel: x-axis label

    Returns (I, koi, num_pieces).
    """
    import matplotlib.pyplot as plt

    I, koi = surf(samp, alp=alp, deg=deg)
    num_pieces = len(I) - 1

    x_plot = np.linspace(0.001, 0.999, 500)
    f_plot = density_at(x_plot, I, koi)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.fill_between(x_plot, f_plot, alpha=0.3, label="SURF density")
    ax.plot(x_plot, f_plot, color="C0", lw=1.5, label="p(x)")
    ax.hist(
        samp,
        bins=min(60, max(20, len(samp) // 200)),
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
        f"α = {alp}\nd = {deg}\n# pieces = {num_pieces}",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.9),
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    return I, koi, num_pieces


def main():
    try:
        import matplotlib.pyplot as plt  # noqa: F401
    except ImportError:
        print("matplotlib is required. Install with: pip install matplotlib")
        return

    np.random.seed(42)
    n = 2**6 - 1
    samp = np.sort(np.random.uniform(0.02, 0.98, n))

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(script_dir, "surf_plot.png")
    I, koi, num_pieces = run_and_plot(
        samp, alp=1.0, deg=1, out_path=out_path,
        title="SURF fit — {num_pieces} pieces",
    )
    print(f"SURF fit: {num_pieces} pieces")
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
