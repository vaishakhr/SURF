"""
Replicate synthetic data experiments from MATLAB (src/synthetic_experiments.m).

Experiments on data from known distributions: beta mixture, gamma mixture,
gaussian mixture. Run SURF, compute L1 difference (discrete), and plot
true density vs SURF estimate.

Parameters: alpha (tuning), degree (polynomial degree), n (power of 2).
Defaults used here: alpha = 1.0, degree = 2. Code supports degree <= 4.
"""

import os
import numpy as np
from scipy.special import beta as beta_fn, gamma as gamma_fn

# Optional matplotlib for plotting
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _has_plt = True
except ImportError:
    _has_plt = False

from surf import surf, regpoly


def surf_estim_at(x, boundaries, piece_coeffs):
    """Evaluate SURF estimator at x (scalar or array). No normalization."""
    x = np.asarray(x, dtype=np.float64)
    idx = np.searchsorted(boundaries, x, side="right") - 1
    idx = np.clip(idx, 0, piece_coeffs.shape[0] - 1)
    out = np.zeros_like(x)
    for i in np.unique(idx):
        mask = idx == i
        out[mask] = regpoly(x[mask], piece_coeffs[i, :])
    return out


def run_beta_mixture(n_trials=10, seed=None):
    """Beta mixture: prob * Beta(beta_alp, beta_bet) + (1-prob) * Beta(beta_alp_1, beta_bet_1)."""
    degree = 2
    alpha = 1.0
    beta_alp, beta_bet = 0.8, 4.0
    beta_alp_1, beta_bet_1 = 2.0, 2.0
    prob = 0.4
    n = 2**9  # n-1 samples
    if seed is not None:
        np.random.seed(seed)

    pres = 16 * n + 1
    Y = np.linspace(0, 1, pres + 1)
    skip = 32  # avoid blow-up at 0

    def true_density(x):
        # Clip to avoid 0^neg or 1^neg at boundaries (beta PDF can blow up at 0/1)
        x = np.clip(np.asarray(x, dtype=np.float64), 1e-12, 1.0 - 1e-12)
        t = prob * (x ** (beta_alp - 1)) * ((1 - x) ** (beta_bet - 1)) / beta_fn(beta_alp, beta_bet)
        t += (1 - prob) * (x ** (beta_alp_1 - 1)) * ((1 - x) ** (beta_bet_1 - 1)) / beta_fn(beta_alp_1, beta_bet_1)
        return t

    l1_total = 0.0
    last_boundaries = last_coeffs = None
    last_Y = last_true = last_estim = None

    for _ in range(n_trials):
        f = (np.random.rand(n - 1) <= prob).astype(np.float64)
        samp = f * np.random.beta(beta_alp, beta_bet, n - 1) + (1 - f) * np.random.beta(beta_alp_1, beta_bet_1, n - 1)
        samp = np.sort(samp)

        boundaries, piece_coeffs = surf(samp, alpha=alpha, degree=degree)
        last_boundaries, last_coeffs = boundaries, piece_coeffs

        density_disc = np.array([true_density(y) for y in Y])
        estim_disc = surf_estim_at(Y, boundaries, piece_coeffs)
        last_Y, last_true, last_estim = Y, density_disc, estim_disc

        l1_total += np.sum(np.abs(density_disc[skip:] - estim_disc[skip:])) / pres

    l1_mean = l1_total / n_trials
    return {
        "l1_mean": l1_mean,
        "Y": last_Y,
        "true_density": last_true,
        "estim_density": last_estim,
        "boundaries": last_boundaries,
        "piece_coeffs": last_coeffs,
        "name": "Beta mixture",
    }


def run_gamma_mixture(n_trials=10, seed=None):
    """Gamma mixture; tail trimmed and rescaled to [0,1]."""
    degree = 2
    alpha = 1.0
    prob = 0.2
    gam_a, gam_b = 4.0, 0.04
    gam_a_1, gam_b_1 = 8.0, 0.06
    n = 2**10
    n_ex = int(np.ceil(n**0.1))
    if seed is not None:
        np.random.seed(seed)

    pres = 16 * n + 1
    Y = np.linspace(0, 1, pres + 1)
    skip = 32

    def true_density(x):
        # x is on original scale (before rescale). Return density on original scale.
        t = prob * (x ** (gam_a - 1)) * np.exp(-x / gam_b) / (gamma_fn(gam_a) * gam_b**gam_a)
        t += (1 - prob) * (x ** (gam_a_1 - 1)) * np.exp(-x / gam_b_1) / (gamma_fn(gam_a_1) * gam_b_1**gam_a_1)
        return t

    l1_total = 0.0
    last_Y = last_true = last_estim = None
    last_boundaries = last_coeffs = None

    for _ in range(n_trials):
        m = n + n_ex
        f = (np.random.rand(m) <= prob).astype(np.float64)
        samp = f * np.random.gamma(gam_a, scale=gam_b, size=m) + (1 - f) * np.random.gamma(gam_a_1, scale=gam_b_1, size=m)
        samp = np.sort(samp)
        scale = samp[n - 1]  # samp[n] in 1-based; samp[n-1] is n-th smallest
        samp = samp[: n - 1] / scale  # samp(1:n-1) in MATLAB

        boundaries, piece_coeffs = surf(samp, alpha=alpha, degree=degree)
        last_boundaries, last_coeffs = boundaries, piece_coeffs

        # True density on rescaled grid: Y*scale is original scale
        density_disc = scale * np.array([true_density(y * scale) for y in Y])
        estim_disc = surf_estim_at(Y, boundaries, piece_coeffs)
        last_Y, last_true, last_estim = Y, density_disc, estim_disc

        l1_total += np.sum(np.abs(density_disc[skip:] - estim_disc[skip:])) / pres + n_ex / n

    l1_mean = l1_total / n_trials
    return {
        "l1_mean": l1_mean,
        "Y": last_Y,
        "true_density": last_true,
        "estim_density": last_estim,
        "boundaries": last_boundaries,
        "piece_coeffs": last_coeffs,
        "name": "Gamma mixture",
    }


def run_gaussian_mixture(n_trials=5, seed=None):
    """Gaussian mixture; outliers trimmed and segment rescaled to [0,1]."""
    degree = 2
    alpha = 1.0
    gauss_m, gauss_sd = 0.4, 0.1
    gauss_m_1, gauss_sd_1 = 0.6, 0.2
    prob = 0.3
    # MATLAB uses n=2^15; smaller n here so it runs in reasonable time
    n = 2**11
    n_ex = int(np.ceil(n**0.1))
    if seed is not None:
        np.random.seed(seed)

    pres = 16 * n + 1
    Y = np.linspace(0, 1, pres + 1)
    skip = 32

    def true_density(x):
        t = prob * np.exp(-((x - gauss_m) ** 2) / (2 * gauss_sd**2)) / np.sqrt(2 * np.pi * gauss_sd**2)
        t += (1 - prob) * np.exp(-((x - gauss_m_1) ** 2) / (2 * gauss_sd_1**2)) / np.sqrt(2 * np.pi * gauss_sd_1**2)
        return t

    l1_total = 0.0
    last_Y = last_true = last_estim = None
    last_boundaries = last_coeffs = None

    for _ in range(n_trials):
        m = n + 2 * n_ex
        f = (np.random.rand(m) <= prob).astype(np.float64)
        samp = f * np.random.normal(gauss_m, gauss_sd, m) + (1 - f) * np.random.normal(gauss_m_1, gauss_sd_1, m)
        samp = np.sort(samp)
        loc1 = samp[n_ex - 1]   # 1-based samp(n_ex)
        loc2 = samp[n + n_ex - 1]
        scale = loc2 - loc1
        samp = (samp[n_ex : n + n_ex - 1] - loc1) / scale  # n_ex+1 : n+n_ex-1 in 1-based -> n_ex to n+n_ex-2 incl

        boundaries, piece_coeffs = surf(samp, alpha=alpha, degree=degree)
        last_boundaries, last_coeffs = boundaries, piece_coeffs

        # Y in [0,1] maps to original scale Y*scale + loc1
        density_disc = scale * np.array([true_density(y * scale + loc1) for y in Y])
        estim_disc = surf_estim_at(Y, boundaries, piece_coeffs)
        last_Y, last_true, last_estim = Y, density_disc, estim_disc

        l1_total += np.sum(np.abs(density_disc[skip:] - estim_disc[skip:])) / pres + 2 * n_ex / n

    l1_mean = l1_total / n_trials
    return {
        "l1_mean": l1_mean,
        "Y": last_Y,
        "true_density": last_true,
        "estim_density": last_estim,
        "boundaries": last_boundaries,
        "piece_coeffs": last_coeffs,
        "name": "Gaussian mixture",
    }


def plot_result(result, out_path):
    """Plot true density vs SURF estimate and save to out_path."""
    if not _has_plt:
        print("matplotlib not available, skip plot")
        return
    plt.figure(figsize=(8, 4))
    plt.plot(result["Y"], result["true_density"], label="True density")
    plt.plot(result["Y"], result["estim_density"], label="SURF estimate")
    plt.xlim(0, 1)
    # Cap y-axis at 99th percentile so beta (and similar) spikes at boundaries don't flatten the plot
    combined = np.concatenate([result["true_density"], result["estim_density"]])
    finite = combined[np.isfinite(combined)]
    y_max = np.percentile(finite, 99) * 1.15 if len(finite) > 0 else 1.0
    plt.ylim(0, max(y_max, 0.1))
    plt.xlabel("x")
    plt.ylabel("density")
    plt.legend()
    plt.title(f"{result['name']} — L1 (mean) = {result['l1_mean']:.6f}")
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()
    print(f"Saved {out_path}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    print("Beta mixture (n=2^9, 10 trials)...")
    r_beta = run_beta_mixture(n_trials=10, seed=42)
    print(f"  L1 mean = {r_beta['l1_mean']:.6f}, # pieces = {len(r_beta['boundaries'])-1}")
    plot_result(r_beta, os.path.join(script_dir, "synthetic_beta.png"))

    print("Gamma mixture (n=2^10, 10 trials)...")
    r_gamma = run_gamma_mixture(n_trials=10, seed=43)
    print(f"  L1 mean = {r_gamma['l1_mean']:.6f}, # pieces = {len(r_gamma['boundaries'])-1}")
    plot_result(r_gamma, os.path.join(script_dir, "synthetic_gamma.png"))

    print("Gaussian mixture (n=2^11, 5 trials)...")
    r_gauss = run_gaussian_mixture(n_trials=5, seed=44)
    print(f"  L1 mean = {r_gauss['l1_mean']:.6f}, # pieces = {len(r_gauss['boundaries'])-1}")
    plot_result(r_gauss, os.path.join(script_dir, "synthetic_gaussian.png"))

    print("Done.")


if __name__ == "__main__":
    main()
