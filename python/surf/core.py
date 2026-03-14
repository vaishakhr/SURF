"""
Efficient Python implementation of SURF (learning univariate distributions).

Optimizations applied:
- Binary search (numpy.searchsorted) for sample counts instead of full scans
- Exact integral of |polynomial| in differ() instead of numerical integration
- Reuse of precomputed single-interval coefficients in maxcrit
- Restrict sample slice before coeffint where possible
- Vectorized coefficient matrix and polynomial evaluation
"""

import numpy as np
from numpy.typing import NDArray

# Node matrix for polynomial degree (rows 0..4 = deg 0..4), first deg+2 cols used
OPND1 = np.array([
    [0, 1, 0, 0, 0, 0],
    [0, 0.5, 1, 0, 0, 0],
    [0, 0.2599, 0.7401, 1, 0, 0],
    [0, 0.1548, 0.5, 0.8452, 1, 0],
    [0, 0.1015, 0.348, 0.652, 0.8985, 1],
], dtype=np.float64)


def regpoly(x: NDArray[np.float64] | float, a: NDArray[np.float64]) -> NDArray[np.float64] | float:
    """Evaluate polynomial with coefficients a at x. a[0] + a[1]*x + a[2]*x**2 + ..."""
    x = np.asarray(x, dtype=np.float64)
    p = np.zeros_like(x)
    for i, c in enumerate(a):
        p = p + c * (x ** i)
    return p if p.shape else float(p.flat[0])


def _matcoeff_row(n1: float, n2: float, n: int, degree: int) -> NDArray[np.float64]:
    """One row of the coefficient matrix: n * (n2^i - n1^i) / i for i in 1..degree+1."""
    i = np.arange(1, degree + 2, dtype=np.float64)
    return n * (np.power(n2, i) - np.power(n1, i)) / i


def _count_in_interval(samples: NDArray[np.float64], a: float, b: float) -> int:
    """Count samples strictly in (a, b) using binary search. samples must be sorted."""
    if samples.size == 0 or a >= b:
        return 0
    left = np.searchsorted(samples, a, side="right")
    right = np.searchsorted(samples, b, side="left")
    return max(0, right - left)


def _coeffint(
    samples: NDArray[np.float64],
    left: float,
    right: float,
    n: int,
    degree: int,
) -> NDArray[np.float64]:
    """
    Coefficients for polynomial fit on interval [left, right].
    Caller may pass samples restricted to [left, right] for efficiency; must be sorted.
    """
    interval_len = right - left
    if interval_len <= 0 or interval_len < 1e-12:
        return np.zeros(degree + 1, dtype=np.float64)
    row = OPND1[degree]
    nodes = row[: degree + 2] * interval_len + left

    count_in_interval = _count_in_interval(samples, left, right)
    normalizer = interval_len * (count_in_interval + 1) / (interval_len * count_in_interval + 1)

    emp_vec = np.zeros(degree + 1, dtype=np.float64)
    for i in range(degree + 1):
        bin_left, bin_right = nodes[i], nodes[i + 1]
        cnt = _count_in_interval(samples, bin_left, bin_right)
        delta_node = (row[i + 1] - row[i]) / interval_len
        emp_vec[i] = (delta_node + cnt) * normalizer * (degree + 1)

    mat = np.zeros((degree + 1, degree + 1), dtype=np.float64)
    scale = float(degree + 1)
    for j in range(degree + 1):
        mat[j, :] = scale * _matcoeff_row(nodes[j], nodes[j + 1], n, degree)

    try:
        coeff = np.linalg.solve(mat, emp_vec)
    except np.linalg.LinAlgError:
        coeff = np.zeros(degree + 1, dtype=np.float64)
    return coeff


def _integral_abs_poly(a: float, b: float, coeffs: NDArray[np.float64]) -> float:
    """
    Compute integral of |p(x)| over [a, b] where p = coeffs[0] + coeffs[1]*x + ...
    Splits at roots of p in (a,b), then integrates p or -p exactly on each subinterval.
    """
    coeffs = np.asarray(coeffs, dtype=np.float64).ravel()
    if coeffs.size == 0 or b <= a:
        return 0.0
    # np.polyval expects descending order (high power first)
    p_desc = coeffs[::-1]
    roots = np.roots(p_desc)
    roots = np.real(roots[np.isreal(roots)])
    in_ab = roots[(roots > a) & (roots < b)]
    split = np.unique(np.asarray([a, *in_ab, b], dtype=np.float64))
    total = 0.0
    for i in range(len(split) - 1):
        x0, x1 = split[i], split[i + 1]
        mid = (x0 + x1) / 2
        sign = 1.0 if np.polyval(p_desc, mid) >= 0 else -1.0
        # Antiderivative of sign*p: P(x) = sign * sum_k c_k x^{k+1}/(k+1)
        prim_asc = np.zeros(coeffs.size + 1)
        for k, c in enumerate(coeffs):
            prim_asc[k + 1] = sign * c / (k + 1)
        P_x1 = np.polyval(prim_asc[::-1], x1)
        P_x0 = np.polyval(prim_asc[::-1], x0)
        total += P_x1 - P_x0
    return float(total)


def _differ(left: float, right: float, coeff_a: NDArray[np.float64], coeff_b: NDArray[np.float64]) -> float:
    """L1 difference between polynomials coeff_a and coeff_b over [left, right] (exact integral)."""
    diff = np.asarray(coeff_a, dtype=np.float64).ravel() - np.asarray(coeff_b, dtype=np.float64).ravel()
    return _integral_abs_poly(left, right, diff)


def _samples_in_interval(samples: NDArray[np.float64], a: float, b: float) -> tuple[NDArray[np.float64], int, int]:
    """Return samples restricted to [a, b] and indices lo, hi such that samples[lo:hi] is that slice."""
    if samples.size == 0:
        return samples, 0, 0
    lo = np.searchsorted(samples, a, side="right")
    hi = np.searchsorted(samples, b, side="left")
    lo, hi = int(lo), int(hi)
    return samples[lo:hi], lo, hi


def _maxcrit(
    samples: NDArray[np.float64],
    boundaries: NDArray[np.float64],
    cumulative_probs: NDArray[np.float64],
    merged_coeff: NDArray[np.float64],
    n: int,
    degree: int,
    alpha: float,
    precomputed_coeffs: NDArray[np.float64] | None,
    precomputed_boundaries: NDArray[np.float64] | None,
) -> float:
    """
    Recursive criterion. If precomputed_coeffs and precomputed_boundaries are given and the
    current block is a single interval matching a segment, use the precomputed coefficient.
    """
    left, right = float(boundaries[0]), float(boundaries[-1])
    num_intervals = boundaries.shape[0] - 1

    if precomputed_coeffs is not None and precomputed_boundaries is not None and num_intervals == 1:
        tol = 1e-12
        for idx in range(precomputed_boundaries.shape[0] - 1):
            if (abs(precomputed_boundaries[idx] - left) <= tol and
                    abs(precomputed_boundaries[idx + 1] - right) <= tol):
                piece_coeff = precomputed_coeffs[idx, :]
                break
        else:
            piece_coeff = _coeffint(samples, left, right, n, degree)
    else:
        piece_coeff = _coeffint(samples, left, right, n, degree)

    total_prob = cumulative_probs[-1]
    criterion = _differ(left, right, merged_coeff, piece_coeff) - alpha * np.sqrt(
        (degree + 1) * total_prob * np.log(n) / n
    )

    if num_intervals <= 1:
        return float(criterion)

    half = 0.5 * total_prob
    mask_left = cumulative_probs <= half
    mask_right = cumulative_probs >= half
    boundaries_left = boundaries[mask_left]
    boundaries_right = boundaries[mask_right]
    probs_left = cumulative_probs[mask_left]
    probs_right = cumulative_probs[mask_right] - half

    samples_left, _, _ = _samples_in_interval(samples, boundaries_left[0], boundaries_left[-1])
    samples_right, _, _ = _samples_in_interval(samples, boundaries_right[0], boundaries_right[-1])

    val_left = _maxcrit(
        samples_left, boundaries_left, probs_left, merged_coeff, n, degree, alpha,
        precomputed_coeffs, precomputed_boundaries,
    )
    val_right = _maxcrit(
        samples_right, boundaries_right, probs_right, merged_coeff, n, degree, alpha,
        precomputed_coeffs, precomputed_boundaries,
    )
    return float(max(criterion, val_left + val_right))


def _merge(
    samples: NDArray[np.float64],
    boundaries: NDArray[np.float64],
    cumulative_probs: NDArray[np.float64],
    n: int,
    degree: int,
    alpha: float,
    precomputed_coeffs: NDArray[np.float64] | None,
    precomputed_boundaries: NDArray[np.float64] | None,
) -> tuple[bool, NDArray[np.float64]]:
    """Decide whether to merge the block; return (should_merge, merged_coeff)."""
    left, right = float(boundaries[0]), float(boundaries[-1])
    samples_in_block, _, _ = _samples_in_interval(samples, left, right)
    merged_coeff = _coeffint(samples_in_block, left, right, n, degree)

    total_prob = cumulative_probs[-1]
    half = 0.5 * total_prob
    mask_left = cumulative_probs <= half
    mask_right = cumulative_probs >= half
    boundaries_left = boundaries[mask_left]
    boundaries_right = boundaries[mask_right]
    probs_left = cumulative_probs[mask_left]
    probs_right = cumulative_probs[mask_right] - half

    samples_left, _, _ = _samples_in_interval(samples, boundaries_left[0], boundaries_left[-1])
    samples_right, _, _ = _samples_in_interval(samples, boundaries_right[0], boundaries_right[-1])

    val = _maxcrit(
        samples_left, boundaries_left, probs_left, merged_coeff, n, degree, alpha,
        precomputed_coeffs, precomputed_boundaries,
    ) + _maxcrit(
        samples_right, boundaries_right, probs_right, merged_coeff, n, degree, alpha,
        precomputed_coeffs, precomputed_boundaries,
    )

    return (val <= 0, merged_coeff)


def surf(
    samples: NDArray[np.float64],
    alpha: float,
    degree: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    SURF algorithm: learn piecewise polynomial density on [0, 1].

    Parameters
    ----------
    samples : array of shape (n-1,) with n = 2^k for some k
        Sorted samples strictly in (0, 1).
    alpha : float
        Tuning parameter (larger → fewer pieces).
    degree : int
        Polynomial degree per piece (1 <= degree <= 4).

    Returns
    -------
    boundaries : array
        Interval boundaries (partition of [0, 1]).
    piece_coeffs : array of shape (num_pieces, degree+1)
        Polynomial coefficients for each piece (constant term first).
    """
    samples = np.asarray(samples, dtype=np.float64).ravel()
    n = samples.size + 1
    if n & (n - 1) != 0:
        raise ValueError("len(samples)+1 must be a power of 2")

    if not (1 <= degree <= 4):
        raise ValueError("degree must be in 1..4")

    boundaries = np.concatenate([[0.0], np.sort(samples), [1.0]])
    cumulative_probs = np.linspace(0, 1, n + 1)
    piece_coeffs = np.zeros((n + 1, degree + 1), dtype=np.float64)

    for i in range(n):
        piece_coeffs[i, :] = _coeffint(samples, boundaries[i], boundaries[i + 1], n, degree)

    num_levels = int(np.log2(n))
    for level in range(num_levels):
        for j in range(2 ** (num_levels - 1 - level)):
            prob_start = j * 2 ** (-(num_levels - 1 - level))
            prob_end = (j + 1) * 2 ** (-(num_levels - 1 - level))
            mask = (cumulative_probs >= prob_start) & (cumulative_probs <= prob_end)
            probs_cand = cumulative_probs[mask] - prob_start
            boundaries_cand = boundaries[mask]

            if boundaries_cand.shape[0] <= 1:
                continue

            should_merge, merged_coeff = _merge(
                samples, boundaries_cand, probs_cand, n, degree, alpha,
                precomputed_coeffs=piece_coeffs, precomputed_boundaries=boundaries,
            )

            if should_merge:
                keep = (cumulative_probs <= prob_start) | (cumulative_probs >= prob_end)
                boundaries = boundaries[keep]
                piece_coeffs = piece_coeffs[keep]
                cumulative_probs = cumulative_probs[keep]
                idx_start = np.argmin(np.abs(cumulative_probs - prob_start))
                piece_coeffs[idx_start, :] = merged_coeff

    num_pieces = boundaries.shape[0] - 1
    piece_coeffs = piece_coeffs[:num_pieces, :]
    return boundaries, piece_coeffs
