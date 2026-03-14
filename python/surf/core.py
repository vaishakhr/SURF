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


def _matcoeff_row(n1: float, n2: float, n: int, deg: int) -> NDArray[np.float64]:
    """One row of the coefficient matrix: n * (n2^i - n1^i) / i for i in 1..deg+1."""
    i = np.arange(1, deg + 2, dtype=np.float64)
    return n * (np.power(n2, i) - np.power(n1, i)) / i


def _count_in_interval(samp: NDArray[np.float64], a: float, b: float) -> int:
    """Count samples strictly in (a, b) using binary search. samp must be sorted."""
    if samp.size == 0 or a >= b:
        return 0
    left = np.searchsorted(samp, a, side="right")   # first index with samp[i] > a
    right = np.searchsorted(samp, b, side="left")   # first index with samp[i] >= b
    return max(0, right - left)


def _coeffint(
    samp: NDArray[np.float64],
    I0: float,
    I1: float,
    n: int,
    deg: int,
) -> NDArray[np.float64]:
    """
    Coefficients for polynomial fit on interval [I0, I1].
    Caller may pass samp restricted to [I0, I1] for efficiency; must be sorted.
    """
    len_i = I1 - I0
    if len_i <= 0 or len_i < 1e-12:
        return np.zeros(deg + 1, dtype=np.float64)
    row = OPND1[deg]
    opnd = row[: deg + 2] * len_i + I0

    sampin_i = _count_in_interval(samp, I0, I1)
    normali = len_i * (sampin_i + 1) / (len_i * sampin_i + 1)

    empvec = np.zeros(deg + 1, dtype=np.float64)
    for i in range(deg + 1):
        bin_left, bin_right = opnd[i], opnd[i + 1]
        cnt = _count_in_interval(samp, bin_left, bin_right)
        delta_node = (row[i + 1] - row[i]) / len_i
        empvec[i] = (delta_node + cnt) * normali * (deg + 1)

    matco = np.zeros((deg + 1, deg + 1), dtype=np.float64)
    scale = float(deg + 1)
    for j in range(deg + 1):
        matco[j, :] = scale * _matcoeff_row(opnd[j], opnd[j + 1], n, deg)

    try:
        r = np.linalg.solve(matco, empvec)
    except np.linalg.LinAlgError:
        r = np.zeros(deg + 1, dtype=np.float64)
    return r


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


def _differ(I0: float, I1: float, c1: NDArray[np.float64], c2: NDArray[np.float64]) -> float:
    """L1 difference between polynomials c1 and c2 over [I0, I1] (exact integral)."""
    diff = np.asarray(c1, dtype=np.float64).ravel() - np.asarray(c2, dtype=np.float64).ravel()
    return _integral_abs_poly(I0, I1, diff)


def _samp_slice(samp: NDArray[np.float64], a: float, b: float) -> tuple[NDArray[np.float64], int, int]:
    """Return samp restricted to [a,b] and indices ilo, ihi such that samp[ilo:ihi] is that slice."""
    if samp.size == 0:
        return samp, 0, 0
    ilo = np.searchsorted(samp, a, side="right")
    ihi = np.searchsorted(samp, b, side="left")
    ilo, ihi = int(ilo), int(ihi)
    return samp[ilo:ihi], ilo, ihi


def _maxcrit(
    samp: NDArray[np.float64],
    I: NDArray[np.float64],
    cumuprobI: NDArray[np.float64],
    maincoeff: NDArray[np.float64],
    n: int,
    deg: int,
    alp: float,
    precomp_koi: NDArray[np.float64] | None,
    precomp_I: NDArray[np.float64] | None,
) -> float:
    """
    Recursive criterion. If precomp_koi and precomp_I are given and current I
    is a single interval matching a segment in precomp_I, use precomp_koi instead of coeffint.
    """
    I0, I1 = float(I[0]), float(I[-1])
    num_intervals = I.shape[0] - 1

    if precomp_koi is not None and precomp_I is not None and num_intervals == 1:
        # Look up precomputed coefficient for this interval
        tol = 1e-12
        for idx in range(precomp_I.shape[0] - 1):
            if abs(precomp_I[idx] - I0) <= tol and abs(precomp_I[idx + 1] - I1) <= tol:
                coeff_i = precomp_koi[idx, :]
                break
        else:
            coeff_i = _coeffint(samp, I0, I1, n, deg)
    else:
        coeff_i = _coeffint(samp, I0, I1, n, deg)

    sumprob = cumuprobI[-1]
    differ_i = _differ(I0, I1, maincoeff, coeff_i) - alp * np.sqrt(
        (deg + 1) * sumprob * np.log(n) / n
    )

    if num_intervals <= 1:
        return float(differ_i)

    half = 0.5 * sumprob
    mask1 = cumuprobI <= half
    mask2 = cumuprobI >= half
    I1_arr = I[mask1]
    I2_arr = I[mask2]
    cumuprobI1 = cumuprobI[mask1]
    cumuprobI2 = cumuprobI[mask2] - half

    samp1, _, _ = _samp_slice(samp, I1_arr[0], I1_arr[-1])
    samp2, _, _ = _samp_slice(samp, I2_arr[0], I2_arr[-1])

    val1 = _maxcrit(
        samp1, I1_arr, cumuprobI1, maincoeff, n, deg, alp,
        precomp_koi, precomp_I,
    )
    val2 = _maxcrit(
        samp2, I2_arr, cumuprobI2, maincoeff, n, deg, alp,
        precomp_koi, precomp_I,
    )
    return float(max(differ_i, val1 + val2))


def _merge(
    samp: NDArray[np.float64],
    I: NDArray[np.float64],
    cumuprobI: NDArray[np.float64],
    n: int,
    deg: int,
    alp: float,
    precomp_koi: NDArray[np.float64] | None,
    precomp_I: NDArray[np.float64] | None,
) -> tuple[bool, NDArray[np.float64]]:
    """Decide whether to merge the block I; return (merge?, merged_coeff)."""
    I0, I1 = float(I[0]), float(I[-1])
    samp_in_I, _, _ = _samp_slice(samp, I0, I1)
    maincoeff = _coeffint(samp_in_I, I0, I1, n, deg)

    sumprob = cumuprobI[-1]
    half = 0.5 * sumprob
    mask1 = cumuprobI <= half
    mask2 = cumuprobI >= half
    I1_arr = I[mask1]
    I2_arr = I[mask2]
    cumuprobI1 = cumuprobI[mask1]
    cumuprobI2 = cumuprobI[mask2] - half

    samp1, _, _ = _samp_slice(samp, I1_arr[0], I1_arr[-1])
    samp2, _, _ = _samp_slice(samp, I2_arr[0], I2_arr[-1])

    val = _maxcrit(
        samp1, I1_arr, cumuprobI1, maincoeff, n, deg, alp,
        precomp_koi, precomp_I,
    ) + _maxcrit(
        samp2, I2_arr, cumuprobI2, maincoeff, n, deg, alp,
        precomp_koi, precomp_I,
    )

    return (val <= 0, maincoeff)


def surf(
    samp: NDArray[np.float64],
    alp: float,
    deg: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    SURF algorithm: learn piecewise polynomial density on [0, 1].

    Parameters
    ----------
    samp : array of shape (n-1,) with n = 2^k for some k
        Sorted samples strictly in (0, 1).
    alp : float
        Tuning parameter (e.g. 0.25).
    deg : int
        Polynomial degree per piece (1 <= deg <= 4).

    Returns
    -------
    I : array
        Interval boundaries (partition of [0, 1]).
    koi : array of shape (len(I)-1, deg+1)
        Polynomial coefficients for each piece (constant term first).
    """
    samp = np.asarray(samp, dtype=np.float64).ravel()
    n = samp.size + 1
    if n & (n - 1) != 0:
        raise ValueError("numel(samp)+1 must be a power of 2")

    if not (1 <= deg <= 4):
        raise ValueError("deg must be in 1..4")

    I = np.concatenate([[0.0], np.sort(samp), [1.0]])
    cumuprobI = np.linspace(0, 1, n + 1)
    koi = np.zeros((n + 1, deg + 1), dtype=np.float64)

    for i in range(n):
        koi[i, :] = _coeffint(samp, I[i], I[i + 1], n, deg)

    D = int(np.log2(n))
    for level in range(D):
        for j in range(2 ** (D - 1 - level)):
            probinit = j * 2 ** (-(D - 1 - level))
            probfin = (j + 1) * 2 ** (-(D - 1 - level))
            mask = (cumuprobI >= probinit) & (cumuprobI <= probfin)
            cumuprobIcand = cumuprobI[mask] - probinit
            Icand = I[mask]

            if Icand.shape[0] <= 1:
                continue

            res, koooi = _merge(
                samp, Icand, cumuprobIcand, n, deg, alp,
                precomp_koi=koi, precomp_I=I,
            )

            if res:
                keep = (cumuprobI <= probinit) | (cumuprobI >= probfin)
                I = I[keep]
                koi = koi[keep]
                cumuprobI = cumuprobI[keep]
                idx_probinit = np.argmin(np.abs(cumuprobI - probinit))
                koi[idx_probinit, :] = koooi

    # koi has one extra row (dummy); return only rows for each interval
    koi = koi[: I.shape[0] - 1, :]
    return I, koi
