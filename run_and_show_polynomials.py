"""
Run SURF on a small random sample and print each piece's polynomial.

Run from repo root: python run_and_show_polynomials.py
"""

import numpy as np
from surf import surf, regpoly


def format_poly(coeffs, degree):
    """Format coefficients as a0 + a1*x + a2*x^2 + ..."""
    terms = []
    for i in range(min(degree + 1, len(coeffs))):
        c = coeffs[i]
        if i == 0:
            terms.append(f"{c:.6g}")
        else:
            terms.append(f"{c:.6g}*x" + (f"^{i}" if i > 1 else ""))
    return " + ".join(terms)


def main():
    np.random.seed(42)
    num_samples = 2**5 - 1  # 31 samples so output stays readable
    samples = np.sort(np.random.uniform(0.05, 0.95, num_samples))

    print("Input: sorted samples (first 10, ..., last 5)")
    print(samples[:10], "...", samples[-5:])
    print()

    for degree in [1, 2]:
        boundaries, piece_coeffs = surf(samples, alpha=0.25, degree=degree)
        num_pieces = len(boundaries) - 1
        print(f"--- degree = {degree} (alpha = 0.25) ---")
        print(f"Number of pieces: {num_pieces}")
        print()
        for k in range(num_pieces):
            left, right = boundaries[k], boundaries[k + 1]
            coeffs = piece_coeffs[k, :]
            poly_str = format_poly(coeffs, degree)
            print(f"  Piece {k+1}: x in [{left:.6f}, {right:.6f}]")
            print(f"           p(x) = {poly_str}")
            print()
        print()


if __name__ == "__main__":
    main()
