"""
Run SURF on sample data and print the piecewise polynomial output.
"""

import numpy as np
from surf import surf, regpoly


def format_poly(coeffs, deg):
    """Format coefficients as a0 + a1*x + a2*x^2 + ..."""
    terms = []
    for i in range(min(deg + 1, len(coeffs))):
        c = coeffs[i]
        if i == 0:
            terms.append(f"{c:.6g}")
        else:
            terms.append(f"{c:.6g}*x" + (f"^{i}" if i > 1 else ""))
    return " + ".join(terms)


def main():
    np.random.seed(42)
    n = 2**5 - 1  # 31 samples so output stays readable
    samp = np.sort(np.random.uniform(0.05, 0.95, n))

    print("Input: sorted samples (first 10, ..., last 5)")
    print(samp[:10], "...", samp[-5:])
    print()

    for deg in [1, 2]:
        I, koi = surf(samp, alp=0.25, deg=deg)
        num_pieces = len(I) - 1
        print(f"--- deg = {deg} (alpha = 0.25) ---")
        print(f"Number of pieces: {num_pieces}")
        print()
        for k in range(num_pieces):
            a, b = I[k], I[k + 1]
            coeffs = koi[k, :]
            poly_str = format_poly(coeffs, deg)
            print(f"  Piece {k+1}: x in [{a:.6f}, {b:.6f}]")
            print(f"           p(x) = {poly_str}")
            print()
        print()


if __name__ == "__main__":
    main()
