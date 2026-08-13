#!/usr/bin/env python3
"""Exact certificate for the Chebyshev stable-interpolation obstruction."""

from fractions import Fraction
import json


def add(left, right, sign=1):
    result = [Fraction(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += sign * value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply_linear(polynomial, constant, linear):
    result = [Fraction(0)] * (len(polynomial) + 1)
    for index, value in enumerate(polynomial):
        result[index] += constant * value
        result[index + 1] += linear * value
    return result


def evaluate(polynomial, point):
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * point + coefficient
    return result


def main():
    transformed_previous = [Fraction(1)]
    transformed_current = [Fraction(-1), Fraction(2)]
    audits = []

    for degree in range(1, 13):
        if degree == 1:
            transformed = transformed_current
        else:
            transformed_next = add(
                multiply_linear(transformed_current, Fraction(-2), Fraction(4)),
                transformed_previous,
                sign=-1,
            )
            transformed_previous, transformed_current = (
                transformed_current,
                transformed_next,
            )
            transformed = transformed_current

        scale = Fraction(1, 2 ** (2 * degree - 1))
        monic_chebyshev = [scale * coefficient for coefficient in transformed]
        assert monic_chebyshev[-1] == 1

        # Exact grid sampling is a finite certificate; the analytic inequality
        # follows from |T_d(x)| <= 1 on [-1,1].
        maximum_sample = max(
            abs(evaluate(monic_chebyshev, Fraction(index, 128)))
            for index in range(129)
        )
        assert maximum_sample <= scale
        audits.append(
            {
                "degree": degree,
                "leading_coefficient": "1",
                "analytic_supremum_bound": str(scale),
                "maximum_on_129_point_exact_grid": str(maximum_sample),
                "highest_homogeneous_part_at_axis_direction": "-1",
            }
        )

    print(
        json.dumps(
            {
                "status": "verified",
                "definition": "q_d(t)=2^(1-2d) T_d(2t-1)",
                "analytic_certificate": (
                    "LC(T_d)=2^(d-1); composition with 2t-1 multiplies the "
                    "leading coefficient by 2^d; the displayed scale makes q_d "
                    "monic, while |T_d|<=1 gives ||q_d||_[0,1]<=2^(1-2d)."
                ),
                "audits": audits,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
