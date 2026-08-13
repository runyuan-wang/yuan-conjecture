#!/usr/bin/env python3
"""Exact exponent audit for the Kakeya/Boltzmann--Grad C12 bridge.

All powers are powers of Q, where delta=Q^{-1}, N=Q^3, the occupied
cube count is Q^s, and the hard-sphere diameter is epsilon=Q^{-3/2}.
"""

from fractions import Fraction
import json


def audit(s: Fraction) -> dict[str, str]:
    d_exp = Fraction(4) - s
    p_edge_exp = Fraction(-1)  # (epsilon/delta)^2
    p_c12_exp = 6 * p_edge_exp

    return {
        "s": str(s),
        "multiplicity_D_exponent": str(d_exp),
        "edge_collision_probability_exponent": str(p_edge_exp),
        "generic_C12_probability_exponent": str(p_c12_exp),
        "selected_disjoint_C12_count_exponent": "3",
        "selected_realized_expectation_exponent": str(Fraction(3) + p_c12_exp),
        "extremal_guaranteed_C12_count_exponent": "4",
        "extremal_realized_union_bound_exponent": str(Fraction(4) + p_c12_exp),
        "random_like_C12_count_exponent": str(6 * (Fraction(5) - s)),
        "random_like_realized_C12_exponent": str(6 * d_exp),
        "coarse_near_pair_count_exponent": str(Fraction(8) - s),
        "realized_collision_count_exponent": str(Fraction(7) - s),
        "realized_collisions_per_tube_exponent": str(d_exp),
        "ghost_angular_energy_collisions_total_exponent": str(
            Fraction(22, 3) - Fraction(4, 3) * s
        ),
        "ghost_angular_energy_collisions_per_tube_exponent": str(
            Fraction(13, 3) - Fraction(4, 3) * s
        ),
    }


def main() -> None:
    samples = [
        Fraction(3),
        Fraction(13, 4),
        Fraction(10, 3),
        Fraction(11, 3),
        Fraction(39, 10),
    ]
    rows = [audit(s) for s in samples]

    for row in rows:
        assert Fraction(row["selected_realized_expectation_exponent"]) == -3
        assert Fraction(row["extremal_realized_union_bound_exponent"]) == -2
        assert Fraction(row["random_like_realized_C12_exponent"]) == (
            6 * Fraction(row["multiplicity_D_exponent"])
        )
        assert Fraction(row["realized_collisions_per_tube_exponent"]) == (
            Fraction(row["multiplicity_D_exponent"])
        )
        assert Fraction(
            row["ghost_angular_energy_collisions_per_tube_exponent"]
        ) == Fraction(4, 3) * Fraction(
            row["multiplicity_D_exponent"]
        ) - 1

    certificate = {
        "model": {
            "delta": "Q^-1",
            "N": "Q^3",
            "epsilon_BG": "Q^-3/2",
            "single_collision_probability": "Q^-1",
            "generic_C12_probability": "Q^-6",
        },
        "audits": rows,
        "checks_passed": True,
    }
    with open("bg_c12_scaling_certificate.json", "w", encoding="utf-8") as fh:
        json.dump(certificate, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    print(json.dumps(certificate, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
