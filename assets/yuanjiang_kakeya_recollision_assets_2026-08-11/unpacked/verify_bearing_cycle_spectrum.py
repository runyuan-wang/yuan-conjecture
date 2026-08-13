#!/usr/bin/env python3
"""Bearing-rigidity spectrum for the strict-time six-worldline C12.

For an encounter edge ij, exact synchronous intersection implies
P_(v_i-v_j)^perp (b_i-b_j)=0.  The stacked operator always kills three
translations and the global time-shift mode h_i=v_i.  A bare C6 has two
additional generic flexes; one generic chord kills them.
"""

from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np


def vec(*xs: F | int) -> np.ndarray:
    return np.array([float(F(x)) for x in xs], dtype=float)


v = np.array(
    [
        vec(F(-3, 10), F(1, 5), F(1, 5)),
        vec(F(1, 5), F(-3, 10), F(-3, 10)),
        vec(F(1, 5), F(-1, 10), F(-1, 5)),
        vec(F(-1, 5), F(-1, 5), F(1, 5)),
        vec(F(79, 198), F(-82, 165), F(-20, 33)),
        vec(F(1, 10), F(-3, 10), F(1, 10)),
    ]
)
b = np.array(
    [
        vec(0, 0, 0),
        vec(0, 0, 0),
        vec(0, F(-1, 5), F(-1, 10)),
        vec(F(4, 5), 0, F(-9, 10)),
        vec(F(-329, 330), F(49, 55), F(167, 110)),
        vec(F(-99, 250), F(99, 200), F(99, 1000)),
    ]
)


def bearing_matrix(edges: list[tuple[int, int]]) -> np.ndarray:
    rows = []
    for i, j in edges:
        g = v[i] - v[j]
        g = g / np.linalg.norm(g)
        p = np.eye(3) - np.outer(g, g)
        block = np.zeros((3, 18))
        block[:, 3 * i : 3 * i + 3] = p
        block[:, 3 * j : 3 * j + 3] = -p
        rows.append(block)
    return np.concatenate(rows, axis=0)


def audit(edges: list[tuple[int, int]], tested_intercepts: np.ndarray) -> dict:
    r = bearing_matrix(edges)
    singular = np.linalg.svd(r, compute_uv=False)
    rank = int(np.sum(singular > 1e-9))
    nullity = 18 - rank
    residual = r @ tested_intercepts.reshape(-1)

    translation = np.zeros((18, 3))
    for i in range(6):
        translation[3 * i : 3 * i + 3, :] = np.eye(3)
    gauge = np.column_stack([translation, v.reshape(-1)])

    components = 1
    cycle_rank = len(edges) - 6 + components
    edge_stalk_dimension = 2 * len(edges)
    h1_dimension = edge_stalk_dimension - rank
    forest_baseline_rank = 2 * (6 - components)
    fresh_cycle_closure_rank = rank - forest_baseline_rank

    return {
        "edges": [list(e) for e in edges],
        "matrix_shape": list(r.shape),
        "rank": rank,
        "nullity": nullity,
        "trivial_gauge_dimension": int(np.linalg.matrix_rank(gauge)),
        "extra_flex_dimension": nullity - int(np.linalg.matrix_rank(gauge)),
        "graph_cycle_rank": cycle_rank,
        "edge_stalk_dimension": edge_stalk_dimension,
        "H0_dimension": nullity,
        "H1_dimension": h1_dimension,
        "forest_baseline_rank": forest_baseline_rank,
        "fresh_cycle_closure_rank": fresh_cycle_closure_rank,
        "closure_rank_plus_H1": fresh_cycle_closure_rank + h1_dimension,
        "twice_graph_cycle_rank": 2 * cycle_rank,
        "intercept_constraint_residual": float(np.linalg.norm(residual)),
        "nonzero_singular_values": [
            float(x) for x in singular if x > 1e-9
        ],
        "translation_gauge_residual": float(np.linalg.norm(r @ translation)),
        "time_shift_gauge_residual": float(np.linalg.norm(r @ v.reshape(-1))),
    }


cycle = [(i, (i + 1) % 6) for i in range(6)]
cycle_audit = audit(cycle, b)

# A generic algebraic chord need not be an actual encounter.  It is included
# only to verify the rigidity rank statement: two scalar flex modes disappear.
chord_audit = audit(cycle + [(0, 3)], b)
two_chord_audit = audit(cycle + [(0, 3), (1, 4)], b)

assert cycle_audit["rank"] == 12
assert cycle_audit["extra_flex_dimension"] == 2
assert cycle_audit["intercept_constraint_residual"] < 1e-12
assert chord_audit["rank"] == 14
assert chord_audit["extra_flex_dimension"] == 0
assert two_chord_audit["rank"] == 14
assert two_chord_audit["H1_dimension"] == 2
assert cycle_audit["translation_gauge_residual"] < 1e-12
assert cycle_audit["time_shift_gauge_residual"] < 1e-12
for item in (cycle_audit, chord_audit, two_chord_audit):
    assert item["closure_rank_plus_H1"] == item["twice_graph_cycle_rank"]

certificate = {
    "status": "verified",
    "verifier": Path(__file__).name,
    "exact_C6_encounter_network": cycle_audit,
    "cycle_plus_generic_chord_rank_test": chord_audit,
    "cycle_plus_two_generic_chords_stress_test": two_chord_audit,
    "scope_warning": (
        "The extra chord is not an encounter in the strict C12 geometry; it "
        "only certifies the generic bearing-rank mechanism.  A Kakeya proof "
        "must extract such independent constraints from actual near events."
    ),
}

Path("bearing_cycle_spectrum_certificate.json").write_text(
    json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(certificate, indent=2))
