#!/usr/bin/env python3
"""Numerical certificate for thickening an exact point-collision C12 into
a finite-radius hard-sphere wiring diagram.

The six affine lanes are the strict-time rational example from
``verify_strict_time_c12_worldlines.py``.  At each prescribed crossing,
equal-mass spheres make a central collision and exchange velocities.  We
solve the linear matching equations for the six contact midpoints and times,
then verify that the outgoing half-ray on every lane meets the next incoming
half-ray and that no unintended overlap occurs on the observed time window.

This is a floating-point *existence certificate for this finite example*, not
a general symbolic hard-sphere theorem.  The matching system and all reported
residuals are explicit and independently reproducible.
"""

from __future__ import annotations

import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import numpy as np


def vec(*xs: F | int) -> np.ndarray:
    return np.array([float(F(x)) for x in xs], dtype=float)


velocities = np.array(
    [
        vec(F(-3, 10), F(1, 5), F(1, 5)),
        vec(F(1, 5), F(-3, 10), F(-3, 10)),
        vec(F(1, 5), F(-1, 10), F(-1, 5)),
        vec(F(-1, 5), F(-1, 5), F(1, 5)),
        vec(F(79, 198), F(-82, 165), F(-20, 33)),
        vec(F(1, 10), F(-3, 10), F(1, 10)),
    ]
)

intercepts = np.array(
    [
        vec(0, 0, 0),
        vec(0, 0, 0),
        vec(0, F(-1, 5), F(-1, 10)),
        vec(F(4, 5), 0, F(-9, 10)),
        vec(F(-329, 330), F(49, 55), F(167, 110)),
        vec(F(-99, 250), F(99, 200), F(99, 1000)),
    ]
)

# Event e joins lane e to lane e+1 (mod 6).
event_times = np.array([0.0, 1.0, 2.0, 3.0, 2.01, 0.99])
event_points = np.array(
    [intercepts[e] + event_times[e] * velocities[e] for e in range(6)]
)


def incident_signs(event: int, lane: int) -> tuple[int, int]:
    """Return (incoming sign, outgoing sign) at an event for one lane.

    For event e, set g=v_e-v_{e+1} and omega=-g/|g|.  Incoming centers are
    m+epsilon*omega/2 on lane e and m-epsilon*omega/2 on lane e+1.  A central
    equal-mass collision swaps velocities, hence the outgoing sign on a fixed
    lane is the opposite incoming sign.
    """

    if lane == event:
        return +1, -1
    if lane == (event + 1) % 6:
        return -1, +1
    raise ValueError("lane is not incident to event")


relative = velocities - np.roll(velocities, -1, axis=0)
omega = -relative / np.linalg.norm(relative, axis=1)[:, None]

# Unknown vector z=(y_0,...,y_5,s_0,...,s_5), where contact midpoint
# m_e=q_e+epsilon*y_e and collision time tau_e=t_e+epsilon*s_e.
# We solve at epsilon=1; scaling then gives every sufficiently small epsilon.
A = np.zeros((18, 24), dtype=float)
b = np.zeros(18, dtype=float)

for lane in range(6):
    events = [(lane - 1) % 6, lane]
    early, late = sorted(events, key=lambda e: event_times[e])
    _, out_sign = incident_signs(early, lane)
    in_sign, _ = incident_signs(late, lane)
    rhs = -0.5 * (in_sign * omega[late] - out_sign * omega[early])

    for coord in range(3):
        row = 3 * lane + coord
        A[row, 3 * late + coord] += 1.0
        A[row, 3 * early + coord] -= 1.0
        A[row, 18 + late] -= velocities[lane, coord]
        A[row, 18 + early] += velocities[lane, coord]
        b[row] = rhs[coord]

solution, _, rank, singular_values = np.linalg.lstsq(A, b, rcond=None)
y = solution[:18].reshape(6, 3)
s = solution[18:]
linear_residual = A @ solution - b

EPSILON = 1.0e-5
midpoints = event_points + EPSILON * y
tau = event_times + EPSILON * s


def contact(event: int, lane: int, incoming: bool) -> np.ndarray:
    incoming_sign, outgoing_sign = incident_signs(event, lane)
    sign = incoming_sign if incoming else outgoing_sign
    return midpoints[event] + 0.5 * EPSILON * sign * omega[event]


lane_data = []
matching_residuals = []
for lane in range(6):
    e0, e1 = (lane - 1) % 6, lane
    early, late = sorted((e0, e1), key=lambda e: tau[e])
    start = contact(early, lane, incoming=False)
    end = contact(late, lane, incoming=True)
    residual = end - start - (tau[late] - tau[early]) * velocities[lane]
    matching_residuals.append(residual)
    lane_data.append((early, late, start, end))


def lane_affine_piece(lane: int, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Return a,v with center=a+t*v on the lane occupied at non-event t."""

    early, late, start, _ = lane_data[lane]
    v = velocities[lane]
    if t < tau[early]:
        p = contact(early, lane, incoming=True)
        return p - tau[early] * v, v
    if t < tau[late]:
        return start - tau[early] * v, v
    p = contact(late, lane, incoming=False)
    return p - tau[late] * v, v


def interval_pair_minimum(a: np.ndarray, g: np.ndarray, lo: float, hi: float):
    if np.dot(g, g) == 0.0:
        t = lo
    else:
        t = float(np.clip(-np.dot(a, g) / np.dot(g, g), lo, hi))
    d = a + t * g
    return float(np.linalg.norm(d)), t


# Check all pair distances on every open interval between collision times.
breaks = [0.0] + sorted(float(t) for t in tau if 0.0 < t < 3.0) + [3.0]
pair_minima = []
for lo, hi in zip(breaks[:-1], breaks[1:]):
    if hi - lo < 1.0e-12:
        continue
    probe = (lo + hi) / 2.0
    pieces = [lane_affine_piece(i, probe) for i in range(6)]
    for i, j in combinations(range(6), 2):
        a = pieces[i][0] - pieces[j][0]
        g = pieces[i][1] - pieces[j][1]
        dist, tmin = interval_pair_minimum(a, g, lo, hi)
        pair_minima.append((dist, i, j, tmin, lo, hi))

# At each collision, verify central approach and exact velocity exchange.
collision_checks = []
for event in range(6):
    i, j = event, (event + 1) % 6
    xi = contact(event, i, incoming=True)
    xj = contact(event, j, incoming=True)
    separation = xi - xj
    g = velocities[i] - velocities[j]
    normal = separation / np.linalg.norm(separation)
    vi_out = velocities[i] - normal * np.dot(normal, g)
    vj_out = velocities[j] + normal * np.dot(normal, g)
    collision_checks.append(
        {
            "event": event,
            "contact_distance": float(np.linalg.norm(separation)),
            "approach_dot": float(np.dot(separation, g)),
            "swap_error_i": float(np.linalg.norm(vi_out - velocities[j])),
            "swap_error_j": float(np.linalg.norm(vj_out - velocities[i])),
        }
    )

unintended = [x for x in pair_minima if x[0] < EPSILON * (1.0 - 1.0e-7)]
minimum_piece_distance = min(x[0] for x in pair_minima)

assert rank == 18
assert np.max(np.abs(linear_residual)) < 1.0e-12
assert np.max(np.abs(matching_residuals)) < 1.0e-12
assert list(np.argsort(tau)) == [0, 5, 1, 2, 4, 3]
assert not unintended
for check in collision_checks:
    assert abs(check["contact_distance"] - EPSILON) < 1.0e-12
    assert check["approach_dot"] < 0.0
    # Contact points are separated by only 1e-5, so forming the unit normal
    # loses a few decimal digits in binary floating point.
    assert check["swap_error_i"] < 1.0e-10
    assert check["swap_error_j"] < 1.0e-10

certificate = {
    "status": "verified",
    "verifier": Path(__file__).name,
    "model": "six equal-mass hard spheres with central binary collisions",
    "sphere_diameter": EPSILON,
    "matching_matrix_shape": list(A.shape),
    "matching_matrix_rank": int(rank),
    "smallest_nonzero_singular_value": float(singular_values[-1]),
    "maximum_linear_system_residual": float(np.max(np.abs(linear_residual))),
    "maximum_lane_matching_residual": float(np.max(np.abs(matching_residuals))),
    "maximum_midpoint_shift_in_diameter_units": float(np.max(np.linalg.norm(y, axis=1))),
    "maximum_time_shift_in_diameter_units": float(np.max(np.abs(s))),
    "strict_collision_order": [int(i) for i in np.argsort(tau)],
    "minimum_center_distance_over_affine_time_pieces": minimum_piece_distance,
    "unintended_overlap_count": len(unintended),
    "collision_checks": collision_checks,
    "certified_consequence": (
        "The strict-time six-line C12 has a finite-radius central-collision "
        "realization for the tested diameter. Every collision exchanges the "
        "two velocities; the six lane segments match to numerical precision; "
        "and no unintended overlap occurs on 0<=t<=3."
    ),
    "scope_warning": (
        "This certifies one isolated bounded-size gadget. It does not embed an "
        "arbitrary Kakeya tube family into one global hard-sphere history."
    ),
}

output = Path(__file__).with_name("hard_sphere_c12_wiring_certificate.json")
output.write_text(json.dumps(certificate, indent=2) + "\n")
print(json.dumps(certificate, indent=2))
