#!/usr/bin/env python3
"""Exact certificate for the strict-time six-worldline C12 counterexample.

All arithmetic is rational.  The program verifies:
  * the six prescribed consecutive intersections and their strict time order;
  * positive separation of every nonconsecutive pair of spacetime segments;
  * positive separation of every event from every nonincident segment;
  * general position of the six velocity vectors.

The segment distance minimizers are found by checking the unconstrained
quadratic minimizer and every face of the parameter square [0, 3]^2.
"""

from fractions import Fraction as F
from itertools import combinations
import json


def V(*xs):
    return tuple(F(x) for x in xs)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(s, a):
    return tuple(s * x for x in a)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def det3(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def clip(x, lo=F(0), hi=F(3)):
    return max(lo, min(hi, x))


def rat(x):
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def spacetime_point(c, v, t):
    return (t,) + add(c, scale(t, v))


def segment_distance_sq(c1, v1, c2, v2):
    a1, a2 = (F(0),) + c1, (F(0),) + c2
    d1, d2 = (F(1),) + v1, (F(1),) + v2
    r = sub(a1, a2)
    aa, bb, ab = dot(d1, d1), dot(d2, d2), dot(d1, d2)
    ar, br = dot(d1, r), dot(d2, r)
    candidates = []

    determinant = aa * bb - ab * ab
    if determinant:
        # [aa,-ab;-ab,bb] [t,s] = [-ar,br]
        t = ((-ar) * bb + ab * br) / determinant
        s = (aa * br - ab * ar) / determinant
        if F(0) <= t <= F(3) and F(0) <= s <= F(3):
            candidates.append((t, s))

    for t in (F(0), F(3)):
        s = clip(dot(d2, add(r, scale(t, d1))) / bb)
        candidates.append((t, s))
    for s in (F(0), F(3)):
        t = clip(-dot(d1, sub(r, scale(s, d2))) / aa)
        candidates.append((t, s))

    def value(ts):
        t, s = ts
        delta = sub(add(r, scale(t, d1)), scale(s, d2))
        return dot(delta, delta)

    t, s = min(candidates, key=value)
    return value((t, s)), t, s


def point_segment_distance_sq(point, c, v):
    a, d = (F(0),) + c, (F(1),) + v
    r = sub(point, a)
    t = clip(dot(r, d) / dot(d, d))
    delta = sub(r, scale(t, d))
    return dot(delta, delta), t


velocities = [
    V(F(-3, 10), F(1, 5), F(1, 5)),
    V(F(1, 5), F(-3, 10), F(-3, 10)),
    V(F(1, 5), F(-1, 10), F(-1, 5)),
    V(F(-1, 5), F(-1, 5), F(1, 5)),
    V(F(79, 198), F(-82, 165), F(-20, 33)),
    V(F(1, 10), F(-3, 10), F(1, 10)),
]

intercepts = [
    V(0, 0, 0),
    V(0, 0, 0),
    V(0, F(-1, 5), F(-1, 10)),
    V(F(4, 5), 0, F(-9, 10)),
    V(F(-329, 330), F(49, 55), F(167, 110)),
    V(F(-99, 250), F(99, 200), F(99, 1000)),
]

# Event i is the intersection of line i and line (i+1) modulo six.
event_times = [F(0), F(1), F(2), F(3), F(201, 100), F(99, 100)]
events = []
for i, t in enumerate(event_times):
    j = (i + 1) % 6
    p_i = spacetime_point(intercepts[i], velocities[i], t)
    p_j = spacetime_point(intercepts[j], velocities[j], t)
    assert p_i == p_j, (i, p_i, p_j)
    events.append(p_i)

assert event_times[0] < event_times[5] < event_times[1] < event_times[2] < event_times[4] < event_times[3]

adjacent_pairs = {tuple(sorted((i, (i + 1) % 6))) for i in range(6)}
nonadjacent = []
for i, j in combinations(range(6), 2):
    if (i, j) not in adjacent_pairs:
        d2, t, s = segment_distance_sq(intercepts[i], velocities[i], intercepts[j], velocities[j])
        nonadjacent.append((d2, i + 1, j + 1, t, s))

nonincident = []
for qi, point in enumerate(events):
    incident = {qi, (qi + 1) % 6}
    for line in range(6):
        if line not in incident:
            d2, t = point_segment_distance_sq(point, intercepts[line], velocities[line])
            nonincident.append((d2, qi, line + 1, t))

triple_dets = []
for inds in combinations(range(6), 3):
    d = det3(*(velocities[i] for i in inds))
    triple_dets.append((abs(d), tuple(i + 1 for i in inds), d))

affine_dets = []
for inds in combinations(range(6), 4):
    i, j, k, ell = inds
    d = det3(
        sub(velocities[j], velocities[i]),
        sub(velocities[k], velocities[i]),
        sub(velocities[ell], velocities[i]),
    )
    affine_dets.append((abs(d), tuple(x + 1 for x in inds), d))

min_pair = min(nonadjacent)
min_event = min(nonincident)
min_triple = min(triple_dets)
min_affine = min(affine_dets)

assert min_pair[0] == F(2017, 55600), min_pair
assert min_event[0] == F(529, 12200), min_event
assert min_triple[0] == F(7, 3000), min_triple
assert min_affine[0] == F(43, 33000), min_affine

result = {
    "status": "verified",
    "strict_event_order": ["Q0", "Q5", "Q1", "Q2", "Q4", "Q3"],
    "event_times": [rat(t) for t in event_times],
    "minimum_nonadjacent_segment_distance_squared": rat(min_pair[0]),
    "minimum_nonadjacent_segment_pair": [min_pair[1], min_pair[2]],
    "minimum_event_to_nonincident_segment_distance_squared": rat(min_event[0]),
    "minimum_event_to_nonincident_line": [f"Q{min_event[1]}", min_event[2]],
    "minimum_absolute_triple_velocity_determinant": rat(min_triple[0]),
    "minimum_absolute_affine_quadruple_velocity_determinant": rat(min_affine[0]),
    "consequence": (
        "For every delta < 1/25, the six event balls give an induced C12; "
        "there is no nonconsecutive tube intersection, no three-tube bush, "
        "and no common affine velocity 2-plane or spacetime 3-hyperplane."
    ),
}

print(json.dumps(result, indent=2, ensure_ascii=False))
