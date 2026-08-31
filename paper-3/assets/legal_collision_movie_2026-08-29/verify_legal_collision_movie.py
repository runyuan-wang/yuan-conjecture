#!/usr/bin/env python3
"""Independent verifier for the global collision-movie certificate.

This file does not import the generator.  It rebuilds the key worldline
families, uses a different variable ordering for the contact-wiring system,
and independently exhausts the finite subset optimizations used in the main
claims.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Line:
    b: tuple[F, F, F]
    v: tuple[F, F, F]


@dataclass(frozen=True)
class Hit:
    i: int
    j: int
    t: F
    label: str


@dataclass(frozen=True)
class Case:
    name: str
    lines: tuple[Line, ...]
    hits: tuple[Hit, ...]


def vec(*xs) -> tuple[F, F, F]:
    return tuple(F(x) for x in xs)  # type: ignore[return-value]


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(c, a):
    return tuple(c * x for x in a)


def at(line: Line, t: F):
    return add(line.b, mul(t, line.v))


def exact(case: Case) -> bool:
    return all(at(case.lines[e.i], e.t) == at(case.lines[e.j], e.t) for e in case.hits)


def c12() -> Case:
    vs = (
        vec(F(-3, 10), F(1, 5), F(1, 5)),
        vec(F(1, 5), F(-3, 10), F(-3, 10)),
        vec(F(1, 5), F(-1, 10), F(-1, 5)),
        vec(F(-1, 5), F(-1, 5), F(1, 5)),
        vec(F(79, 198), F(-82, 165), F(-20, 33)),
        vec(F(1, 10), F(-3, 10), F(1, 10)),
    )
    bs = (
        vec(0, 0, 0),
        vec(0, 0, 0),
        vec(0, F(-1, 5), F(-1, 10)),
        vec(F(4, 5), 0, F(-9, 10)),
        vec(F(-329, 330), F(49, 55), F(167, 110)),
        vec(F(-99, 250), F(99, 200), F(99, 1000)),
    )
    ts = (F(0), F(1), F(2), F(3), F(201, 100), F(99, 100))
    return Case(
        "strict_C12_reference",
        tuple(Line(bs[i], vs[i]) for i in range(6)),
        tuple(Hit(i, (i + 1) % 6, ts[i], f"E{i}") for i in range(6)),
    )


def polygon9() -> Case:
    times = [0, 5, 1, 7, 3, 8, 2, 6, 4]
    points = [
        (
            F(times[k], 8),
            vec(k * k + 2 * k + 1, k**3 - 3 * k + 2, 2 * k * k - 5 * k + 7),
        )
        for k in range(9)
    ]
    lines = []
    for i in range(9):
        t0, x0 = points[(i - 1) % 9]
        t1, x1 = points[i]
        v = mul(F(1, 1) / (t1 - t0), sub(x1, x0))
        lines.append(Line(sub(x0, mul(t0, v)), v))
    return Case(
        "generic_rational_C18",
        tuple(lines),
        tuple(Hit(i, (i + 1) % 9, points[i][0], f"E{i}") for i in range(9)),
    )


def hinge8() -> Case:
    lines = [Line(vec(0, 0, 0), vec(0, 0, 0))]
    hits = []
    for j in range(8):
        t = F(j + 1, 9)
        v = vec(j + 1, (j + 1) * (j + 2), 2 * j * j + 3)
        lines.append(Line(mul(-t, v), v))
        hits.append(Hit(0, j + 1, t, f"E{j}"))
    return Case("time_separated_hinge_star", tuple(lines), tuple(hits))


def bush5() -> Case:
    lines = tuple(
        Line(vec(0, 0, 0), vec(i + 1, (i + 1) ** 2 + 1, (i + 1) ** 3 - 2 * i))
        for i in range(5)
    )
    hits = tuple(Hit(i, j, F(0), f"E{i}_{j}") for i, j in itertools.combinations(range(5), 2))
    return Case("simultaneous_bush_K5", lines, hits)


def weave6() -> Case:
    p = [-15, 1, 4, 5, 16, 29]
    q = [19, 18, 7, -10, -16, -18]
    lines = tuple(Line(vec(p[i], 0, 0), vec(q[i] - p[i], 0, 0)) for i in range(6))
    hits = []
    for i, j in itertools.combinations(range(6), 2):
        t = F(p[j] - p[i], 1) / (lines[i].v[0] - lines[j].v[0])
        hits.append(Hit(i, j, t, f"E{i}_{j}"))
    return Case("planar_weave_K6", lines, tuple(sorted(hits, key=lambda e: (e.t, e.label))))


def regulus4(sheared: bool) -> Case:
    ss = [F(1, 20), F(7, 50), F(31, 100), F(13, 20)]
    us = [F(1, 50), F(11, 100), F(2, 5), F(17, 20)]

    def sh(a):
        x, y, z = a
        if not sheared:
            return (x, y, z)
        return (x + F(2, 3) * y + F(1, 5) * z, y + F(1, 4) * z, z)

    lines = []
    for s in ss:
        lines.append(Line(sh(vec(-s, -s * s, 0)), sh(vec(1, 2 * s, 0))))
    for u in us:
        lines.append(Line(sh(vec(u, -u * u, 0)), sh(vec(-1, 2 * u, 0))))
    hits = []
    for i, s in enumerate(ss):
        for j, u in enumerate(us):
            hits.append(Hit(i, 4 + j, (s + u) / 2, f"E{i}_{j}"))
    tag = "sheared_" if sheared else ""
    return Case(f"{tag}regulus_K4_4", tuple(lines), tuple(sorted(hits, key=lambda e: (e.t, e.label))))


def chronology(case: Case, chosen: tuple[int, ...]) -> bool:
    per_line = [[] for _ in case.lines]
    for k in chosen:
        e = case.hits[k]
        per_line[e.i].append(e.t)
        per_line[e.j].append(e.t)
    return all(len(xs) == len(set(xs)) for xs in per_line)


def degrees(case: Case, chosen: tuple[int, ...]) -> list[int]:
    out = [0] * len(case.lines)
    for k in chosen:
        e = case.hits[k]
        out[e.i] += 1
        out[e.j] += 1
    return out


def signs(e: Hit, lane: int) -> tuple[int, int]:
    if lane == e.i:
        return +1, -1
    if lane == e.j:
        return -1, +1
    raise AssertionError


def system(case: Case, chosen: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    """Use [all midpoint coordinates | all times], unlike the generator."""

    m = len(chosen)
    where = {event: local for local, event in enumerate(chosen)}
    omega = {}
    for event in chosen:
        e = case.hits[event]
        g = np.asarray(case.lines[e.i].v, float) - np.asarray(case.lines[e.j].v, float)
        omega[event] = -g / np.linalg.norm(g)

    incident = [[] for _ in case.lines]
    for event in chosen:
        e = case.hits[event]
        incident[e.i].append(event)
        incident[e.j].append(event)
    for lane in range(len(incident)):
        incident[lane].sort(key=lambda k: (case.hits[k].t, case.hits[k].label))

    rows = 3 * sum(max(0, len(xs) - 1) for xs in incident)
    a = np.zeros((rows, 4 * m))
    b = np.zeros(rows)
    r = 0
    for lane, events in enumerate(incident):
        v = np.asarray(case.lines[lane].v, float)
        for early, late in zip(events, events[1:]):
            le, ll = where[early], where[late]
            _in_e, out_e = signs(case.hits[early], lane)
            in_l, _out_l = signs(case.hits[late], lane)
            rhs = -0.5 * (in_l * omega[late] - out_e * omega[early])
            for c in range(3):
                a[r, 3 * ll + c] += 1
                a[r, 3 * le + c] -= 1
                a[r, 3 * m + ll] -= v[c]
                a[r, 3 * m + le] += v[c]
                b[r] = rhs[c]
                r += 1
    return a, b


def audit(a: np.ndarray, b: np.ndarray) -> dict:
    if not len(b):
        return {"rows": 0, "cols": a.shape[1], "rank": 0, "left_null": 0, "residual": 0.0, "ok": True}
    x, _, rank, singular = np.linalg.lstsq(a, b, rcond=1e-12)
    residual = float(np.linalg.norm(a @ x - b) / max(1.0, np.linalg.norm(b)))
    nz = singular[singular > max(a.shape) * np.finfo(float).eps * singular[0]]
    return {
        "rows": a.shape[0],
        "cols": a.shape[1],
        "rank": int(rank),
        "left_null": int(a.shape[0] - rank),
        "residual": residual,
        "sigma_min": float(nz[-1]) if len(nz) else 0.0,
        "ok": residual <= 2e-9,
    }


def feasible(case: Case, chosen: tuple[int, ...], cap: int | None) -> bool:
    if not chronology(case, chosen):
        return False
    if cap is not None and max(degrees(case, chosen), default=0) > cap:
        return False
    a, b = system(case, chosen)
    return bool(audit(a, b)["ok"])


def optimum(case: Case, cap: int | None) -> int:
    m = len(case.hits)
    for size in range(m, -1, -1):
        if any(feasible(case, chosen, cap) for chosen in itertools.combinations(range(m), size)):
            return size
    raise AssertionError


def run() -> dict:
    cases = [c12(), polygon9(), hinge8(), bush5(), weave6(), regulus4(False), regulus4(True)]
    assert all(exact(case) for case in cases)
    full = {}
    for case in cases:
        chosen = tuple(range(len(case.hits)))
        if chronology(case, chosen):
            a, b = system(case, chosen)
            full[case.name] = audit(a, b)
        else:
            full[case.name] = {"ok": False, "reason": "simultaneous reuse"}

    assert full["strict_C12_reference"]["rank"] == 18
    assert full["strict_C12_reference"]["residual"] < 1e-12
    assert full["generic_rational_C18"]["rank"] == 27
    assert full["generic_rational_C18"]["residual"] < 1e-11
    assert full["time_separated_hinge_star"]["rank"] == 21
    assert full["time_separated_hinge_star"]["residual"] < 1e-12
    assert full["planar_weave_K6"]["rank"] == 52
    assert full["planar_weave_K6"]["left_null"] == 20
    assert full["planar_weave_K6"]["residual"] < 1e-10
    assert full["regulus_K4_4"]["rank"] == 59
    assert full["regulus_K4_4"]["left_null"] == 13
    assert full["regulus_K4_4"]["residual"] > 1e-5
    assert full["sheared_regulus_K4_4"]["rank"] == 59
    assert full["sheared_regulus_K4_4"]["residual"] > 1e-4

    optimized = {
        "simultaneous_bush_K5_no_cap": optimum(cases[3], None),
        "planar_weave_K6_cap_2": optimum(cases[4], 2),
        "planar_weave_K6_cap_4": optimum(cases[4], 4),
        "regulus_K4_4_no_cap": optimum(cases[5], None),
        "regulus_K4_4_cap_2": optimum(cases[5], 2),
        "regulus_K4_4_cap_4": optimum(cases[5], 4),
    }
    assert optimized == {
        "simultaneous_bush_K5_no_cap": 2,
        "planar_weave_K6_cap_2": 6,
        "planar_weave_K6_cap_4": 12,
        "regulus_K4_4_no_cap": 12,
        "regulus_K4_4_cap_2": 8,
        "regulus_K4_4_cap_4": 12,
    }

    stability = {}
    for case in (cases[4], cases[5], cases[6]):
        a, b = system(case, tuple(range(len(case.hits))))
        samples = []
        for rcond in (1e-8, 1e-10, 1e-12, 1e-14):
            x, _, rank, _ = np.linalg.lstsq(a, b, rcond=rcond)
            samples.append(
                {
                    "rcond": rcond,
                    "rank": int(rank),
                    "relative_residual": float(np.linalg.norm(a @ x - b) / max(1.0, np.linalg.norm(b))),
                }
            )
        assert len({sample["rank"] for sample in samples}) == 1
        stability[case.name] = samples
    return {
        "status": "PASS",
        "arithmetic": "exact Fraction intersection checks; independent float64 SVD/lstsq wiring",
        "independence": "does not import legal_collision_movie_compiler.py and uses a different variable ordering",
        "full_network_checks": full,
        "exhaustive_optima": optimized,
        "rcond_stability_sweep": stability,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
