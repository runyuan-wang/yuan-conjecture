#!/usr/bin/env python3
"""Verify the conditional-copy / collision-loop proof of the 7/4 bound.

For a finitely supported joint distribution U=(X,Y), put

    C = X+Y,  D = X+2Y,  Z = X-Y.

The Katz--Tao four-atom collision molecule is made probabilistic in two steps:

1. V=(X,Y,Y') where (X,Y) and (X,Y') are conditionally independent given X.
2. T=(V0,V1) where V0 and V1 are conditionally independent given
   W=(X+2Y,Y').

Writing V_i=(a_i,b_i,b'_i), every state of T obeys

    a0+2b0 = a1+2b1,    b'0 = b'1.

The cut/encoding map

    h(T) = (a0+b0, a0+b'0, b1)

determines z1=a1-b'1.  Together, h(T) and q1=(a1,b'1) determine all
of T.  Entropy bookkeeping gives the universal weighted inequality

    4 H(X-Y) <= 2 H(X) + 2 H(Y) + 2 H(X+Y) + H(X+2Y),

and hence the four-slope 7/4 upper bound.

This script constructs the copied distributions exactly (up to floating-point
arithmetic), checks every identity, and stress-tests the inequality.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import log2
import random
from typing import Callable, Hashable, Mapping

from entropy_cycle_lab import CERTIFICATE_26


Atom = tuple[int, int]
VState = tuple[int, int, int]
TState = tuple[VState, VState]
Distribution = dict[Hashable, float]


def normalize(raw: Mapping[Hashable, float]) -> Distribution:
    total = sum(float(p) for p in raw.values() if p > 0)
    if total <= 0:
        raise ValueError("positive total mass required")
    return {state: float(p) / total for state, p in raw.items() if p > 0}


def entropy(dist: Mapping[Hashable, float]) -> float:
    return -sum(p * log2(p) for p in dist.values() if p > 0)


def pushforward(
    dist: Mapping[Hashable, float], fn: Callable[[Hashable], Hashable]
) -> Distribution:
    result: dict[Hashable, float] = defaultdict(float)
    for state, p in dist.items():
        result[fn(state)] += p
    return dict(result)


def joint_pushforward(
    dist: Mapping[Hashable, float],
    *fns: Callable[[Hashable], Hashable],
) -> Distribution:
    return pushforward(dist, lambda state: tuple(fn(state) for fn in fns))


def l1_distance(
    left: Mapping[Hashable, float], right: Mapping[Hashable, float]
) -> float:
    keys = set(left) | set(right)
    return sum(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys)


def conditional_copy_over_x(u: Mapping[Atom, float]) -> dict[VState, float]:
    """Return P(x,y,y')=P(x,y)P(x,y')/P_X(x)."""
    by_x: dict[int, list[tuple[int, float]]] = defaultdict(list)
    px: dict[int, float] = defaultdict(float)
    for (x, y), p in u.items():
        by_x[x].append((y, p))
        px[x] += p

    result: dict[VState, float] = {}
    for x, fibre in by_x.items():
        for y, p in fibre:
            for yp, pp in fibre:
                # Divide first to avoid underflow for very sparse test laws.
                mass = (p / px[x]) * pp
                if mass > 0:
                    result[(x, y, yp)] = mass
    return result


def conditional_copy_over_w(v: Mapping[VState, float]) -> dict[TState, float]:
    """Copy V conditionally independently over W=(x+2y,y')."""
    by_w: dict[tuple[int, int], list[tuple[VState, float]]] = defaultdict(list)
    pw: dict[tuple[int, int], float] = defaultdict(float)
    for state, p in v.items():
        x, y, yp = state
        w = (x + 2 * y, yp)
        by_w[w].append((state, p))
        pw[w] += p

    result: dict[TState, float] = {}
    for w, fibre in by_w.items():
        for v0, p0 in fibre:
            for v1, p1 in fibre:
                # Divide first to avoid underflow for very sparse test laws.
                mass = (p0 / pw[w]) * p1
                if mass > 0:
                    result[(v0, v1)] = mass
    return result


def h_map(t: Hashable) -> tuple[int, int, int]:
    v0, v1 = t  # type: ignore[misc]
    a0, b0, bp0 = v0
    _a1, b1, _bp1 = v1
    return (a0 + b0, a0 + bp0, b1)


def q1_map(t: Hashable) -> Atom:
    _v0, v1 = t  # type: ignore[misc]
    a1, _b1, bp1 = v1
    return (a1, bp1)


def z_from_h(h: tuple[int, int, int]) -> int:
    c0, cp0, b1 = h
    return 2 * c0 - cp0 - 2 * b1


@dataclass(frozen=True)
class Audit:
    name: str
    support_u: int
    support_v: int
    support_t: int
    h_u: float
    h_x: float
    h_y: float
    h_c: float
    h_d: float
    h_z: float
    h_v: float
    h_w: float
    h_t: float
    h_h: float
    h_t_given_h: float
    h_u_given_z: float
    copy_x_error: float
    copy_w_error: float
    q1_marginal_error: float
    reconstruction_entropy_error: float
    reconstruction_collisions: int
    z_formula_failures: int
    weighted_slack: float
    max_ratio: float


def audit(name: str, raw_u: Mapping[Atom, float], tolerance: float = 2e-10) -> Audit:
    u = normalize(raw_u)
    x = pushforward(u, lambda atom: atom[0])
    y = pushforward(u, lambda atom: atom[1])
    c = pushforward(u, lambda atom: atom[0] + atom[1])
    d = pushforward(u, lambda atom: atom[0] + 2 * atom[1])
    z = pushforward(u, lambda atom: atom[0] - atom[1])

    v = conditional_copy_over_x(u)  # type: ignore[arg-type]
    # The proof copies over the pair W=(X+2Y,Y'), not merely X+2Y.
    w = pushforward(v, lambda state: (state[0] + 2 * state[1], state[2]))
    t = conditional_copy_over_w(v)
    h = pushforward(t, h_map)
    q1_dist = pushforward(t, q1_map)
    h_q1 = joint_pushforward(t, h_map, q1_map)

    hu = entropy(u)
    hx, hy, hc, hd, hz = map(entropy, (x, y, c, d, z))
    hv, hw, ht, hh = map(entropy, (v, w, t, h))

    # Verify that (h,q1) uniquely reconstructs T on the support.
    reconstruction: dict[tuple[Hashable, Hashable], TState] = {}
    reconstruction_collisions = 0
    z_formula_failures = 0
    for state in t:
        h_value = h_map(state)
        q1_atom = q1_map(state)
        key = (h_value, q1_atom)
        old = reconstruction.setdefault(key, state)
        if old != state:
            reconstruction_collisions += 1
        if z_from_h(h_value) != q1_atom[0] - q1_atom[1]:
            z_formula_failures += 1

    ht_given_h = ht - hh
    hu_given_z = hu - hz
    hq1_given_h = entropy(h_q1) - hh
    copy_x_error = abs(hv - (2 * hu - hx))
    copy_w_error = abs(ht - (2 * hv - hw))
    q1_marginal_error = l1_distance(q1_dist, u)
    reconstruction_entropy_error = abs(ht_given_h - hq1_given_h)
    weighted_rhs = 2 * hx + 2 * hy + 2 * hc + hd
    weighted_slack = weighted_rhs - 4 * hz
    max_input = max(hx, hy, hc, hd)
    max_ratio = hz / max_input if max_input else 1.0

    # All assertions correspond to explicit steps in the proof.
    assert copy_x_error <= tolerance
    assert copy_w_error <= tolerance
    assert q1_marginal_error <= tolerance
    assert reconstruction_entropy_error <= tolerance
    assert reconstruction_collisions == 0
    assert z_formula_failures == 0
    assert ht_given_h <= hu_given_z + tolerance
    assert ht + tolerance >= 4 * hu - 2 * hx - hd - hy
    assert ht <= 2 * hc + hy + hu_given_z + tolerance
    assert weighted_slack >= -tolerance
    assert max_ratio <= 7 / 4 + tolerance

    return Audit(
        name=name,
        support_u=len(u),
        support_v=len(v),
        support_t=len(t),
        h_u=hu,
        h_x=hx,
        h_y=hy,
        h_c=hc,
        h_d=hd,
        h_z=hz,
        h_v=hv,
        h_w=hw,
        h_t=ht,
        h_h=hh,
        h_t_given_h=ht_given_h,
        h_u_given_z=hu_given_z,
        copy_x_error=copy_x_error,
        copy_w_error=copy_w_error,
        q1_marginal_error=q1_marginal_error,
        reconstruction_entropy_error=reconstruction_entropy_error,
        reconstruction_collisions=reconstruction_collisions,
        z_formula_failures=z_formula_failures,
        weighted_slack=weighted_slack,
        max_ratio=max_ratio,
    )


def print_audit(result: Audit) -> None:
    print(f"\n[{result.name}]")
    print(
        f"supports U/V/T = {result.support_u}/{result.support_v}/{result.support_t}"
    )
    print(
        "H(U), H(Z), H(U|Z) = "
        f"{result.h_u:.12f}, {result.h_z:.12f}, {result.h_u_given_z:.12f}"
    )
    print(
        "H(X), H(Y), H(C), H(D) = "
        f"{result.h_x:.12f}, {result.h_y:.12f}, "
        f"{result.h_c:.12f}, {result.h_d:.12f}"
    )
    print(
        "copy identity errors = "
        f"{result.copy_x_error:.3e}, {result.copy_w_error:.3e}"
    )
    print(
        "q1 marginal / reconstruction entropy errors = "
        f"{result.q1_marginal_error:.3e}, "
        f"{result.reconstruction_entropy_error:.3e}"
    )
    print(
        "H(T|h), H(U|Z) = "
        f"{result.h_t_given_h:.12f}, {result.h_u_given_z:.12f}"
    )
    print(
        "reconstruction collisions / z failures = "
        f"{result.reconstruction_collisions}/{result.z_formula_failures}"
    )
    print(
        "weighted slack, H(Z)/max input = "
        f"{result.weighted_slack:.12f}, {result.max_ratio:.12f}"
    )


def certificate_distribution() -> dict[Atom, float]:
    return {(x, y): p for x, y, p in CERTIFICATE_26}


def random_distribution(
    rng: random.Random, side: int, concentration_log10: float
) -> dict[Atom, float]:
    alpha = 10.0**concentration_log10
    result: dict[Atom, float] = {}
    for x in range(side):
        for y in range(side):
            # Gamma draws give a Dirichlet distribution after normalization.
            result[(x, y)] = rng.gammavariate(alpha, 1.0)
    return result


def stress_test(trials: int = 500, seed: int = 20260811) -> tuple[float, Audit]:
    rng = random.Random(seed)
    worst: Audit | None = None
    for trial in range(trials):
        side = rng.choice((2, 3, 4))
        concentration = rng.uniform(-2.5, 1.0)
        result = audit(
            f"random-{trial}", random_distribution(rng, side, concentration)
        )
        if worst is None or result.weighted_slack < worst.weighted_slack:
            worst = result
    assert worst is not None
    return worst.weighted_slack, worst


def main() -> None:
    star = {(0, y): 1.0 for y in range(5)}
    dense = {
        (x, y): (1 + x + 2 * y) ** 1.3
        for x in range(4)
        for y in range(4)
    }

    for result in (
        audit("single-X star", star),
        audit("26-point public certificate", certificate_distribution()),
        audit("dense non-injective 4x4", dense),
    ):
        print_audit(result)

    minimum_slack, worst = stress_test()
    print("\n[random stress test]")
    print("trials = 500, deterministic seed = 20260811")
    print(f"minimum weighted slack = {minimum_slack:.12e}")
    print(
        "worst trial H(Z)/max input = "
        f"{worst.max_ratio:.12f}; all proof assertions passed"
    )


if __name__ == "__main__":
    main()
