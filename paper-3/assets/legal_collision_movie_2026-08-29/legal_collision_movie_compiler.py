#!/usr/bin/env python3
"""Global finite-radius collision-movie compilation audit.

The input is a finite family of genuine affine worldlines

    x_i(t) = b_i + t v_i,  b_i,v_i in R^3,

together with exact pairwise spacetime intersections.  Each selected
intersection is interpreted as a central equal-mass collision, so the two
velocity tokens swap lanes.  Every event is locally legal.  The non-local
question is whether all contact offsets can be chosen consistently along the
shared original lanes when the point particles are thickened to spheres of
diameter epsilon.

At event e, write its contact midpoint and time as

    m_e = q_e + epsilon*y_e,  tau_e = t_e + epsilon*s_e.

Continuity along every lane between consecutive selected events gives a
linear system A(y,s)=c(omega).  The audit reports its rank, left-null
compatibility dimension, residual, conditioning, and the largest event subset
that is simultaneously wireable.  The code also applies explicit congestion
budgets to distinguish a physical movie from a low-reuse Kakeya extraction.

This is a finite mechanism test, not a Kakeya theorem.  It deliberately does
not compute local concurrence rank, velocity rectangles, graph spectra, or
artificial collision networks.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class Lane:
    name: str
    root: tuple[Q, Q, Q]
    velocity: tuple[Q, Q, Q]


@dataclass(frozen=True)
class Event:
    name: str
    lane_i: int
    lane_j: int
    time: Q
    point: tuple[Q, Q, Q]


@dataclass(frozen=True)
class Family:
    name: str
    geometry_class: str
    lanes: tuple[Lane, ...]
    events: tuple[Event, ...]
    notes: str


def qv(*xs: int | Q) -> tuple[Q, Q, Q]:
    return tuple(Q(x) for x in xs)  # type: ignore[return-value]


def qadd(a: Sequence[Q], b: Sequence[Q]) -> tuple[Q, ...]:
    return tuple(x + y for x, y in zip(a, b))


def qsub(a: Sequence[Q], b: Sequence[Q]) -> tuple[Q, ...]:
    return tuple(x - y for x, y in zip(a, b))


def qscale(c: Q, a: Sequence[Q]) -> tuple[Q, ...]:
    return tuple(c * x for x in a)


def qposition(lane: Lane, time: Q) -> tuple[Q, Q, Q]:
    return qadd(lane.root, qscale(time, lane.velocity))  # type: ignore[return-value]


def as_float(a: Sequence[Q]) -> np.ndarray:
    return np.asarray([float(x) for x in a], dtype=float)


def exact_intersection(lane_i: Lane, lane_j: Lane, time: Q) -> tuple[Q, Q, Q]:
    left = qposition(lane_i, time)
    right = qposition(lane_j, time)
    if left != right:
        raise ValueError(f"declared event is not exact: {lane_i.name}, {lane_j.name}, t={time}")
    return left


def make_event(name: str, i: int, j: int, time: Q, lanes: Sequence[Lane]) -> Event:
    return Event(name, i, j, time, exact_intersection(lanes[i], lanes[j], time))


def strict_c12_family() -> Family:
    velocities = (
        qv(Q(-3, 10), Q(1, 5), Q(1, 5)),
        qv(Q(1, 5), Q(-3, 10), Q(-3, 10)),
        qv(Q(1, 5), Q(-1, 10), Q(-1, 5)),
        qv(Q(-1, 5), Q(-1, 5), Q(1, 5)),
        qv(Q(79, 198), Q(-82, 165), Q(-20, 33)),
        qv(Q(1, 10), Q(-3, 10), Q(1, 10)),
    )
    roots = (
        qv(0, 0, 0),
        qv(0, 0, 0),
        qv(0, Q(-1, 5), Q(-1, 10)),
        qv(Q(4, 5), 0, Q(-9, 10)),
        qv(Q(-329, 330), Q(49, 55), Q(167, 110)),
        qv(Q(-99, 250), Q(99, 200), Q(99, 1000)),
    )
    times = (Q(0), Q(1), Q(2), Q(3), Q(201, 100), Q(99, 100))
    lanes = tuple(Lane(f"L{i}", roots[i], velocities[i]) for i in range(6))
    events = tuple(make_event(f"E{i}", i, (i + 1) % 6, times[i], lanes) for i in range(6))
    return Family(
        "strict_C12_reference",
        "wide generic fixed cycle",
        lanes,
        events,
        "Previously verified positive control; included only to validate the generalized compiler.",
    )


def generic_polygon_cycle(n: int = 9) -> Family:
    """A rational exact cycle made from generic spacetime event points."""

    # The non-monotone time labels prevent the polygon order from becoming a
    # physical directed cycle.  Spatial coordinates are deterministic integer
    # polynomials chosen to avoid accidental collinearities.
    permutation = [0, 5, 1, 7, 3, 8, 2, 6, 4]
    if n != len(permutation):
        raise ValueError("this audited generic cycle uses n=9")
    points: list[tuple[Q, tuple[Q, Q, Q]]] = []
    for k, tk in enumerate(permutation):
        t = Q(tk, 8)
        x = qv(k * k + 2 * k + 1, k * k * k - 3 * k + 2, 2 * k * k - 5 * k + 7)
        points.append((t, x))

    lanes: list[Lane] = []
    for i in range(n):
        t0, x0 = points[(i - 1) % n]
        t1, x1 = points[i]
        velocity = qscale(Q(1, 1) / (t1 - t0), qsub(x1, x0))
        root = qsub(x0, qscale(t0, velocity))
        lanes.append(Lane(f"L{i}", root, velocity))
    lane_tuple = tuple(lanes)
    events = tuple(
        make_event(f"E{i}", i, (i + 1) % n, points[i][0], lane_tuple)
        for i in range(n)
    )
    return Family(
        "generic_rational_C18",
        "wide generic fixed cycle",
        lane_tuple,
        events,
        "Nine genuine rational R4 lines; every event is an exact intersection of consecutive lanes.",
    )


def time_hinge_star(k: int = 8) -> Family:
    """One shared lane meets k generic transverse lanes at separated times."""

    lanes: list[Lane] = [Lane("spine", qv(0, 0, 0), qv(0, 0, 0))]
    events: list[Event] = []
    for j in range(k):
        t = Q(j + 1, k + 1)
        velocity = qv(j + 1, (j + 1) * (j + 2), 2 * j * j + 3)
        root = qscale(-t, velocity)
        lanes.append(Lane(f"leaf_{j}", root, velocity))
    lane_tuple = tuple(lanes)
    for j in range(k):
        t = Q(j + 1, k + 1)
        events.append(make_event(f"E{j}", 0, j + 1, t, lane_tuple))
    return Family(
        "time_separated_hinge_star",
        "hinge tree / shared lane",
        lane_tuple,
        tuple(events),
        "All collisions on the shared lane are strictly time separated; leaf directions are non-coplanar.",
    )


def simultaneous_bush(n: int = 5) -> Family:
    """All lines pass through one event; candidate collisions are all pairs."""

    lanes = tuple(
        Lane(
            f"L{i}",
            qv(0, 0, 0),
            qv(i + 1, (i + 1) ** 2 + 1, (i + 1) ** 3 - 2 * i),
        )
        for i in range(n)
    )
    events = tuple(
        make_event(f"E{i}_{j}", i, j, Q(0), lanes)
        for i, j in itertools.combinations(range(n), 2)
    )
    return Family(
        "simultaneous_bush_K5",
        "single spacetime bush",
        lanes,
        events,
        "All pair events are local but a binary hard-sphere movie may not reuse one particle at the same instant.",
    )


def planar_weave(n: int = 6) -> Family:
    """A dense one-spatial-dimensional line arrangement with distinct crossings."""

    if n != 6:
        raise ValueError("this audited weave uses n=6")
    left = [-15, 1, 4, 5, 16, 29]
    right = [19, 18, 7, -10, -16, -18]
    lanes = tuple(
        Lane(f"L{i}", qv(left[i], 0, 0), qv(right[i] - left[i], 0, 0))
        for i in range(n)
    )
    events: list[Event] = []
    for i, j in itertools.combinations(range(n), 2):
        vi = lanes[i].velocity[0]
        vj = lanes[j].velocity[0]
        t = Q(left[j] - left[i], 1) / (vi - vj)
        events.append(make_event(f"E{i}_{j}", i, j, t, lanes))
    return Family(
        "planar_weave_K6",
        "coherent plany carrier",
        lanes,
        tuple(sorted(events, key=lambda e: (e.time, e.name))),
        "All 15 crossings are distinct and lie in 0<t<1; this is the hard-rod escape branch.",
    )


def regulus_grid(p: int = 4, q: int = 4, shear: bool = False) -> Family:
    """Two rulings of a rational quadric; every A_s meets every B_u."""

    ss = [Q(1, 20), Q(7, 50), Q(31, 100), Q(13, 20)]
    us = [Q(1, 50), Q(11, 100), Q(2, 5), Q(17, 20)]
    if p != 4 or q != 4:
        raise ValueError("this audited regulus uses 4+4 lanes")

    def transform(v: tuple[Q, Q, Q]) -> tuple[Q, Q, Q]:
        if not shear:
            return v
        x, y, z = v
        return (x + Q(2, 3) * y + Q(1, 5) * z, y + Q(1, 4) * z, z)

    lanes: list[Lane] = []
    for i, s in enumerate(ss):
        lanes.append(Lane(f"A{i}", transform(qv(-s, -s * s, 0)), transform(qv(1, 2 * s, 0))))
    for j, u in enumerate(us):
        lanes.append(Lane(f"B{j}", transform(qv(u, -u * u, 0)), transform(qv(-1, 2 * u, 0))))
    lane_tuple = tuple(lanes)
    events: list[Event] = []
    for i, s in enumerate(ss):
        for j, u in enumerate(us):
            t = (s + u) / 2
            events.append(make_event(f"E{i}_{j}", i, p + j, t, lane_tuple))
    tag = "sheared_" if shear else ""
    return Family(
        f"{tag}regulus_K4_4",
        "ruled quadric carrier",
        lane_tuple,
        tuple(sorted(events, key=lambda e: (e.time, e.name))),
        "A determinant-one rational spatial shear is applied." if shear else "Exact rational doubly ruled grid.",
    )


def regulus_scaling_family(n: int) -> Family:
    """A deterministic K_{n,n} regulus with pairwise-distinct event labels."""

    denominator = 2 * n * (n + 1) + 1
    ss = [Q(i + 1, denominator) for i in range(n)]
    us = [Q(n * (j + 1), denominator) for j in range(n)]
    lanes: list[Lane] = []
    for i, s in enumerate(ss):
        lanes.append(Lane(f"A{i}", qv(-s, -s * s, 0), qv(1, 2 * s, 0)))
    for j, u in enumerate(us):
        lanes.append(Lane(f"B{j}", qv(u, -u * u, 0), qv(-1, 2 * u, 0)))
    lane_tuple = tuple(lanes)
    events: list[Event] = []
    for i, s in enumerate(ss):
        for j, u in enumerate(us):
            events.append(make_event(f"E{i}_{j}", i, n + j, (s + u) / 2, lane_tuple))
    return Family(
        f"regulus_scaling_K{n}_{n}",
        "ruled quadric carrier",
        lane_tuple,
        tuple(sorted(events, key=lambda e: (e.time, e.name))),
        "Used only for the deterministic scaling probe; no maximum-subset claim is made for n>4.",
    )


def event_grid_indices(event: Event) -> tuple[int, int]:
    left, right = event.name.split("_")
    return int(left[1:]), int(right)


def regulus_scaling_probe() -> list[dict]:
    results: list[dict] = []
    for n in range(2, 9):
        family = regulus_scaling_family(n)
        full_selection = tuple(range(len(family.events)))
        full_ok, full_audit = subset_wireable(family, full_selection, degree_cap=None)
        skeletons = []
        for degree in range(2, min(n, 5) + 1):
            allowed_offsets = set(range(degree))
            selected = tuple(
                event_index
                for event_index, event in enumerate(family.events)
                if (event_grid_indices(event)[1] - event_grid_indices(event)[0]) % n
                in allowed_offsets
            )
            ok, audit = subset_wireable(family, selected, degree_cap=None)
            skeletons.append(
                {
                    "circulant_degree": degree,
                    "selected_events": len(selected),
                    "retention_lower_bound_if_wireable": len(selected) / len(family.events),
                    "wireable": ok,
                    "matching_audit": audit,
                }
            )
        results.append(
            {
                "n_per_ruling": n,
                "lane_count": len(family.lanes),
                "candidate_event_count": len(family.events),
                "full_grid_wireable": full_ok,
                "full_grid_audit": full_audit,
                "tested_circulant_skeletons": skeletons,
                "warning": "Wireable skeleton fractions are certified lower bounds, not optimal retention for n>4.",
            }
        )
    return results


def incident_sign(event: Event, lane: int, incoming: bool) -> int:
    if lane == event.lane_i:
        return +1 if incoming else -1
    if lane == event.lane_j:
        return -1 if incoming else +1
    raise ValueError("lane is not incident to event")


def selected_degrees(family: Family, selected: Sequence[int]) -> list[int]:
    degrees = [0] * len(family.lanes)
    for event_index in selected:
        event = family.events[event_index]
        degrees[event.lane_i] += 1
        degrees[event.lane_j] += 1
    return degrees


def chronology_valid(family: Family, selected: Sequence[int]) -> bool:
    by_lane: list[list[Q]] = [[] for _ in family.lanes]
    for event_index in selected:
        event = family.events[event_index]
        by_lane[event.lane_i].append(event.time)
        by_lane[event.lane_j].append(event.time)
    return all(len(times) == len(set(times)) for times in by_lane)


def token_identity_audit(family: Family, selected: Sequence[int]) -> dict:
    """Trace particle labels through central swaps in chronological order."""

    if not chronology_valid(family, selected):
        return {
            "valid": False,
            "reason": "one lane is asked to undergo two binary collisions at the same time",
            "extra_global_identity_entropy_bits": None,
        }
    lane_token = list(range(len(family.lanes)))
    ordered = sorted(selected, key=lambda e: (family.events[e].time, family.events[e].name))
    for event_index in ordered:
        event = family.events[event_index]
        lane_token[event.lane_i], lane_token[event.lane_j] = (
            lane_token[event.lane_j],
            lane_token[event.lane_i],
        )
    return {
        "valid": sorted(lane_token) == list(range(len(family.lanes))),
        "chronological_event_names": [family.events[e].name for e in ordered],
        "final_token_on_lane": lane_token,
        "local_route_count_per_event": 1,
        "global_route_count": 1,
        "extra_global_identity_entropy_bits": 0,
        "interpretation": "For central swaps, labels are a deterministic permutation; identity adds no nonlocal constraint.",
    }


def build_matching_system(family: Family, selected: Sequence[int]) -> tuple[np.ndarray, np.ndarray, list[float]]:
    selected = tuple(selected)
    local_index = {event_index: k for k, event_index in enumerate(selected)}
    m = len(selected)
    omegas: list[np.ndarray] = []
    for event_index in selected:
        event = family.events[event_index]
        vi = as_float(family.lanes[event.lane_i].velocity)
        vj = as_float(family.lanes[event.lane_j].velocity)
        relative = vi - vj
        norm = float(np.linalg.norm(relative))
        if norm <= 1.0e-14:
            raise ValueError("parallel equal-velocity lanes cannot make a central collision")
        omegas.append(-relative / norm)

    lane_events: list[list[int]] = [[] for _ in family.lanes]
    for event_index in selected:
        event = family.events[event_index]
        lane_events[event.lane_i].append(event_index)
        lane_events[event.lane_j].append(event_index)
    for lane_index in range(len(lane_events)):
        lane_events[lane_index].sort(key=lambda e: (family.events[e].time, family.events[e].name))

    row_count = 3 * sum(max(0, len(xs) - 1) for xs in lane_events)
    a = np.zeros((row_count, 4 * m), dtype=float)
    rhs = np.zeros(row_count, dtype=float)
    row = 0
    time_gaps: list[float] = []
    for lane_index, event_indices in enumerate(lane_events):
        velocity = as_float(family.lanes[lane_index].velocity)
        for early_index, late_index in zip(event_indices, event_indices[1:]):
            early = family.events[early_index]
            late = family.events[late_index]
            if late.time <= early.time:
                raise ValueError("selected events are not strictly ordered on a lane")
            time_gaps.append(float(late.time - early.time))
            early_local = local_index[early_index]
            late_local = local_index[late_index]
            omega_early = omegas[early_local]
            omega_late = omegas[late_local]
            out_sign = incident_sign(early, lane_index, incoming=False)
            in_sign = incident_sign(late, lane_index, incoming=True)
            local_rhs = -0.5 * (in_sign * omega_late - out_sign * omega_early)
            for coord in range(3):
                a[row, 4 * late_local + coord] += 1.0
                a[row, 4 * early_local + coord] -= 1.0
                a[row, 4 * late_local + 3] -= velocity[coord]
                a[row, 4 * early_local + 3] += velocity[coord]
                rhs[row] = local_rhs[coord]
                row += 1
    return a, rhs, time_gaps


def matrix_audit(a: np.ndarray, rhs: np.ndarray) -> dict[str, float | int | bool]:
    rows, cols = a.shape
    if rows == 0:
        return {
            "rows": 0,
            "cols": int(cols),
            "rank": 0,
            "left_null_dimension": 0,
            "relative_residual": 0.0,
            "absolute_residual": 0.0,
            "solution_norm": 0.0,
            "smallest_nonzero_singular_value": 0.0,
            "condition_number_nonzero": 1.0,
            "wireable": True,
        }
    solution, _, rank, singular_values = np.linalg.lstsq(a, rhs, rcond=1.0e-12)
    residual_vector = a @ solution - rhs
    absolute_residual = float(np.linalg.norm(residual_vector))
    rhs_norm = float(np.linalg.norm(rhs))
    relative_residual = absolute_residual / max(rhs_norm, 1.0)
    nonzero = singular_values[singular_values > max(a.shape) * np.finfo(float).eps * singular_values[0]]
    smallest = float(nonzero[-1]) if len(nonzero) else 0.0
    condition = float(nonzero[0] / nonzero[-1]) if len(nonzero) else 1.0
    return {
        "rows": int(rows),
        "cols": int(cols),
        "rank": int(rank),
        "left_null_dimension": int(rows - rank),
        "relative_residual": relative_residual,
        "absolute_residual": absolute_residual,
        "solution_norm": float(np.linalg.norm(solution)),
        "smallest_nonzero_singular_value": smallest,
        "condition_number_nonzero": condition,
        "wireable": bool(relative_residual <= 2.0e-9),
    }


def subset_wireable(family: Family, selected: Sequence[int], degree_cap: int | None) -> tuple[bool, dict]:
    if not chronology_valid(family, selected):
        return False, {"reason": "simultaneous reuse"}
    degrees = selected_degrees(family, selected)
    if degree_cap is not None and max(degrees, default=0) > degree_cap:
        return False, {"reason": "degree cap"}
    a, rhs, gaps = build_matching_system(family, selected)
    audit = matrix_audit(a, rhs)
    audit["minimum_time_gap"] = min(gaps) if gaps else None
    return bool(audit["wireable"]), audit


def maximum_wireable_subset(family: Family, degree_cap: int | None) -> dict:
    m = len(family.events)
    checked = 0
    for size in range(m, -1, -1):
        witnesses: list[tuple[int, ...]] = []
        witness_audit: dict | None = None
        for subset in itertools.combinations(range(m), size):
            checked += 1
            feasible, audit = subset_wireable(family, subset, degree_cap)
            if feasible:
                witnesses.append(subset)
                if witness_audit is None:
                    witness_audit = audit
                if len(witnesses) >= 8:
                    break
        if witnesses:
            chosen = witnesses[0]
            return {
                "maximum_events": size,
                "event_retention": size / m if m else 1.0,
                "incidence_retention": size / m if m else 1.0,
                "checked_subsets_until_first_optimum": checked,
                "witness_event_indices": list(chosen),
                "witness_event_names": [family.events[i].name for i in chosen],
                "witness_audit": witness_audit,
                "at_least_this_many_optimal_witnesses_found": len(witnesses),
            }
    raise AssertionError("empty subset should always be feasible")


def exact_geometry_audit(family: Family) -> dict:
    exact_residuals = []
    repeated_pair = False
    seen_pairs: set[tuple[int, int]] = set()
    for event in family.events:
        pair = tuple(sorted((event.lane_i, event.lane_j)))
        repeated_pair = repeated_pair or pair in seen_pairs
        seen_pairs.add(pair)
        pi = qposition(family.lanes[event.lane_i], event.time)
        pj = qposition(family.lanes[event.lane_j], event.time)
        exact_residuals.append(pi == event.point and pj == event.point)

    velocities = [lane.velocity for lane in family.lanes]
    distinct_directions = len(set(velocities)) == len(velocities)
    all_times_in_unit_window = all(Q(0) <= event.time <= Q(1) for event in family.events)
    time_multiplicities: dict[Q, int] = {}
    for event in family.events:
        time_multiplicities[event.time] = time_multiplicities.get(event.time, 0) + 1

    declared_pairs = {tuple(sorted((event.lane_i, event.lane_j))) for event in family.events}
    unintended_exact_pairs = []
    for i, j in itertools.combinations(range(len(family.lanes)), 2):
        if (i, j) in declared_pairs:
            continue
        db = qsub(family.lanes[i].root, family.lanes[j].root)
        dv = qsub(family.lanes[i].velocity, family.lanes[j].velocity)
        candidate_times = [-db[k] / dv[k] for k in range(3) if dv[k] != 0]
        if not candidate_times:
            intersects = all(x == 0 for x in db)
            meeting_time = None
        else:
            meeting_time = candidate_times[0]
            intersects = all(db[k] + meeting_time * dv[k] == 0 for k in range(3))
        if intersects:
            unintended_exact_pairs.append(
                {
                    "lanes": [family.lanes[i].name, family.lanes[j].name],
                    "time": str(meeting_time) if meeting_time is not None else "coincident lines",
                }
            )
    return {
        "all_declared_intersections_exact_over_Q": all(exact_residuals),
        "distinct_velocity_vectors": distinct_directions,
        "no_repeated_lane_pair": not repeated_pair,
        "all_event_times_in_0_1": all_times_in_unit_window,
        "distinct_event_time_count": len(time_multiplicities),
        "maximum_event_time_multiplicity": max(time_multiplicities.values(), default=0),
        "unintended_exact_line_intersection_count": len(unintended_exact_pairs),
        "unintended_exact_line_intersections": unintended_exact_pairs,
    }


def digest_family(family: Family) -> str:
    payload = {
        "name": family.name,
        "lanes": [
            {"name": lane.name, "root": list(map(str, lane.root)), "velocity": list(map(str, lane.velocity))}
            for lane in family.lanes
        ],
        "events": [
            {
                "name": event.name,
                "lanes": [event.lane_i, event.lane_j],
                "time": str(event.time),
                "point": list(map(str, event.point)),
            }
            for event in family.events
        ],
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def audit_family(family: Family) -> dict:
    all_events = tuple(range(len(family.events)))
    chronology = chronology_valid(family, all_events)
    if chronology:
        full_a, full_rhs, full_gaps = build_matching_system(family, all_events)
        full = matrix_audit(full_a, full_rhs)
        full["minimum_time_gap"] = min(full_gaps) if full_gaps else None
    else:
        full = {
            "wireable": False,
            "reason": "one or more lanes are reused at the same time",
        }

    unconstrained = maximum_wireable_subset(family, degree_cap=None)
    cap_two = maximum_wireable_subset(family, degree_cap=2)
    cap_four = maximum_wireable_subset(family, degree_cap=4)
    return {
        "name": family.name,
        "geometry_class": family.geometry_class,
        "notes": family.notes,
        "lane_count": len(family.lanes),
        "candidate_event_count": len(family.events),
        "candidate_incidence_count": 2 * len(family.events),
        "family_sha256": digest_family(family),
        "geometry_audit": exact_geometry_audit(family),
        "discrete_identity_audit": token_identity_audit(family, all_events),
        "full_network": full,
        "maximum_compilable": {
            "no_lane_reuse_cap": unconstrained,
            "at_most_two_events_per_original_lane": cap_two,
            "at_most_four_events_per_original_lane": cap_four,
        },
    }


def build_certificate() -> dict:
    families = (
        strict_c12_family(),
        generic_polygon_cycle(),
        time_hinge_star(),
        simultaneous_bush(),
        planar_weave(),
        regulus_grid(shear=False),
        regulus_grid(shear=True),
    )
    results = [audit_family(family) for family in families]

    # Key separation: lane continuity fixes the coefficient matrix A.  The
    # hard-sphere law changes only c(omega).  Thus rank loss is not a new
    # collision entropy; only a nonzero left-null projection of c is.
    collision_specific_obstructions = []
    for result in results:
        full = result["full_network"]
        collision_specific_obstructions.append(
            {
                "family": result["name"],
                "chronology_valid": "reason" not in full,
                "left_null_dimension": full.get("left_null_dimension"),
                "relative_obstruction": full.get("relative_residual"),
                "collision_specific_obstruction_detected": bool(
                    full.get("relative_residual", 0.0) > 2.0e-9
                ),
            }
        )

    return {
        "metadata": {
            "title": "Global legal collision-movie compilation audit",
            "date": "2026-08-29",
            "arithmetic": "exact rational incidence geometry; double-precision SVD/lstsq for unit-normal wiring",
            "scope": "finite central-swap hard-sphere compiler; not an asymptotic Kakeya extraction theorem",
            "randomness": "none",
        },
        "model": {
            "worldline": "x_i(t)=b_i+t v_i in R^3, hence an affine line in R^4",
            "local_event": "exact pair concurrence followed by a central equal-mass velocity swap",
            "global_constraint": "finite-radius contact half-rays must match on every shared original lane",
            "matching_system": "A(y_e,s_e)=c(omega_e)",
            "retention": "maximum selected candidate events divided by all candidate events; event and incidence fractions coincide here",
            "collision_specific_entropy_proxy": "dimension and norm of the projection of c(omega) onto the left nullspace of A",
            "congestion_sensitivity": "reported with no cap, degree cap 4, and degree cap 2 per original lane",
        },
        "families": results,
        "regulus_scaling_probe": regulus_scaling_probe(),
        "collision_specific_obstruction_summary": collision_specific_obstructions,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    payload = json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
