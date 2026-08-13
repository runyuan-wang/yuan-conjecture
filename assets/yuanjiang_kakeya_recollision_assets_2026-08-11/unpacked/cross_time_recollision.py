#!/usr/bin/env python3
"""Numerical stress test for a cross-time recollision rank of line families in R^4.

Each spacetime line is represented as x(t) = a + t v, with x,v in R^3.
For every pair, we find its closest encounter on t in [0,1].  Encounters below
the spatial tolerance are colored by a spacetime cell.  The statistic

    beta(global encounter graph) - sum beta(each spacetime-cell graph)

counts independent cycles that require encounters in more than one local
spacetime event.  Thus a single bush/star should score zero, while genuinely
cross-time incidence structures should score positively.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict
import argparse
import math
import numpy as np


class UnionFind:
    def __init__(self, n: int):
        self.parent = np.arange(n, dtype=np.int64)
        self.size = np.ones(n, dtype=np.int64)
        self.components = n

    def find(self, x: int) -> int:
        parent = self.parent
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = int(parent[x])
        return x

    def union(self, x: int, y: int) -> bool:
        rx = self.find(x)
        ry = self.find(y)
        if rx == ry:
            return False
        if self.size[rx] < self.size[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        self.size[rx] += self.size[ry]
        self.components -= 1
        return True


@dataclass
class Family:
    name: str
    a: np.ndarray
    v: np.ndarray

    @property
    def n(self) -> int:
        return int(self.a.shape[0])


@dataclass
class EncounterStats:
    n: int
    edges: int
    event_cells: int
    global_beta: int
    local_beta: int
    cross_beta: int
    global_components: int
    time_span_cells: int

    @property
    def cross_per_line(self) -> float:
        return self.cross_beta / self.n

    @property
    def local_fraction(self) -> float:
        if self.global_beta == 0:
            return 0.0
        return self.local_beta / self.global_beta


def direction_grid(m: int) -> np.ndarray:
    q = np.linspace(-1.0, 1.0, m)
    return np.stack(np.meshgrid(q, q, q, indexing="ij"), axis=-1).reshape(-1, 3)


def make_star(m: int) -> Family:
    v = direction_grid(m)
    t0 = 0.5
    center = np.array([0.17, -0.11, 0.07])
    a = center[None, :] - t0 * v
    return Family(f"star(m={m})", a, v)


def make_worst_affine(m: int) -> Family:
    v = direction_grid(m)
    roots = np.array([
        0.5 - 0.5 / math.sqrt(2.0),
        0.5,
        0.5 + 0.5 / math.sqrt(2.0),
    ])
    a = -v * roots[None, :]
    return Family(f"three_affine_collapses(m={m})", a, v)


def make_random(m: int, seed: int = 7) -> Family:
    rng = np.random.default_rng(seed)
    v = direction_grid(m)
    a = rng.uniform(-0.8, 0.8, size=v.shape)
    return Family(f"random_selector(m={m})", a, v)


def make_smooth_fold(m: int) -> Family:
    v = direction_grid(m)
    a = np.empty_like(v)
    a[:, 0] = -0.42 * v[:, 0] ** 2 + 0.08 * np.sin(math.pi * v[:, 1])
    a[:, 1] = -0.55 * v[:, 1] ** 2 + 0.07 * np.sin(math.pi * v[:, 2])
    a[:, 2] = -0.68 * v[:, 2] ** 2 + 0.06 * np.sin(math.pi * v[:, 0])
    return Family(f"smooth_fold(m={m})", a, v)


def make_two_bush_grid(r: int, seed: int = 11) -> Family:
    """All r^2 lines connect one of r points at t1 to one of r points at t2."""
    rng = np.random.default_rng(seed)
    t1, t2 = 0.27, 0.73
    p = rng.uniform(-0.8, 0.8, size=(r, 3))
    q = rng.uniform(-0.8, 0.8, size=(r, 3))
    aa = []
    vv = []
    for i in range(r):
        for j in range(r):
            vel = (q[j] - p[i]) / (t2 - t1)
            base = p[i] - t1 * vel
            aa.append(base)
            vv.append(vel)
    return Family(f"two_bush_grid(r={r})", np.asarray(aa), np.asarray(vv))


def make_regulus(m: int) -> Family:
    """Two rulings of z=xy after choosing a generic time coordinate.

    Ruling A_s: a=(-s,-s^2,0), v=(1,2s,0)
    Ruling B_s: a=( s,-s^2,0), v=(-1,2s,0)
    Every A_s and B_u meet at t=(s+u)/2.
    """
    s = np.linspace(0.05, 0.95, m)
    a1 = np.column_stack((-s, -s**2, np.zeros_like(s)))
    v1 = np.column_stack((np.ones_like(s), 2 * s, np.zeros_like(s)))
    a2 = np.column_stack((s, -s**2, np.zeros_like(s)))
    v2 = np.column_stack((-np.ones_like(s), 2 * s, np.zeros_like(s)))
    return Family(f"regulus({m}+{m})", np.vstack((a1, a2)), np.vstack((v1, v2)))


def make_book(pages: int, per_page: int) -> Family:
    """Rotating pages through the central time-axis, plus that central line."""
    aa = [np.zeros(3)]
    vv = [np.zeros(3)]
    taus = np.linspace(0.12, 0.88, per_page)
    lambdas = np.linspace(0.55, 1.55, per_page)
    for p in range(pages):
        phi = 2.0 * math.pi * p / pages
        # Non-coplanar page normals in spatial R^3.
        w = np.array([math.cos(phi), math.sin(phi), 0.35 * math.sin(2 * phi)])
        w /= np.linalg.norm(w)
        for tau, lam in zip(taus, np.roll(lambdas, 2 * p % per_page)):
            vel = lam * w
            aa.append(-tau * vel)
            vv.append(vel)
    return Family(
        f"rotating_book(pages={pages},per={per_page})",
        np.asarray(aa),
        np.asarray(vv),
    )


def graph_beta(n: int, edges: list[tuple[int, int]]) -> tuple[int, int]:
    uf = UnionFind(n)
    for i, j in edges:
        uf.union(i, j)
    beta = len(edges) - n + uf.components
    return int(beta), int(uf.components)


def encounter_stats(
    family: Family,
    tolerance: float,
    time_cell: float,
    space_cell: float,
    min_angle: float = 0.0,
) -> EncounterStats:
    a = family.a
    v = family.v
    n = family.n
    tol2 = tolerance * tolerance
    edges: list[tuple[int, int]] = []
    event_edges: dict[tuple[int, int, int, int], list[tuple[int, int]]] = defaultdict(list)
    time_bins = set()

    for i in range(n - 1):
        da = a[i + 1 :] - a[i]
        dv = v[i + 1 :] - v[i]
        speed2 = np.einsum("ij,ij->i", dv, dv)
        angle_ok = speed2 >= min_angle * min_angle
        movable = speed2 > 1e-24
        raw_t = np.zeros_like(speed2)
        raw_t[movable] = -np.einsum("ij,ij->i", da[movable], dv[movable]) / speed2[movable]
        t = np.clip(raw_t, 0.0, 1.0)
        sep = da + t[:, None] * dv
        d2 = np.einsum("ij,ij->i", sep, sep)
        mask = (d2 <= tol2) & angle_ok
        js = np.flatnonzero(mask) + i + 1
        for local_idx, j in zip(np.flatnonzero(mask), js):
            tt = float(t[local_idx])
            xi = a[i] + tt * v[i]
            xj = a[j] + tt * v[j]
            midpoint = 0.5 * (xi + xj)
            tb = int(math.floor(tt / time_cell + 1e-12))
            sb = tuple(np.floor(midpoint / space_cell + 1e-12).astype(int).tolist())
            color = (tb, sb[0], sb[1], sb[2])
            edge = (i, int(j))
            edges.append(edge)
            event_edges[color].append(edge)
            time_bins.add(tb)

    global_beta, components = graph_beta(n, edges)
    local_beta = 0
    for group in event_edges.values():
        beta, _ = graph_beta(n, group)
        local_beta += beta
    cross_beta = global_beta - local_beta
    if cross_beta < 0:
        raise RuntimeError("Cross beta became negative; event edge partition is inconsistent")
    return EncounterStats(
        n=n,
        edges=len(edges),
        event_cells=len(event_edges),
        global_beta=global_beta,
        local_beta=local_beta,
        cross_beta=cross_beta,
        global_components=components,
        time_span_cells=len(time_bins),
    )


def run_suite(delta: float, grid_m: int) -> None:
    # Event cells are intentionally a few tube radii wide.  We also discard
    # near-parallel pairs here; those belong to a separate angular scale.
    tolerance = delta
    time_cell = 4.0 * delta
    space_cell = 4.0 * delta
    min_angle = 0.18
    families = [
        make_random(grid_m),
        make_star(grid_m),
        make_worst_affine(grid_m),
        make_smooth_fold(grid_m),
        make_two_bush_grid(16),
        make_book(10, 18),
        make_regulus(55),
    ]
    print(
        f"delta={delta:g}  time_cell={time_cell:g}  "
        f"space_cell={space_cell:g}  min_angle={min_angle:g}"
    )
    print(
        "family                              N       E   events  beta_all  "
        "beta_local  beta_cross  cross/N  local%"
    )
    for fam in families:
        stats = encounter_stats(
            fam,
            tolerance=tolerance,
            time_cell=time_cell,
            space_cell=space_cell,
            min_angle=min_angle,
        )
        print(
            f"{fam.name:34s} {stats.n:5d} {stats.edges:7d} "
            f"{stats.event_cells:7d} {stats.global_beta:9d} "
            f"{stats.local_beta:11d} {stats.cross_beta:11d} "
            f"{stats.cross_per_line:8.3f} {100*stats.local_fraction:7.2f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta", type=float, default=0.01)
    parser.add_argument("--grid-m", type=int, default=7)
    args = parser.parse_args()
    run_suite(args.delta, args.grid_m)


if __name__ == "__main__":
    main()
