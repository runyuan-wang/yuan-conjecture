#!/usr/bin/env python3
"""Entropy/cycle diagnostics for the four-slope arithmetic Kakeya problem.

The input layers are X, Y, X+Y, and X+2Y.  The output X-Y is required to be
injective on the support, so H(X-Y) equals the entropy of the support atom.

For the cycle statistic, make a bipartite incidence graph:

    support atom -- projection value at one input slope.

Each projection fibre is one collision event.  Its circuit rank is a canonical
version of the cross-layer recollision count: cycles contained in a single
fibre have already been collapsed to one event node.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
import json
from dataclasses import dataclass
from math import log2
from typing import Callable, Iterable, Sequence
from urllib.request import urlopen


Point = tuple[int, int, float]
Projection = tuple[str, Callable[[int, int], int]]


INPUT_PROJECTIONS: tuple[Projection, ...] = (
    ("X", lambda x, y: x),
    ("Y", lambda x, y: y),
    ("X+Y", lambda x, y: x + y),
    ("X+2Y", lambda x, y: x + 2 * y),
)
OUTPUT_PROJECTION: Projection = ("X-Y", lambda x, y: x - y)


# Sebastian Griego's May 2026 certificate, as displayed in teorth/
# optimizationproblems PR #70.  The rounded masses are normalized below.
CERTIFICATE_26: tuple[Point, ...] = (
    (1, 7, 1.30506564e-12),
    (3, 5, 9.25745153e-7),
    (3, 6, 1.52498825e-6),
    (3, 7, 6.57289640e-8),
    (5, 4, 0.001244422280),
    (5, 5, 0.006113255152),
    (5, 6, 0.000037354259),
    (7, 2, 0.016223522296),
    (7, 3, 0.169721258229),
    (7, 4, 0.008420798270),
    (7, 5, 0.011728626004),
    (9, 1, 0.248375293086),
    (9, 2, 0.046862624648),
    (9, 3, 0.243496670072),
    (11, 0, 0.010100758400),
    (11, 1, 0.188770387301),
    (11, 2, 0.018228788494),
    (13, -1, 0.006653795672),
    (13, 0, 0.002961930500),
    (13, 1, 0.020758390703),
    (15, -2, 0.000003025952),
    (15, -1, 0.000275897370),
    (15, 0, 0.000020680540),
    (17, -3, 6.99564607e-12),
    (17, -2, 1.51075546e-9),
    (17, -1, 2.78909933e-9),
)


CERTIFICATE_95_URL = (
    "https://gist.githubusercontent.com/463464q435q43/"
    "1b3aa1b58016076808dc336979b3a063/raw/"
    "251682e88835f1af20bc7b5af8e5c219f2596274/"
    "certificate_3c_v3_95pt_20260610.json"
)


@dataclass(frozen=True)
class Diagnostics:
    support_size: int
    entropy_output: float
    input_entropies: dict[str, float]
    ratio: float
    value_counts: dict[str, int]
    conditional_losses: dict[str, float]
    entropy_cycle_surplus: float
    entropy_imbalance_penalty: float
    effective_entropy_surplus: float
    incidence_vertices: int
    incidence_edges: int
    components: int
    circuit_rank: int
    girth: int | None
    six_cycles: int
    six_cycle_color_types: dict[tuple[str, ...], int]
    output_injective: bool


def normalize(points: Iterable[Point]) -> list[Point]:
    kept = [(x, y, float(p)) for x, y, p in points if p > 0]
    total = sum(p for _, _, p in kept)
    if total <= 0:
        raise ValueError("positive total mass required")
    return [(x, y, p / total) for x, y, p in kept]


def load_certificate_95() -> list[Point]:
    """Load Mosaic Intelligence's public June 2026 exact certificate."""
    with urlopen(CERTIFICATE_95_URL, timeout=30) as response:
        payload = json.load(response)
    result: list[Point] = []
    for x, y, (numerator, denominator) in payload["distribution_x_y_numden"]:
        result.append((int(x), int(y), int(numerator) / int(denominator)))
    if len(result) != 95:
        raise ValueError(f"expected 95 support points, got {len(result)}")
    return result


def entropy(masses: Iterable[float]) -> float:
    return -sum(p * log2(p) for p in masses if p > 0)


def pushforward(points: Sequence[Point], fn: Callable[[int, int], int]) -> dict[int, float]:
    result: dict[int, float] = defaultdict(float)
    for x, y, p in points:
        result[fn(x, y)] += p
    return dict(result)


def incidence_graph(
    points: Sequence[Point],
) -> tuple[dict[str, set[str]], dict[str, int]]:
    graph: dict[str, set[str]] = defaultdict(set)
    value_counts: dict[str, int] = {}
    for label, fn in INPUT_PROJECTIONS:
        values = {fn(x, y) for x, y, _ in points}
        value_counts[label] = len(values)
        for i, (x, y, _) in enumerate(points):
            atom = f"a:{i}"
            event = f"e:{label}:{fn(x, y)}"
            graph[atom].add(event)
            graph[event].add(atom)
    return dict(graph), value_counts


def component_count(graph: dict[str, set[str]]) -> int:
    unseen = set(graph)
    count = 0
    while unseen:
        count += 1
        start = unseen.pop()
        stack = [start]
        while stack:
            node = stack.pop()
            for neighbor in graph[node]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
    return count


def graph_girth(graph: dict[str, set[str]]) -> int | None:
    best: int | None = None
    for start in graph:
        dist = {start: 0}
        parent = {start: None}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            if best is not None and 2 * dist[node] + 1 >= best:
                continue
            for neighbor in graph[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    queue.append(neighbor)
                elif parent[node] != neighbor:
                    length = dist[node] + dist[neighbor] + 1
                    if best is None or length < best:
                        best = length
    return best


def collision_edges(points: Sequence[Point]) -> dict[tuple[int, int], str]:
    """Return the input-layer color of every colliding atom pair."""
    edges: dict[tuple[int, int], str] = {}
    for label, fn in INPUT_PROJECTIONS:
        fibres: dict[int, list[int]] = defaultdict(list)
        for i, (x, y, _) in enumerate(points):
            fibres[fn(x, y)].append(i)
        for fibre in fibres.values():
            for offset, i in enumerate(fibre):
                for j in fibre[offset + 1 :]:
                    pair = (min(i, j), max(i, j))
                    if pair in edges:
                        raise ValueError("two distinct projections identify the same atom pair")
                    edges[pair] = label
    return edges


def six_cycle_types(points: Sequence[Point]) -> Counter[tuple[str, ...]]:
    """Count length-six cycles in the atom--projection-value incidence graph.

    Such a cycle is equivalently a triangle of atoms whose three pairwise
    collisions occur in three different input layers.
    """
    edges = collision_edges(points)
    result: Counter[tuple[str, ...]] = Counter()
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            label_ij = edges.get((i, j))
            if label_ij is None:
                continue
            for k in range(j + 1, n):
                label_ik = edges.get((i, k))
                label_jk = edges.get((j, k))
                if label_ik is None or label_jk is None:
                    continue
                labels = (label_ij, label_ik, label_jk)
                if len(set(labels)) == 3:
                    result[tuple(sorted(labels))] += 1
    return result


def six_cycle_participation(points: Sequence[Point]) -> list[int]:
    edges = collision_edges(points)
    counts = [0] * len(points)
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            label_ij = edges.get((i, j))
            if label_ij is None:
                continue
            for k in range(j + 1, len(points)):
                labels = (label_ij, edges.get((i, k)), edges.get((j, k)))
                if None not in labels and len(set(labels)) == 3:
                    counts[i] += 1
                    counts[j] += 1
                    counts[k] += 1
    return counts


def diagnose(points: Iterable[Point]) -> Diagnostics:
    pts = normalize(points)
    output_label, output_fn = OUTPUT_PROJECTION
    output_dist = pushforward(pts, output_fn)
    output_entropy = entropy(output_dist.values())
    input_entropies = {
        label: entropy(pushforward(pts, fn).values())
        for label, fn in INPUT_PROJECTIONS
    }
    denominator = max(input_entropies.values())
    entropy_sum = sum(input_entropies.values())
    entropy_cycle_surplus = 3 * output_entropy - entropy_sum
    imbalance_penalty = 4 * denominator - entropy_sum
    graph, value_counts = incidence_graph(pts)
    cycle_types = six_cycle_types(pts)
    vertices = len(graph)
    edges = sum(map(len, graph.values())) // 2
    components = component_count(graph)
    rank = edges - vertices + components
    return Diagnostics(
        support_size=len(pts),
        entropy_output=output_entropy,
        input_entropies=input_entropies,
        ratio=output_entropy / denominator,
        value_counts=value_counts,
        conditional_losses={
            label: output_entropy - h for label, h in input_entropies.items()
        },
        entropy_cycle_surplus=entropy_cycle_surplus,
        entropy_imbalance_penalty=imbalance_penalty,
        effective_entropy_surplus=entropy_cycle_surplus - imbalance_penalty,
        incidence_vertices=vertices,
        incidence_edges=edges,
        components=components,
        circuit_rank=rank,
        girth=graph_girth(graph),
        six_cycles=sum(cycle_types.values()),
        six_cycle_color_types=dict(cycle_types),
        output_injective=len(output_dist) == len(pts),
    )


def collision_pair_colors(points: Iterable[Point]) -> Counter[int]:
    """How many input layers identify each unordered atom pair.

    Distinct linear projections should never identify the same two distinct
    lattice points twice.  A value above one would signal a coding error or
    duplicate support points.
    """
    pts = normalize(points)
    colors: Counter[tuple[int, int]] = Counter()
    for _, fn in INPUT_PROJECTIONS:
        fibres: dict[int, list[int]] = defaultdict(list)
        for i, (x, y, _) in enumerate(pts):
            fibres[fn(x, y)].append(i)
        for fibre in fibres.values():
            for offset, i in enumerate(fibre):
                for j in fibre[offset + 1 :]:
                    colors[(i, j)] += 1
    return Counter(colors.values())


def pruned(points: Sequence[Point], threshold: float) -> list[Point]:
    return [(x, y, p) for x, y, p in points if p >= threshold]


def star_baseline(n: int = 26) -> list[Point]:
    # One layer (Y) collapses completely; the other three layers are injective.
    return [(i, 0, 1.0 / n) for i in range(n)]


def uniformized(points: Sequence[Point]) -> list[Point]:
    mass = 1.0 / len(points)
    return [(x, y, mass) for x, y, _ in points]


def print_diagnostics(name: str, result: Diagnostics) -> None:
    counts = ", ".join(f"{k}:{v}" for k, v in result.value_counts.items())
    entropies = ", ".join(f"{k}:{v:.12f}" for k, v in result.input_entropies.items())
    losses = ", ".join(f"{k}:{v:.12f}" for k, v in result.conditional_losses.items())
    print(name)
    print(
        f"  N={result.support_size} output_H={result.entropy_output:.12f} "
        f"ratio={result.ratio:.12f} injective={result.output_injective}"
    )
    print(f"  input_H=[{entropies}]")
    print(f"  losses=[{losses}]")
    print(
        "  entropy_loop_decomposition "
        f"raw={result.entropy_cycle_surplus:.12f} "
        f"imbalance={result.entropy_imbalance_penalty:.12f} "
        f"effective={result.effective_entropy_surplus:.12f} "
        f"effective/H={result.effective_entropy_surplus / result.entropy_output:.12f}"
    )
    print(f"  projection_values=[{counts}]")
    print(
        f"  incidence V={result.incidence_vertices} E={result.incidence_edges} "
        f"C={result.components} rho={result.circuit_rank} girth={result.girth} "
        f"six_cycles={result.six_cycles}"
    )
    if result.six_cycle_color_types:
        print(f"  six_cycle_types={result.six_cycle_color_types}")


def print_removal_analysis(points: Sequence[Point], limit: int = 12) -> None:
    pts = normalize(points)
    base = diagnose(pts)
    participation = six_cycle_participation(pts)
    rows = []
    for i, (x, y, p) in enumerate(pts):
        reduced = pts[:i] + pts[i + 1 :]
        result = diagnose(reduced)
        rows.append(
            {
                "x": x,
                "y": y,
                "p": p,
                "cycles": participation[i],
                "weighted_cycles": p * participation[i],
                "delta_rho": base.circuit_rank - result.circuit_rank,
                "ratio_loss": base.ratio - result.ratio,
            }
        )

    print("single-atom deletion diagnostic (weights not re-optimized)")
    print("  sorted by ratio loss")
    for row in sorted(rows, key=lambda item: item["ratio_loss"], reverse=True)[:limit]:
        print(
            "  "
            f"({row['x']:>3},{row['y']:>3}) p={row['p']:.6g} "
            f"cycles={row['cycles']:>3} p*cycles={row['weighted_cycles']:.6g} "
            f"delta_rho={row['delta_rho']:>2} ratio_loss={row['ratio_loss']:+.9g}"
        )

    try:
        from scipy.stats import spearmanr

        target = [row["ratio_loss"] for row in rows]
        print("  Spearman correlation with ratio loss")
        for key in ("p", "cycles", "weighted_cycles", "delta_rho"):
            statistic = float(spearmanr([row[key] for row in rows], target).statistic)
            print(f"    {key:>15}: {statistic:+.6f}")
    except (ImportError, ValueError):
        pass


def main() -> None:
    print_diagnostics("single-star baseline", diagnose(star_baseline()))
    print(f"  pair collision color multiplicities={dict(collision_pair_colors(star_baseline()))}")
    print()

    print_diagnostics("26-point certificate", diagnose(CERTIFICATE_26))
    print(f"  pair collision color multiplicities={dict(collision_pair_colors(CERTIFICATE_26))}")
    print_removal_analysis(CERTIFICATE_26)
    print()

    print_diagnostics(
        "same 26-point support with uniform weights",
        diagnose(uniformized(CERTIFICATE_26)),
    )
    print()

    for threshold in (1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 1e-2):
        subset = pruned(CERTIFICATE_26, threshold)
        print_diagnostics(f"certificate pruned at p >= {threshold:g}", diagnose(subset))
        print()

    certificate_95 = load_certificate_95()
    print_diagnostics("95-point June 2026 record certificate", diagnose(certificate_95))
    print()
    print_diagnostics(
        "same 95-point support with uniform weights",
        diagnose(uniformized(certificate_95)),
    )
    print()
    for threshold in (1e-100, 1e-50, 1e-20, 1e-12, 1e-9, 1e-6, 1e-3):
        subset = pruned(certificate_95, threshold)
        print_diagnostics(f"95-point certificate pruned at p >= {threshold:g}", diagnose(subset))
        print()


if __name__ == "__main__":
    main()
