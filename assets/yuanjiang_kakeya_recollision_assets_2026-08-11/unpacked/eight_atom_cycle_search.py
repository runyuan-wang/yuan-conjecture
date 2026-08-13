#!/usr/bin/env python3
"""Search a genuinely new finite family beyond the known Katz--Tao 7/4 loop.

This is *not* another verification of the known four-atom proof.  It asks two
bounded questions suggested by the Yuanjiang collision-cycle conjecture:

1. Is any two-layer four-atom conditional-copy schema, using only X, Y, X+Y,
   X+2Y, strictly cheaper than 7/4 per atom?
2. Can one more conditional-copy layer turn *any* cost-7 four-atom seed into an
   eight-atom molecule strictly cheaper than 7/4 per atom?

The first search exhausts 1,725 four-atom generation schemas and 313,762
encodings.  It finds 16 fixed-label constraint sets at cost 7 and none at cost
at most 6.  For every one of those 16 seeds, the second search takes two
conditionally independent copies over a tuple W.  Every component of W is one
of the 16 allowed values

    projection r in {X,Y,X+Y,X+2Y} of atom i in {0,1,2,3}.

Choosing k components adds k cross-copy collision constraints and costs at
most k input entropies.  We then choose m individual input-projection values
as a cut/encoding h.  Linear algebra tests whether:

1. h determines X-Y for a collection of anchor atoms; and
2. h plus the full anchor atoms determines all eight atoms.

If so, conditional entropy pays for the anchor ambiguity and yields

    8 H(X-Y) <= (8+k+m) max_r H(X+rY).

Thus a strict improvement on 7/4 would require 8+k+m <= 13.  Across all seeds,
the program exhausts 110,160 choices of W and 8,931,424 encodings in that
strict-improvement budget.  It also verifies the known four-atom and doubled
7/4 baselines.

All rank tests are exact over Q, despite being computed modulo P=1,000,000,007.
Every original row has Euclidean norm at most sqrt(10), so Hadamard's bound
puts every relevant (at most 16 by 16) nonzero integer minor below 10^8 < P.
Reduction modulo P therefore cannot change any rank used by the search.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from functools import lru_cache
from itertools import combinations
from time import perf_counter
from typing import Iterable, Sequence


P = 1_000_000_007

SLOPES: dict[str, tuple[int, int]] = {
    "X": (1, 0),
    "Y": (0, 1),
    "C=X+Y": (1, 1),
    "D=X+2Y": (1, 2),
}
TARGET = (1, -1)

Edge = tuple[int, int, str]
Observation = tuple[int, str]


@lru_cache(maxsize=100_000)
def inverse(value: int) -> int:
    return pow(value % P, P - 2, P)


class RREF:
    """A tiny exact row-space basis over GF(P), maintained in reduced form."""

    __slots__ = ("dimension", "rows")

    def __init__(self, dimension: int, rows: Iterable[Sequence[int]] = ()) -> None:
        self.dimension = dimension
        self.rows: dict[int, list[int]] = {}
        for row in rows:
            self.add(row)

    def copy(self) -> "RREF":
        result = RREF(self.dimension)
        result.rows = {pivot: row.copy() for pivot, row in self.rows.items()}
        return result

    @property
    def rank(self) -> int:
        return len(self.rows)

    def reduce(self, vector: Sequence[int]) -> list[int]:
        work = [value % P for value in vector]
        for pivot in sorted(self.rows):
            coefficient = work[pivot]
            if coefficient:
                row = self.rows[pivot]
                work = [
                    (left - coefficient * right) % P
                    for left, right in zip(work, row)
                ]
        return work

    def contains(self, vector: Sequence[int]) -> bool:
        return not any(self.reduce(vector))

    def add(self, vector: Sequence[int]) -> bool:
        work = self.reduce(vector)
        try:
            pivot = next(index for index, value in enumerate(work) if value)
        except StopIteration:
            return False

        scale = inverse(work[pivot])
        work = [(value * scale) % P for value in work]

        # Back-eliminate the new pivot so reduce() is a canonical quotient map.
        for old_pivot, old_row in list(self.rows.items()):
            coefficient = old_row[pivot]
            if coefficient:
                self.rows[old_pivot] = [
                    (left - coefficient * right) % P
                    for left, right in zip(old_row, work)
                ]
        self.rows[pivot] = work
        return True


def projection_row(n: int, atom: int, coefficients: tuple[int, int]) -> list[int]:
    row = [0] * (2 * n)
    row[2 * atom], row[2 * atom + 1] = coefficients
    return row


def collision_row(n: int, left: int, right: int, label: str) -> list[int]:
    coefficients = SLOPES[label]
    row = projection_row(n, left, coefficients)
    other = projection_row(n, right, coefficients)
    return [a - b for a, b in zip(row, other)]


def constraint_rows(n: int, edges: Sequence[Edge]) -> list[list[int]]:
    return [collision_row(n, left, right, label) for left, right, label in edges]


KT4_EDGES: tuple[Edge, ...] = (
    (0, 1, "X"),
    (2, 3, "X"),
    (0, 2, "D=X+2Y"),
    (1, 3, "Y"),
)


def doubled_seed_edges(
    seed_edges: Sequence[Edge], cross_components: Sequence[Observation]
) -> list[Edge]:
    edges = list(seed_edges)
    edges.extend((left + 4, right + 4, label) for left, right, label in seed_edges)
    edges.extend((atom, atom + 4, label) for atom, label in cross_components)
    return edges


def doubled_kt_edges(cross_components: Sequence[Observation]) -> list[Edge]:
    return doubled_seed_edges(KT4_EDGES, cross_components)


@dataclass(frozen=True)
class QuotientModel:
    n: int
    constraint_rank: int
    quotient_dimension: int
    observation_labels: tuple[Observation, ...]
    observation_rows: tuple[tuple[int, ...], ...]
    target_rows: tuple[tuple[int, ...], ...]
    x_rows: tuple[tuple[int, ...], ...]
    y_rows: tuple[tuple[int, ...], ...]


def quotient_model(n: int, edges: Sequence[Edge]) -> QuotientModel:
    constraints = RREF(2 * n, constraint_rows(n, edges))
    free_columns = [column for column in range(2 * n) if column not in constraints.rows]

    def quotient(row: Sequence[int]) -> tuple[int, ...]:
        remainder = constraints.reduce(row)
        return tuple(remainder[column] for column in free_columns)

    # Equal observables modulo the constraints have identical values on every
    # supported state.  Keep one representative because all cost one unit.
    representatives: dict[tuple[int, ...], Observation] = {}
    for atom in range(n):
        for label, coefficients in SLOPES.items():
            vector = quotient(projection_row(n, atom, coefficients))
            if any(vector):
                representatives.setdefault(vector, (atom, label))

    observation_rows = tuple(representatives)
    observation_labels = tuple(representatives[row] for row in observation_rows)
    target_rows = tuple(
        quotient(projection_row(n, atom, TARGET)) for atom in range(n)
    )
    x_rows = tuple(
        quotient(projection_row(n, atom, (1, 0))) for atom in range(n)
    )
    y_rows = tuple(
        quotient(projection_row(n, atom, (0, 1))) for atom in range(n)
    )
    return QuotientModel(
        n=n,
        constraint_rank=constraints.rank,
        quotient_dimension=len(free_columns),
        observation_labels=observation_labels,
        observation_rows=observation_rows,
        target_rows=target_rows,
        x_rows=x_rows,
        y_rows=y_rows,
    )


@dataclass(frozen=True)
class EncodingCheck:
    feasible: bool
    reconstruction_rank: int
    anchors: tuple[int, ...]


def check_encoding(model: QuotientModel, indices: Sequence[int]) -> EncodingCheck:
    basis = RREF(
        model.quotient_dimension,
        (model.observation_rows[index] for index in indices),
    )
    anchors = tuple(
        atom
        for atom, target in enumerate(model.target_rows)
        if basis.contains(target)
    )
    reconstruction = basis.copy()
    for atom in anchors:
        reconstruction.add(model.x_rows[atom])
        reconstruction.add(model.y_rows[atom])
    return EncodingCheck(
        feasible=reconstruction.rank == model.quotient_dimension,
        reconstruction_rank=reconstruction.rank,
        anchors=anchors,
    )


@dataclass(frozen=True)
class EncodingSearch:
    feasible: bool
    tested: int
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    best_rank: int
    best_encoding: tuple[Observation, ...]
    best_anchors: tuple[int, ...]
    quotient_dimension: int
    constraint_rank: int


def search_encodings(model: QuotientModel, max_size: int) -> EncodingSearch:
    tested = 0
    best_rank = -1
    best_indices: tuple[int, ...] = ()
    best_anchors: tuple[int, ...] = ()
    count = len(model.observation_rows)

    for size in range(max_size + 1):
        for indices in combinations(range(count), size):
            tested += 1
            result = check_encoding(model, indices)
            if (result.reconstruction_rank, len(result.anchors)) > (
                best_rank,
                len(best_anchors),
            ):
                best_rank = result.reconstruction_rank
                best_indices = indices
                best_anchors = result.anchors
            if result.feasible:
                labels = tuple(model.observation_labels[index] for index in indices)
                return EncodingSearch(
                    feasible=True,
                    tested=tested,
                    encoding=labels,
                    anchors=result.anchors,
                    best_rank=best_rank,
                    best_encoding=labels,
                    best_anchors=result.anchors,
                    quotient_dimension=model.quotient_dimension,
                    constraint_rank=model.constraint_rank,
                )

    return EncodingSearch(
        feasible=False,
        tested=tested,
        encoding=(),
        anchors=(),
        best_rank=best_rank,
        best_encoding=tuple(model.observation_labels[index] for index in best_indices),
        best_anchors=best_anchors,
        quotient_dimension=model.quotient_dimension,
        constraint_rank=model.constraint_rank,
    )


def format_observations(observations: Sequence[Observation]) -> str:
    return ", ".join(f"{label}[{atom}]" for atom, label in observations) or "(none)"


def verify_baselines() -> None:
    print("[baseline: known four-atom molecule]")
    model4 = quotient_model(4, KT4_EDGES)
    result4 = search_encodings(model4, max_size=3)
    assert result4.feasible
    print(
        f"constraint rank / quotient dimension = "
        f"{result4.constraint_rank}/{result4.quotient_dimension}"
    )
    print(f"minimal encoding found = {format_observations(result4.encoding)}")
    print(f"anchors = {result4.anchors}")
    print("cost = 4 collision components + 3 encoding components = 7; ratio = 7/4")

    print("\n[baseline: two independent four-atom molecules]")
    model8 = quotient_model(8, doubled_kt_edges(()))
    known_encoding: tuple[Observation, ...] = (
        (0, "C=X+Y"),
        (1, "C=X+Y"),
        (2, "Y"),
        (4, "C=X+Y"),
        (5, "C=X+Y"),
        (6, "Y"),
    )
    label_to_index = {
        label: index for index, label in enumerate(model8.observation_labels)
    }
    indices = tuple(label_to_index[label] for label in known_encoding)
    check = check_encoding(model8, indices)
    assert check.feasible
    print(f"encoding = {format_observations(known_encoding)}")
    print(f"anchors = {check.anchors}")
    print("cost = 8 collision components + 6 encoding components = 14; ratio = 7/4")


@dataclass(frozen=True)
class FourAtomSchema:
    first_copy: tuple[str, ...]
    second_copy: tuple[Observation, ...]
    collision_cost: int
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    total_cost: int


def four_atom_edges(
    first_copy: Sequence[str], second_copy: Sequence[Observation]
) -> list[Edge]:
    edges: list[Edge] = []
    for label in first_copy:
        edges.append((0, 1, label))
        edges.append((2, 3, label))
    edges.extend((atom, atom + 2, label) for atom, label in second_copy)
    return edges


def general_four_atom_search() -> list[FourAtomSchema]:
    """Exhaust every two-layer four-atom copy schema with total cost <= 7."""
    first_universe = tuple(SLOPES)
    second_universe = tuple(
        (atom, label) for atom in range(2) for label in SLOPES
    )
    schemas = 0
    encodings = 0
    feasible: list[FourAtomSchema] = []
    strict: list[FourAtomSchema] = []
    started = perf_counter()

    print("\n[all two-layer four-atom schemas: target total cost <= 7]")
    for first_size in range(4):
        for first in combinations(first_universe, first_size):
            max_second_size = 7 - 2 * first_size
            for second_size in range(min(8, max_second_size) + 1):
                collision_cost = 2 * first_size + second_size
                max_encoding_size = 7 - collision_cost
                for second in combinations(second_universe, second_size):
                    schemas += 1
                    model = quotient_model(4, four_atom_edges(first, second))
                    result = search_encodings(model, max_encoding_size)
                    encodings += result.tested
                    if not result.feasible:
                        continue
                    schema = FourAtomSchema(
                        first_copy=tuple(first),
                        second_copy=tuple(second),
                        collision_cost=collision_cost,
                        encoding=result.encoding,
                        anchors=result.anchors,
                        total_cost=collision_cost + len(result.encoding),
                    )
                    feasible.append(schema)
                    if schema.total_cost < 7:
                        strict.append(schema)

    elapsed = perf_counter() - started
    distribution = Counter(
        (schema.collision_cost, len(schema.encoding), schema.total_cost)
        for schema in feasible
    )
    print(f"schemas exhausted = {schemas}")
    print(f"encodings tested = {encodings}")
    print(f"feasible schemas with total cost <=7 = {len(feasible)}")
    print(f"strict improvements with total cost <=6 = {len(strict)}")
    print(f"elapsed seconds = {elapsed:.3f}")
    print("feasible cost distribution (collisions, encoding, total -> count):")
    for key, count in sorted(distribution.items()):
        print(f"  {key} -> {count}")

    for index, schema in enumerate(feasible[:8], start=1):
        print(f"representative {index}:")
        print(f"  first copy W1 = {', '.join(schema.first_copy) or '(none)'}")
        print(f"  second copy W2 = {format_observations(schema.second_copy)}")
        print(f"  encoding h = {format_observations(schema.encoding)}")
        print(f"  anchors = {schema.anchors}; total cost = {schema.total_cost}")

    if strict:
        print("WARNING: a putative four-atom improvement was found and needs audit")
    else:
        print(
            "Result: Katz--Tao's 7/4 cost is optimal inside this entire "
            "two-layer, individual-projection copy/encoding family."
        )
    return feasible


@dataclass(frozen=True)
class GlobalNearMiss:
    reconstruction_rank: int
    quotient_dimension: int
    cross_components: tuple[Observation, ...]
    encoding: tuple[Observation, ...]
    anchors: tuple[int, ...]
    total_cost: int


@dataclass(frozen=True)
class SeedSearchResult:
    candidate: GlobalNearMiss | None
    near_miss: GlobalNearMiss
    graph_count: int
    encoding_count: int
    elapsed: float


def search_one_eight_atom_seed(seed_edges: Sequence[Edge]) -> SeedSearchResult:
    cross_universe = tuple(
        (atom, label) for atom in range(4) for label in SLOPES
    )
    graph_count = 0
    encoding_count = 0
    near_miss: GlobalNearMiss | None = None
    started = perf_counter()

    for cross_size in range(6):
        max_encoding_size = 5 - cross_size

        for cross in combinations(cross_universe, cross_size):
            graph_count += 1
            model = quotient_model(8, doubled_seed_edges(seed_edges, cross))
            result = search_encodings(model, max_encoding_size)
            encoding_count += result.tested

            candidate = GlobalNearMiss(
                reconstruction_rank=result.best_rank,
                quotient_dimension=result.quotient_dimension,
                cross_components=tuple(cross),
                encoding=result.best_encoding,
                anchors=result.best_anchors,
                total_cost=8 + cross_size + len(result.best_encoding),
            )
            if near_miss is None or (
                candidate.reconstruction_rank - candidate.quotient_dimension,
                len(candidate.anchors),
                -candidate.total_cost,
            ) > (
                near_miss.reconstruction_rank - near_miss.quotient_dimension,
                len(near_miss.anchors),
                -near_miss.total_cost,
            ):
                near_miss = candidate

            if result.feasible:
                elapsed = perf_counter() - started
                return SeedSearchResult(
                    candidate=candidate,
                    near_miss=candidate,
                    graph_count=graph_count,
                    encoding_count=encoding_count,
                    elapsed=elapsed,
                )

    assert near_miss is not None
    elapsed = perf_counter() - started
    return SeedSearchResult(
        candidate=None,
        near_miss=near_miss,
        graph_count=graph_count,
        encoding_count=encoding_count,
        elapsed=elapsed,
    )


def strict_improvement_search() -> None:
    print("\n[strict-improvement search: target total cost <= 13]")
    result = search_one_eight_atom_seed(KT4_EDGES)
    if result.candidate is not None:
        candidate = result.candidate
        print("STRICT IMPROVEMENT CANDIDATE FOUND")
        print(f"cross W = {format_observations(candidate.cross_components)}")
        print(f"encoding h = {format_observations(candidate.encoding)}")
        print(f"anchors = {candidate.anchors}")
        print(
            f"cost = {candidate.total_cost}; resulting ratio = "
            f"{candidate.total_cost}/8 = {candidate.total_cost/8:.6f}"
        )
        return

    near_miss = result.near_miss
    print("\nNO STRICT IMPROVEMENT IN THE EXHAUSTED FAMILY")
    print(f"graphs exhausted = {result.graph_count}")
    print(f"encodings tested = {result.encoding_count}")
    print(f"elapsed seconds = {result.elapsed:.3f}")
    print("best near-miss within the <=13 budget:")
    print(f"  cross W = {format_observations(near_miss.cross_components)}")
    print(f"  encoding h = {format_observations(near_miss.encoding)}")
    print(f"  anchors = {near_miss.anchors}")
    print(
        f"  reconstruction rank = {near_miss.reconstruction_rank}/"
        f"{near_miss.quotient_dimension}"
    )
    print(f"  accounted total cost = {near_miss.total_cost}")
    print(
        "Interpretation: one extra conditional-copy layer of this restricted form "
        "cannot beat 7/4; this is a finite-family no-go result, not a general no-go theorem."
    )


def all_seed_eight_atom_search(feasible: Sequence[FourAtomSchema]) -> None:
    """Test an extra copy layer over every four-atom cost-7 seed found above."""
    # Different generation histories can yield the same colored constraint set.
    unique: dict[tuple[Edge, ...], FourAtomSchema] = {}
    for schema in feasible:
        edges = tuple(sorted(four_atom_edges(schema.first_copy, schema.second_copy)))
        unique.setdefault(edges, schema)

    print("\n[all cost-7 seeds, one extra copy layer: target total cost <= 13]")
    print(
        f"feasible generation schemas = {len(feasible)}; "
        f"distinct fixed-label constraint sets = {len(unique)}"
    )
    total_graphs = 0
    total_encodings = 0
    global_near_miss: GlobalNearMiss | None = None
    global_near_seed = 0
    started = perf_counter()

    for seed_index, (edges, schema) in enumerate(unique.items(), start=1):
        result = search_one_eight_atom_seed(edges)
        total_graphs += result.graph_count
        total_encodings += result.encoding_count
        candidate = result.near_miss
        if global_near_miss is None or (
            candidate.reconstruction_rank - candidate.quotient_dimension,
            len(candidate.anchors),
            -candidate.total_cost,
        ) > (
            global_near_miss.reconstruction_rank
            - global_near_miss.quotient_dimension,
            len(global_near_miss.anchors),
            -global_near_miss.total_cost,
        ):
            global_near_miss = candidate
            global_near_seed = seed_index

        print(
            f"seed {seed_index:02d}/{len(unique)}: "
            f"W1={','.join(schema.first_copy)}; "
            f"W2={format_observations(schema.second_copy)}; "
            f"graphs/encodings={result.graph_count}/{result.encoding_count}; "
            f"elapsed={result.elapsed:.3f}s; "
            f"{'CANDIDATE' if result.candidate else 'no candidate'}",
            flush=True,
        )
        if result.candidate is not None:
            found = result.candidate
            print("STRICT IMPROVEMENT CANDIDATE FOUND; stop for independent audit")
            print(f"seed edges = {edges}")
            print(f"cross W = {format_observations(found.cross_components)}")
            print(f"encoding h = {format_observations(found.encoding)}")
            print(f"anchors = {found.anchors}")
            print(
                f"cost = {found.total_cost}; ratio = "
                f"{found.total_cost}/8 = {found.total_cost/8:.6f}"
            )
            return

    assert global_near_miss is not None
    elapsed = perf_counter() - started
    print("NO STRICT IMPROVEMENT FOR ANY COST-7 FOUR-ATOM SEED")
    print(f"seed constraint sets exhausted = {len(unique)}")
    print(f"graphs exhausted = {total_graphs}")
    print(f"encodings tested = {total_encodings}")
    print(f"elapsed seconds = {elapsed:.3f}")
    print(f"global near-miss came from seed {global_near_seed}")
    print(
        f"  reconstruction rank = {global_near_miss.reconstruction_rank}/"
        f"{global_near_miss.quotient_dimension}; "
        f"anchors = {global_near_miss.anchors}; "
        f"cost = {global_near_miss.total_cost}"
    )
    print(f"  cross W = {format_observations(global_near_miss.cross_components)}")
    print(f"  encoding h = {format_observations(global_near_miss.encoding)}")
    print(
        "Interpretation: within this fully exhausted two-layer seed family, "
        "adding one more individual-projection conditional-copy layer cannot "
        "improve 7/4.  This remains a bounded-family no-go result."
    )


def main() -> None:
    print("Yuanjiang eight-atom collision-cycle search")
    print(f"rank field prime P = {P}; Hadamard safety bound = 100,000,000")
    verify_baselines()
    feasible = general_four_atom_search()
    all_seed_eight_atom_search(feasible)


if __name__ == "__main__":
    main()
