#!/usr/bin/env python3
"""Finite verification of the rooted-C12 greedy deletion lemma.

The general lemma is elementary: repeatedly choose a C12 in the current
bipartite graph and delete one of its edges.  The chosen cycles are triangular
against their deleted root edges, hence linearly independent over every field;
the final graph is C12-free.

This script checks the construction and its GF(2) rank on K_{6,6}.  It is a
reproducibility aid, not a substitute for the one-paragraph general proof.
"""

from __future__ import annotations

import json
from pathlib import Path


LEFT = tuple(range(6))
RIGHT = tuple(range(6, 12))
initial_edges = {(u, v) for u in LEFT for v in RIGHT}


def adjacency(edges):
    out = {v: set() for v in LEFT + RIGHT}
    for u, v in edges:
        out[u].add(v)
        out[v].add(u)
    return out


def find_c12(edges):
    """Return one simple 12-cycle as an edge set, or None."""

    adj = adjacency(edges)
    start = LEFT[0]

    def dfs(path):
        current = path[-1]
        if len(path) == 12:
            if start not in adj[current]:
                return None
            cyc = set()
            for a, b in zip(path, path[1:] + [start]):
                cyc.add((a, b) if a in LEFT else (b, a))
            return cyc
        for nxt in sorted(adj[current]):
            if nxt in path:
                continue
            found = dfs(path + [nxt])
            if found is not None:
                return found
        return None

    # A C12 in this 6+6 graph must use every left vertex, hence it contains 0.
    return dfs([start])


edges = set(initial_edges)
cycles = []
roots = []
while True:
    cycle = find_c12(edges)
    if cycle is None:
        break
    root = min(cycle)
    cycles.append(cycle)
    roots.append(root)
    edges.remove(root)


def gf2_rank(rows):
    rows = list(rows)
    rank = 0
    pivot = 0
    while pivot < len(rows):
        bit = rows[pivot].bit_length() - 1
        if bit < 0:
            pivot += 1
            continue
        for j in range(len(rows)):
            if j != pivot and ((rows[j] >> bit) & 1):
                rows[j] ^= rows[pivot]
        rank += 1
        pivot += 1
    return rank


edge_order = sorted(initial_edges)
edge_index = {e: i for i, e in enumerate(edge_order)}
bit_rows = []
for cycle in cycles:
    row = 0
    for edge in cycle:
        row |= 1 << edge_index[edge]
    bit_rows.append(row)

# In deletion order, cycle i cannot contain an earlier root and contains its own.
triangular = all(
    roots[i] in cycles[i] and all(roots[j] not in cycles[i] for j in range(i))
    for i in range(len(cycles))
)

certificate = {
    "status": "verified",
    "verifier": Path(__file__).name,
    "graph": "K_{6,6}",
    "initial_edge_count": len(initial_edges),
    "deleted_root_count": len(roots),
    "residual_edge_count": len(edges),
    "rooted_cycle_count": len(cycles),
    "gf2_cycle_rank": gf2_rank(bit_rows),
    "triangular_against_root_order": triangular,
    "residual_has_c12": find_c12(edges) is not None,
    "identity_deleted_plus_residual": len(roots) + len(edges),
    "certified_consequence": (
        "The greedy roots are distinct, their selected C12 incidence vectors "
        "are linearly independent, and deletion stops at a C12-free graph."
    ),
}

assert triangular
assert certificate["gf2_cycle_rank"] == len(cycles)
assert certificate["residual_has_c12"] is False
assert len(roots) + len(edges) == len(initial_edges)

output = Path(__file__).with_name("rooted_c12_deletion_certificate.json")
output.write_text(json.dumps(certificate, indent=2) + "\n")
print(json.dumps(certificate, indent=2))
