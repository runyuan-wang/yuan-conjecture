#!/usr/bin/env python3
"""Holonomy certificate for the six-collision hard-sphere C12 gadget.

Each equal-mass binary collision is an orthogonal reflection on the total
velocity vector in R^(3N).  We follow physical particle labels through the
six central collisions, multiply the reflection matrices, and inspect the
fixed and rotating subspaces after quotienting the three momentum modes.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from verify_hard_sphere_c12_wiring import omega, tau, velocities


N = 6
DIM = 3 * N
order = [int(i) for i in np.argsort(tau)]

# Before the first collision, particle p occupies lane p.
particle_on_lane = list(range(N))
lane_of_particle = list(range(N))

V0 = velocities.reshape(-1).copy()
V = V0.copy()
H = np.eye(DIM)
events = []

for event in order:
    lane_i = event
    lane_j = (event + 1) % N
    particle_i = particle_on_lane[lane_i]
    particle_j = particle_on_lane[lane_j]

    normal = np.zeros(DIM)
    normal[3 * particle_i : 3 * particle_i + 3] = omega[event] / np.sqrt(2.0)
    normal[3 * particle_j : 3 * particle_j + 3] = -omega[event] / np.sqrt(2.0)
    reflection = np.eye(DIM) - 2.0 * np.outer(normal, normal)

    before = V.copy()
    V = reflection @ V
    H = reflection @ H

    # Central equal-mass collision swaps the two occupied lanes.
    particle_on_lane[lane_i], particle_on_lane[lane_j] = (
        particle_on_lane[lane_j],
        particle_on_lane[lane_i],
    )
    lane_of_particle[particle_i], lane_of_particle[particle_j] = (
        lane_j,
        lane_i,
    )

    expected = before.copy()
    expected[3 * particle_i : 3 * particle_i + 3] = before[
        3 * particle_j : 3 * particle_j + 3
    ]
    expected[3 * particle_j : 3 * particle_j + 3] = before[
        3 * particle_i : 3 * particle_i + 3
    ]
    events.append(
        {
            "event": event,
            "particle_pair": [particle_i, particle_j],
            "incoming_lane_pair": [lane_i, lane_j],
            "reflection_swap_residual": float(np.linalg.norm(V - expected)),
        }
    )

# Uniform velocity translations are fixed by every collision reflection.
momentum_basis = np.zeros((DIM, 3))
for p in range(N):
    momentum_basis[3 * p : 3 * p + 3, :] = np.eye(3) / np.sqrt(N)

assert np.linalg.norm(H.T @ H - np.eye(DIM)) < 1e-12
assert np.linalg.norm(H @ momentum_basis - momentum_basis) < 1e-12
assert all(e["reflection_swap_residual"] < 1e-10 for e in events)

# Build an orthonormal complement to the momentum modes.
q_full, _ = np.linalg.qr(
    np.concatenate([momentum_basis, np.eye(DIM)], axis=1), mode="complete"
)
# QR may reorder/sign the first columns; use an SVD nullspace instead.
_, _, vh = np.linalg.svd(momentum_basis.T)
relative_basis = vh[3:].T
H_relative = relative_basis.T @ H @ relative_basis

singular_I_minus_H = np.linalg.svd(np.eye(DIM) - H, compute_uv=False)
singular_relative = np.linalg.svd(
    np.eye(DIM - 3) - H_relative, compute_uv=False
)
tol = 1e-9
fixed_dim_total = int(np.sum(singular_I_minus_H < tol))
fixed_dim_relative = int(np.sum(singular_relative < tol))

eigenvalues = np.linalg.eigvals(H_relative)
angles = sorted(
    float(abs(np.angle(z)))
    for z in eigenvalues
    if abs(np.angle(z)) > 1e-8
)

certificate = {
    "status": "verified",
    "verifier": Path(__file__).name,
    "ambient_velocity_dimension": DIM,
    "chronological_event_order": order,
    "event_reflections": events,
    "final_particle_on_lane": particle_on_lane,
    "holonomy_orthogonality_residual": float(
        np.linalg.norm(H.T @ H - np.eye(DIM))
    ),
    "holonomy_determinant": float(np.linalg.det(H)),
    "momentum_fixed_residual": float(
        np.linalg.norm(H @ momentum_basis - momentum_basis)
    ),
    "fixed_subspace_dimension_total": fixed_dim_total,
    "fixed_subspace_dimension_after_momentum_quotient": fixed_dim_relative,
    "smallest_nonzero_singular_value_of_I_minus_H_relative": float(
        min(x for x in singular_relative if x >= tol)
    ),
    "nonzero_eigenangle_magnitudes_radians": angles,
    "initial_to_final_velocity_residual": float(np.linalg.norm(H @ V0 - V)),
    "interpretation": (
        "The collision history carries a nontrivial orthogonal holonomy. "
        "The three uniform-momentum modes are necessarily fixed; all further "
        "fixed modes are genuine degeneracies of this particular cycle."
    ),
}

Path("collision_holonomy_certificate.json").write_text(
    json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(certificate, indent=2))
