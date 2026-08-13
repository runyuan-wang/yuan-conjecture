#!/usr/bin/env python3
"""Independent PSITIP audit of the four-slope 7/4 entropy cycle.

This is deliberately not a new entropy solver.  It uses PSITIP 1.1.7 with
OR-Tools/GLOP to check three things:

1. Pairwise invertibility of the five linear forms alone does *not* let the
   Shannon LP prove the weighted 7/4 inequality.
2. Adding the two conditional-copy laws and the two exact cut identities from
   ``weighted_cycle_cutting.py`` does let the same LP prove it.
3. Two deliberately false conclusions are rejected, so the positive result is
   not caused by an inconsistent assumption set.

The arithmetic cut identities are independently checked on distributions by
``weighted_cycle_cutting.py``.  Here they enter as functional dependencies;
PSITIP audits only the information-inequality bookkeeping after that point.

Minimal reproducible environment:

    python -m venv --system-site-packages /tmp/psitip-audit
    /tmp/psitip-audit/bin/pip install --no-deps psitip==1.1.7
    /tmp/psitip-audit/bin/pip install ortools==9.15.6755
    /tmp/psitip-audit/bin/python psitip_cycle_audit.py
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from itertools import combinations
from time import perf_counter


def package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not-installed"


def main() -> None:
    try:
        from psitip import H, I, PsiOpts, eqdist, markov, rv, universe
    except ImportError as exc:  # pragma: no cover - dependency gate
        raise SystemExit(
            "PSITIP is required; use the pinned installation commands in the docstring"
        ) from exc

    PsiOpts.setting(solver="ortools.GLOP")
    print(f"psitip={package_version('psitip')} ortools={package_version('ortools')}")

    # A tiny official-style sanity check before the project-specific model.
    A, B = rv("sanity_A", "sanity_B")
    assert bool(I(A & B) >= 0)

    # First isolate what a single arithmetic atom tells an abstract entropy LP.
    X0, Y0, C0, D0, Z0 = rv("fd_X", "fd_Y", "fd_C", "fd_D", "fd_Z")
    forms = [X0, Y0, C0, D0, Z0]
    all_forms = X0 + Y0 + C0 + D0 + Z0
    fd_only = universe()
    for left, right in combinations(forms, 2):
        # Any two nonparallel forms determine all five.
        fd_only &= H(all_forms | left + right) == 0
    weighted_fd_target = 4 * H(Z0) <= (
        2 * H(X0) + 2 * H(Y0) + 2 * H(C0) + H(D0)
    )
    fd_only_proves = bool(fd_only >> weighted_fd_target)

    # Now encode the actual Katz--Tao conditional-copy molecule.
    X, Y, Yp, X1, Y1, C, Cp, D, Z, Zp = rv(
        "X", "Y", "Yp", "X1", "Y1", "C", "Cp", "D", "Z", "Zp"
    )
    region = universe()
    region &= H(C + D + Z | X + Y) == 0

    # Yp is a conditional copy of Y over X, and Cp carries the same sum map.
    region &= eqdist([X, Yp, Cp], [X, Y, C])
    region &= markov(Y, X, Yp)

    # (X1,Y1,Yp) is a conditional copy over W=(D,Yp).
    region &= eqdist([D, Yp, X1, Y1], [D, Yp, X, Y])
    region &= markov(X + Y, D + Yp, X1 + Y1)

    # Zp is the same difference map on q1=(X1,Yp).
    region &= eqdist([X1, Yp, Zp], [X, Y, Z])

    # h=(C,Cp,Y1) reveals Zp; (h,q1) reconstructs the whole copied state.
    region &= H(Zp | C + Cp + Y1) == 0
    region &= H(X + Y + Yp + X1 + Y1 | C + Cp + Y1 + X1 + Yp) == 0

    weighted_target = 4 * H(Z) <= (
        2 * H(X) + 2 * H(Y) + 2 * H(C) + H(D)
    )

    started = perf_counter()
    cycle_proves = bool(region >> weighted_target)
    proof_seconds = perf_counter() - started

    # The star law X=constant, Z=-Y refutes H(Z)<=H(X).  The second check
    # ensures the homogeneous cone has not accidentally collapsed to zero.
    false_hz_le_hx_accepted = bool(region >> (H(Z) <= H(X)))
    zero_cone_accepted = bool(region >> (H(X + Y) <= 0))

    assert not fd_only_proves
    assert cycle_proves
    assert not false_hz_le_hx_accepted
    assert not zero_cone_accepted

    print(f"fd_only_7_over_4_proved={fd_only_proves}")
    print(f"copy_cycle_7_over_4_proved={cycle_proves}")
    print(f"false_HZ_le_HX_accepted={false_hz_le_hx_accepted}")
    print(f"zero_cone_accepted={zero_cone_accepted}")
    print(f"proof_seconds={proof_seconds:.6f}")
    print("RESULT: independent entropy-bookkeeping audit passed")


if __name__ == "__main__":
    main()
