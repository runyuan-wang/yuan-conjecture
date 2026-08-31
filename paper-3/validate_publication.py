#!/usr/bin/env python3
"""Validate the public Paper-3 package.

Static checks use only the Python standard library.  --run additionally executes
both numerical programs and therefore requires NumPy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
ASSET = ROOT / "assets" / "legal_collision_movie_2026-08-29"

FROZEN = {
    "main.tex": "dc805326178f57bda5a37f436d73cf37ab2e06d8e2409c1bcd5b1dbac0f0bb91",
    "main.pdf": "7c2b63a6241acc3efb93781ebf522afa493bb68ba3953f51f55b4a73ee0cf081",
}

FAMILIES = [
    "strict_C12_reference",
    "generic_rational_C18",
    "time_separated_hinge_star",
    "simultaneous_bush_K5",
    "planar_weave_K6",
    "regulus_K4_4",
    "sheared_regulus_K4_4",
]

OPTIMA = {
    "simultaneous_bush_K5_no_cap": 2,
    "planar_weave_K6_cap_2": 6,
    "planar_weave_K6_cap_4": 12,
    "regulus_K4_4_no_cap": 12,
    "regulus_K4_4_cap_2": 8,
    "regulus_K4_4_cap_4": 12,
}

FORBIDDEN_PUBLIC_MARKERS = (
    "/Users/",
    ".lingtai/",
    "changshengwu4802",
    "505150827",
    "parent-review-",
    "PARENT_ACCEPTANCE",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def check_manifest() -> None:
    manifest = ASSET / "SHA256SUMS_legal_collision_movie_2026-08-29.txt"
    rows = []
    for raw in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = raw.split(maxsplit=1)
        name = name.strip()
        rows.append(name)
        path = ASSET / name
        require(path.is_file(), f"manifest file missing: {name}")
        require(sha256(path) == digest, f"manifest hash mismatch: {name}")
    require(len(rows) == 6 and len(set(rows)) == 6, "companion manifest must name exactly six unique files")


def check_frozen_artifacts() -> None:
    for rel, digest in FROZEN.items():
        path = ROOT / rel
        require(path.is_file(), f"missing frozen artifact: {rel}")
        require(sha256(path) == digest, f"frozen hash mismatch: {rel}")


def check_public_manifest() -> None:
    manifest = ROOT / "SHA256SUMS.txt"
    require(manifest.is_file(), "missing public SHA256SUMS.txt")
    listed: dict[str, str] = {}
    for raw in manifest.read_text(encoding="utf-8").splitlines():
        digest, rel = raw.split(maxsplit=1)
        rel = rel.strip()
        require(rel != "SHA256SUMS.txt", "public manifest must not include itself")
        require(rel not in listed, f"duplicate public-manifest path: {rel}")
        path = ROOT / rel
        require(path.is_file(), f"public-manifest file missing: {rel}")
        require(sha256(path) == digest, f"public-manifest hash mismatch: {rel}")
        listed[rel] = digest
    expected = {
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS.txt"
        and "__pycache__" not in path.parts
        and path.suffix != ".pyc"
    }
    require(set(listed) == expected, "public manifest file set is incomplete or contains extras")


def check_claims(certificate: dict, verifier: dict) -> None:
    require(certificate["metadata"]["randomness"] == "none", "compiler must declare no randomness")
    require([row["name"] for row in certificate["families"]] == FAMILIES, "unexpected family inventory/order")
    require(len(certificate["regulus_scaling_probe"]) == 7, "expected seven scaling-probe rows")

    families = {row["name"]: row for row in certificate["families"]}
    require(families["strict_C12_reference"]["full_network"]["wireable"] is True, "C12 control must pass")
    require(families["generic_rational_C18"]["full_network"]["rank"] == 27, "C18 rank mismatch")
    require(families["time_separated_hinge_star"]["full_network"]["rank"] == 21, "hinge rank mismatch")
    require(families["planar_weave_K6"]["full_network"]["wireable"] is True, "planar weave must pass")
    require(families["planar_weave_K6"]["full_network"]["left_null_dimension"] == 20, "planar left-null mismatch")
    require(families["regulus_K4_4"]["full_network"]["wireable"] is False, "regulus must fail")
    require(families["regulus_K4_4"]["full_network"]["rank"] == 59, "regulus rank mismatch")
    require(families["regulus_K4_4"]["full_network"]["left_null_dimension"] == 13, "regulus left-null mismatch")
    require(families["sheared_regulus_K4_4"]["full_network"]["wireable"] is False, "sheared regulus must fail")

    require(verifier["status"] == "PASS", "independent verifier status mismatch")
    require(verifier["exhaustive_optima"] == OPTIMA, "finite subset optima mismatch")
    checks = verifier["full_network_checks"]
    require(checks["strict_C12_reference"]["ok"] is True and checks["strict_C12_reference"]["rank"] == 18, "verifier C12 mismatch")
    require(checks["planar_weave_K6"]["ok"] is True and checks["planar_weave_K6"]["left_null"] == 20, "verifier planar mismatch")
    require(checks["regulus_K4_4"]["ok"] is False and checks["regulus_K4_4"]["rank"] == 59, "verifier regulus mismatch")
    require(checks["sheared_regulus_K4_4"]["ok"] is False, "verifier sheared-regulus mismatch")


def check_manuscript_boundaries() -> None:
    tex = (ROOT / "main.tex").read_text(encoding="utf-8")
    required = (
        "not a theorem about general hard-sphere dynamics",
        "No construction maps arbitrary Kakeya sets",
        "not a proof bridge",
        "An interface problem remains open",
        "not a continuum theorem",
    )
    for phrase in required:
        require(phrase in tex, f"missing scientific-boundary phrase: {phrase}")


def check_public_hygiene() -> None:
    text_suffixes = {".md", ".tex", ".py", ".json", ".txt"}
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in text_suffixes:
            continue
        if path.name == "validate_publication.py":
            continue
        text = path.read_text(encoding="utf-8", errors="strict")
        rel = path.relative_to(ROOT)
        for marker in FORBIDDEN_PUBLIC_MARKERS:
            require(marker not in text, f"forbidden internal marker {marker!r} in {rel}")


def categorical_signature(certificate: dict, verifier: dict) -> dict:
    families = {}
    for row in certificate["families"]:
        full = row["full_network"]
        families[row["name"]] = {
            "wireable": full.get("wireable"),
            "rank": full.get("rank"),
            "left_null_dimension": full.get("left_null_dimension"),
            "reason": full.get("reason"),
        }
    return {
        "families": families,
        "scaling_n": [row["n_per_ruling"] for row in certificate["regulus_scaling_probe"]],
        "verifier_status": verifier["status"],
        "verifier_optima": verifier["exhaustive_optima"],
        "verifier_ranks": {
            name: {"ok": row.get("ok"), "rank": row.get("rank"), "left_null": row.get("left_null"), "reason": row.get("reason")}
            for name, row in verifier["full_network_checks"].items()
        },
    }


def run_programs(frozen_certificate: dict, frozen_verifier: dict) -> None:
    compiler = ASSET / "legal_collision_movie_compiler.py"
    verifier = ASSET / "verify_legal_collision_movie.py"
    with tempfile.TemporaryDirectory(prefix="paper3-validation-") as tmp:
        tmpdir = Path(tmp)
        generated_certificate = tmpdir / "certificate.json"
        generated_verifier = tmpdir / "verifier.json"
        subprocess.run([sys.executable, str(compiler), "--output", str(generated_certificate)], check=True)
        subprocess.run([sys.executable, str(verifier), "--output", str(generated_verifier)], check=True)
        fresh_certificate = read_json(generated_certificate)
        fresh_verifier = read_json(generated_verifier)
    require(
        categorical_signature(fresh_certificate, fresh_verifier)
        == categorical_signature(frozen_certificate, frozen_verifier),
        "fresh executable rerun changed categorical results",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="also execute compiler and independent verifier (requires NumPy)")
    args = parser.parse_args()

    check_frozen_artifacts()
    check_manifest()
    check_public_manifest()
    certificate = read_json(ASSET / "legal_collision_movie_certificate_2026-08-29.json")
    verifier = read_json(ASSET / "legal_collision_movie_verifier_2026-08-29.json")
    check_claims(certificate, verifier)
    check_manuscript_boundaries()
    check_public_hygiene()
    print("PASS: static Paper-3 publication validation")

    if args.run:
        run_programs(certificate, verifier)
        print("PASS: fresh compiler/verifier categorical agreement")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
