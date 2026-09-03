#!/usr/bin/env python3
"""Deterministic source, certificate, citation, and manuscript checks.

Only the Python standard library is used.  The script reads but never writes
the accepted experiment package supplied with --source.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import struct
import sys
from pathlib import Path

CANDIDATE_HASH_DOMAIN = b"elastic-rectangle-candidate-v1\0uint32be(+,+,-,-)\0"
PRIMES = (2_147_483_647, 2_147_483_629)
CITATION_ALLOWLIST = {
    "BobylevCercignani1999DVM",
    "BernhoffVinerean2016Mixtures",
    "GoelMontenegroTetali2006SpectralProfile",
}

EXPECTED = {
    2: {
        "candidate": "collision_rectangle_expander/results/candidates_m2/m2_residual_greedy_d6_t0_5bcda3b714a4.json",
        "cert": "PARENT_RERUNS/certificates/m2.json",
        "sha": "5bcda3b714a4874bb5efe9b0814a21d2b0e72d5d41a645b73a2150fb360b22e5",
        "n": 125, "rows": 188, "dmin": 4, "dmax": 10, "rank": 120,
        "stored_gap": 0.13756732710431335,
        "fresh_gap": 0.13756732710431338,
    },
    3: {
        "candidate": "collision_rectangle_expander/results/candidates_m3/m3_residual_greedy_d6_t0_ba712eefee09.json",
        "cert": "PARENT_RERUNS/certificates/m3.json",
        "sha": "ba712eefee0910e23a3ca852badf083f98aa5ca71af84ae8b36e09343d959edc",
        "n": 343, "rows": 515, "dmin": 3, "dmax": 10, "rank": 338,
        "stored_gap": 0.10883600589808415,
        "fresh_gap": 0.10883600589808418,
    },
    4: {
        "candidate": "collision_rectangle_expander/results/candidates_m4/m4_residual_greedy_d6_t0_f5a7cb6aead6.json",
        "cert": "PARENT_RERUNS/certificates/m4.json",
        "sha": "f5a7cb6aead6eb38f73dadce7a757fbf5fc6e2b31475a784b9653e748e6a597c",
        "n": 729, "rows": 1094, "dmin": 3, "dmax": 10, "rank": 724,
        "stored_gap": 0.10730444046320077,
        "fresh_gap": 0.107304440462551,
    },
    5: {
        "candidate": "collision_rectangle_expander/results/candidates_m5/m5_residual_greedy_d6_t0_8cfe7af58e85.json",
        "cert": "PARENT_RERUNS/certificates/m5.json",
        "sha": "8cfe7af58e859d13e076f97ea9e77b5fe30648a3e2c8ffb69ff3991ae3f1f0aa",
        "n": 1331, "rows": 1997, "dmin": 3, "dmax": 10, "rank": 1326,
        "stored_gap": 0.10446152042043529,
        "fresh_gap": 0.10446152042018168,
    },
}

PUBLIC_FILES = [
    "PAPER_PLAN.md",
    "README.md",
    "SOURCE_MATRIX.md",
    "main.tex",
    "SUMMARY.zh-CN.md",
    "REPORT.md",
    "references.bib",
]


class Audit:
    def __init__(self) -> None:
        self.checks = 0
        self.notes: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            raise AssertionError(message)

    def note(self, message: str) -> None:
        self.notes.append(message)
        print(f"NOTE {message}")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def coords(m: int, idx: int) -> tuple[int, int, int]:
    side = 2 * m + 1
    xq, z = divmod(idx, side)
    x, y = divmod(xq, side)
    return x - m, y - m, z - m


def add(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def sub(a: tuple[int, int, int], b: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def dot(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    return sum(x * y for x, y in zip(a, b))


def hash_rows(rows: list[list[int]]) -> str:
    h = hashlib.sha256(CANDIDATE_HASH_DOMAIN)
    for row in rows:
        h.update(struct.pack(">4I", *row))
    return h.hexdigest()


def rank_restricted_minor_mod2(rows: list[list[int]], row_ids: list[int], col_ids: list[int]) -> int:
    colmap = {c: j for j, c in enumerate(col_ids)}
    pivots: dict[int, int] = {}
    for rid in row_ids:
        bits = 0
        for c in rows[rid]:
            if c in colmap:
                bits ^= 1 << colmap[c]
        while bits:
            pivot = bits.bit_length() - 1
            if pivot in pivots:
                bits ^= pivots[pivot]
            else:
                pivots[pivot] = bits
                break
    return len(pivots)


def audit_candidate_and_certificate(audit: Audit, source: Path, m: int, spec: dict) -> None:
    obj = load_json(source / spec["candidate"])
    rows = obj["rectangles"]
    n = spec["n"]
    audit.require(obj["m"] == m, f"candidate m mismatch for m={m}")
    audit.require(obj["vertex_count"] == n == (2 * m + 1) ** 3, f"vertex count mismatch m={m}")
    audit.require(len(rows) == spec["rows"], f"row count mismatch m={m}")
    audit.require(len({tuple(r) for r in rows}) == len(rows), f"duplicate row m={m}")

    degrees = [0] * n
    modes = [[1, *coords(m, i), dot(coords(m, i), coords(m, i))] for i in range(n)]
    for rid, row in enumerate(rows):
        audit.require(len(row) == 4 and all(type(c) is int for c in row), f"malformed row {rid} m={m}")
        i, j, k, ell = row
        audit.require(0 <= i < j < n and 0 <= k < ell < n and (i, j) < (k, ell),
                      f"noncanonical row {rid} m={m}")
        audit.require(len(set(row)) == 4, f"degenerate row {rid} m={m}")
        u, v, r, s = (coords(m, c) for c in row)
        audit.require(add(u, v) == add(r, s), f"momentum failure row {rid} m={m}")
        audit.require(dot(sub(u, v), sub(u, v)) == dot(sub(r, s), sub(r, s)),
                      f"diagonal-length failure row {rid} m={m}")
        p, q = sub(r, u), sub(s, u)
        audit.require(p != (0, 0, 0) and q != (0, 0, 0) and dot(p, q) == 0,
                      f"orthogonal-side failure row {rid} m={m}")
        residual = [modes[i][c] + modes[j][c] - modes[k][c] - modes[ell][c] for c in range(5)]
        audit.require(residual == [0] * 5, f"known-mode residual row {rid} m={m}")
        for c in row:
            degrees[c] += 1

    digest = hash_rows(rows)
    audit.require(digest == obj["candidate_sha256"] == spec["sha"], f"row hash mismatch m={m}")
    dmin, dmax = min(degrees), max(degrees)
    audit.require(dmin == spec["dmin"] and dmax == spec["dmax"], f"degree mismatch m={m}")
    audit.require(0 not in degrees, f"degree-zero vertex from raw rows m={m}")
    bs = obj["basic_statistics"]
    audit.require(bs["min_incidence_degree"] == dmin, f"stored dmin mismatch m={m}")
    audit.require(bs["max_incidence_degree"] == dmax, f"stored dmax mismatch m={m}")
    audit.require(bs["zero_degree_vertices"] == 0, f"zero degree m={m}")
    audit.require(math.isclose(obj["normalized_spectrum"]["gap_singular_value"], spec["stored_gap"], rel_tol=0, abs_tol=1e-15),
                  f"stored gap mismatch m={m}")

    cert = load_json(source / spec["cert"])
    audit.require(cert["schema"] == "exact_elastic_rectangle_kernel_certificate.v1", f"certificate schema m={m}")
    audit.require(cert["m"] == m and cert["candidate_sha256"] == digest, f"certificate identity m={m}")
    audit.require(cert["vertex_count"] == n and cert["rectangle_count"] == len(rows), f"certificate size m={m}")
    audit.require(cert["candidate_rows_unique"] is True, f"certificate uniqueness m={m}")
    audit.require(cert["integer_rows_formula_and_canonical_pool_membership_all_verified"] is True,
                  f"certificate row audit m={m}")
    audit.require(cert["known_mode_integer_residual_max_abs"] == 0, f"certificate mode residual m={m}")
    audit.require(abs(cert["known_mode_anchor_determinant_exact"]) == 2, f"anchor determinant m={m}")
    audit.require(cert["rank_over_Q"] == cert["rank_upper_bound_over_Q_from_known_kernel"] == spec["rank"] == n - 5,
                  f"rational rank m={m}")
    audit.require(cert["kernel_dimension_over_Q"] == 5 and cert["known_modes_span_entire_rational_kernel"] is True,
                  f"rational kernel m={m}")
    screens = cert["large_prime_rank_screens"]
    audit.require([(x["prime"], x["rank"], x["full_n_minus_5"]) for x in screens]
                  == [(p, n - 5, True) for p in PRIMES], f"prime screens m={m}")
    audit.require(cert["large_primes_are_independent_screens_not_the_rational_proof"] is True,
                  f"prime-screen interpretation m={m}")

    minor = cert["odd_maximal_minor_certificate"]
    row_ids = minor["minor_row_indices"]
    col_ids = minor["minor_column_indices"]
    audit.require(minor["minor_order"] == len(row_ids) == len(col_ids) == n - 5, f"minor order m={m}")
    audit.require(len(set(row_ids)) == len(row_ids) and len(set(col_ids)) == len(col_ids), f"minor index duplicates m={m}")
    audit.require(all(0 <= r < len(rows) for r in row_ids) and all(0 <= c < n for c in col_ids),
                  f"minor index range m={m}")
    rank2 = rank_restricted_minor_mod2(rows, row_ids, col_ids)
    audit.require(rank2 == n - 5 == minor["minor_rank_mod_2_recheck"], f"independent minor rank m={m}")
    audit.require(minor["minor_determinant_mod_2"] == 1, f"odd minor determinant m={m}")
    audit.note(f"m={m}: {n} vertices, {len(rows)} rows, degree {min(degrees)}--{max(degrees)}, odd minor rank {rank2}")


def audit_gap_rerun(audit: Audit, source: Path) -> None:
    obj = load_json(source / "PARENT_RERUNS/independent_gap_recompute.json")
    audit.require(obj["schema"] == "parent_independent_gap_recompute.v1", "gap-rerun schema")
    audit.require(obj["normalization"] == "A=(1/2) C D^(-1/2)", "gap-rerun normalization")
    audit.require(obj["raw_json_only"] is True and obj["all_checks_passed"] is True, "gap-rerun global pass")
    results = {r["m"]: r for r in obj["results"]}
    audit.require(set(results) == set(EXPECTED), "gap-rerun m set")
    for m, spec in EXPECTED.items():
        row = results[m]
        audit.require(row["pass"] is True and row["warning_messages"] == [], f"fresh gap pass/warnings m={m}")
        audit.require(row["vertex_count"] == spec["n"] and row["rectangle_count"] == spec["rows"], f"fresh gap size m={m}")
        audit.require(row["degree_min"] == spec["dmin"] and row["degree_max"] == spec["dmax"], f"fresh gap degrees m={m}")
        audit.require(math.isclose(row["stored_gap_sigma"], spec["stored_gap"], rel_tol=0, abs_tol=1e-15), f"rerun stored gap m={m}")
        audit.require(math.isclose(row["fresh_gap_sigma"], spec["fresh_gap"], rel_tol=0, abs_tol=1e-15), f"fresh gap value m={m}")
        audit.require(row["fresh_operator_residual"] < 6e-10, f"fresh operator residual m={m}")
        audit.require(row["fresh_constraint_residual"] < 1e-8, f"fresh constraint residual m={m}")
        audit.require(row["fresh_relative_sigma_error"] < 2e-5, f"fresh relative gap error m={m}")


def moment_values(m: int) -> tuple[int, int, int, int]:
    s0 = 2 * m + 1
    s2 = m * (m + 1) * (2 * m + 1) // 3
    s4 = m * (m + 1) * (2 * m + 1) * (3 * m * m + 3 * m - 1) // 15
    g = 2 * s0 * (s0 * s4 - s2 * s2)
    return s0, s2, s4, g


def audit_theorem_arithmetic(audit: Audit, source: Path) -> None:
    rerun = load_json(source / "PARENT_RERUNS/bounded_span_validation.json")
    audit.require(rerun["schema"] == "bounded_span_theorem_validation.v1", "theorem-rerun schema")
    audit.require(rerun["all_checks_passed"] is True, "theorem-rerun pass")
    audit.require(len(rerun["exact_moment_and_orthogonality_checks"]) == 10, "theorem-rerun moment count")
    audit.require(len(rerun["block_candidate_inequality_checks"]) == 4, "theorem-rerun block count")
    audit.require(all(x["equal"] and x["unweighted_mode_inner_products"] == [0] * 5
                      for x in rerun["exact_moment_and_orthogonality_checks"]), "theorem-rerun orthogonality")
    audit.require(all(x["inequality_chain_verified"] for x in rerun["block_candidate_inequality_checks"]),
                  "theorem-rerun inequality chain")

    for m in range(1, 11):
        s0, _, _, formula = moment_values(m)
        closed = 2 * (2 * m + 1) ** 3 * m * (m + 1) * (4 * m * m + 4 * m - 3) // 45
        brute = 0
        inner = [0] * 5
        for x in range(-m, m + 1):
            for y in range(-m, m + 1):
                for z in range(-m, m + 1):
                    g = x * x - y * y
                    k = [1, x, y, z, x * x + y * y + z * z]
                    brute += g * g
                    for j in range(5):
                        inner[j] += g * k[j]
        audit.require(formula == closed == brute, f"independent G_m identity m={m}")
        audit.require(inner == [0] * 5, f"independent orthogonality m={m}")
        audit.require(s0 ** 3 == (2 * m + 1) ** 3, f"vertex count identity m={m}")


def citation_keys(tex: str) -> set[str]:
    keys: set[str] = set()
    for group in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}", tex):
        keys.update(k.strip() for k in group.split(",") if k.strip())
    return keys


def bib_keys(bib: str) -> set[str]:
    return set(re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", bib))


def audit_manuscript(audit: Audit, root: Path) -> None:
    texts: dict[str, str] = {}
    for name in PUBLIC_FILES:
        path = root / name
        audit.require(path.is_file() and path.stat().st_size > 0, f"missing/empty public artifact {name}")
        text = path.read_text(encoding="utf-8")
        texts[name] = text
        audit.require("\x00" not in text, f"NUL in {name}")
        audit.require("<<<<<<<" not in text and ">>>>>>>" not in text and "=======" not in text, f"conflict marker in {name}")
        audit.require(not any(line.endswith((" ", "\t")) for line in text.splitlines()), f"trailing whitespace in {name}")
        audit.require(not re.search(r"/Users/|\\Users\\|\.lingtai|daemons/em-[0-9]+", text), f"private path/run id in {name}")

    tex = texts["main.tex"]
    plan = texts["PAPER_PLAN.md"]
    readme = texts["README.md"]
    matrix = texts["SOURCE_MATRIX.md"]
    zh = texts["SUMMARY.zh-CN.md"]
    report = texts["REPORT.md"]
    bib = texts["references.bib"]

    audit.require(tex.count("{") == tex.count("}"), "main.tex brace count")
    audit.require("\\begin{document}" in tex and "\\end{document}" in tex, "main.tex document markers")
    audit.require(("\\begin{theorem}" in tex or "\\begin{proposition}" in tex)
                  and "\\begin{proof}" in tex and "\\end{proof}" in tex,
                  "main theorem/proof structure")
    audit.require("A_m=\\frac12 C_mD_m^{-1/2}" in tex, "normalization missing from main.tex")
    audit.require("\\sum_{R\\in\\mathcal R_m}|p_R|^2|q_R|^2" in tex, "weighted-span theorem missing")
    audit.require("\\rho=\\min\\{1,c/\\sqrt{135D}\\}" in tex, "macroscopic threshold rho missing")
    audit.require("\\eta=\\min\\{1,c^2/(540D)\\}" in tex, "macroscopic threshold eta missing")
    audit.require("literature-boundary" in tex.lower(), "literature-boundary section missing")
    audit.require("unrestricted nonlocal" in tex.lower() and "remains open" in tex.lower(), "open boundary missing")
    audit.require("finite evidence only" in tex, "finite-evidence boundary missing")

    english_disclosure = (
        "Every sentence of manuscript prose and the complete manuscript structure were "
        "generated by LingTai AI under Runyuan Wang's direction and authorization."
    )
    chinese_disclosure = (
        "本文全部正文文字与整体结构均由 LingTai AI 在王润圆的指导和授权下生成"
    )
    english_playful_disclaimer = (
        "This work is shared purely for personal enjoyment and playful exploration with AI. "
        "If anything is incorrect, please excuse it; criticism and corrections are warmly welcomed."
    )
    chinese_playful_disclaimer = (
        "这纯粹是自娱自乐，和AI一起玩耍，如果有不对之处，请见谅，请批评指正"
    )
    normalized_tex = re.sub(r"\s+", " ", tex)
    audit.require(english_disclosure in normalized_tex,
                  "exact English all-AI-writing disclosure missing")
    audit.require("AI production is not evidence for any claim" in normalized_tex,
                  "English evidence boundary missing from disclosure")
    audit.require("Runyuan Wang retains final scientific judgment and responsibility" in normalized_tex,
                  "English responsibility boundary missing")
    audit.require(english_playful_disclaimer in normalized_tex,
                  "English personal-enjoyment and correction disclaimer missing")
    audit.require(chinese_disclosure in zh,
                  "exact Chinese all-AI-writing disclosure missing")
    audit.require("AI 生成身份本身不构成任何结论的证据" in zh,
                  "Chinese evidence boundary missing")
    audit.require("王润圆保留最终科学判断" in zh and "承担最终责任" in zh,
                  "Chinese responsibility boundary missing")
    audit.require(chinese_playful_disclaimer in zh,
                  "Chinese personal-enjoyment disclaimer missing from summary")
    audit.require(chinese_playful_disclaimer in readme,
                  "Chinese personal-enjoyment disclaimer missing from README")
    audit.require("English and Chinese disclosures state exactly" in plan,
                  "AI disclosure acceptance gate missing from PAPER_PLAN.md")

    for token in ("[R]", "[E]", "[N]", "[O]"):
        audit.require(token in plan and token in matrix and token in zh, f"evidence token {token} absent")
    audit.require("SECOND_PAPER_" in report, "final verdict missing from REPORT.md")
    verdicts = re.findall(r"SECOND_PAPER_(?:DRAFT_READY_FOR_PARENT_REVIEW|DRAFT_NEEDS_MATERIAL_CORRECTION|BLOCKED)", report)
    audit.require(len(verdicts) == 1, "REPORT.md must contain exactly one authorized verdict")
    audit.require(report.rstrip().endswith(verdicts[0]), "REPORT.md verdict must be the final line")

    cites = citation_keys(tex)
    entries = bib_keys(bib)
    audit.require(cites <= CITATION_ALLOWLIST, f"citation set mismatch: {sorted(cites)}")
    audit.require(CITATION_ALLOWLIST <= entries, f"bibliography set mismatch: {sorted(entries)}")
    audit.require(cites <= entries, "unresolved citation keys")
    audit.require(all(key in matrix for key in cites), "cited key absent from SOURCE_MATRIX.md")

    combined = "\n".join(texts.values()).lower()
    positive_overclaims = [
        r"this (?:paper|result) (?:settles|solves|disproves) (?:the )?(?:yuan|unrestricted)",
        r"liminf[^\n]{0,40}>\s*0[^\n]{0,80}(?:therefore|hence|proves)",
        r"(?:this|our) (?:paper|result|method|bound)[^\n]{0,40}(?:is|gives|achieves|sets)[^\n]{0,30}(?:first[- ]ever|world[- ]first|current[- ]best|state[- ]of[- ]the[- ]art)",
        r"our (?:search|attacks?) (?:is|are|was|were) exhaustive",
        r"no (?:prior|previous) work exists",
    ]
    for pattern in positive_overclaims:
        audit.require(re.search(pattern, combined) is None, f"positive overclaim pattern: {pattern}")

    for value in ("0.137567327104", "0.108836005898", "0.107304440463", "0.104461520420"):
        audit.require(value in tex and value in zh, f"finite gap value {value} not synchronized")
    for row in ("125", "343", "729", "1331", "188", "515", "1094", "1997"):
        audit.require(row in tex and row in zh, f"finite count {row} not synchronized")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True,
                        help="accepted collision-rectangle experiment root")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent,
                        help="second-paper draft root")
    args = parser.parse_args()
    source = args.source.resolve()
    root = args.root.resolve()
    audit = Audit()
    audit.require(source.is_dir(), "source directory does not exist")
    audit.require(root.is_dir(), "draft root does not exist")

    for m, spec in EXPECTED.items():
        audit_candidate_and_certificate(audit, source, m, spec)
    audit_gap_rerun(audit, source)
    audit_theorem_arithmetic(audit, source)
    audit_manuscript(audit, root)

    print(f"PASS {audit.checks} deterministic assertions; {len(audit.notes)} candidate notes")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FAIL {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
