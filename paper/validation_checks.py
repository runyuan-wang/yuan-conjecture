#!/usr/bin/env python3
"""Deterministic editorial checks for the bounded draft package.

Usage (from draft/):
  python3 validation_checks.py --project .. --source "$SOURCE_ROOT"
The source root is read-only and is selected by the frozen evidence matrices.
This script does not execute any source-package research computation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

EXPECTED_MATRIX_HASHES = {
    "CLAIM_EVIDENCE_MATRIX.md": "92959c933457c120cd9b5bd25bf41cdae4ff84737b3efa1927cdfc44ec05c683",
    "INTERNAL_EVIDENCE_MANIFEST.md": "35dd8c394a1677954d517e57507cb9a8046f62ef1f2b70aa4db6cd06e39e6fe6",
}
PUBLIC_FILES = ("main.tex", "references.bib", "SUMMARY.zh-CN.md", "paper.html")
ASSET_PREFIX = "../assets/yuanjiang_kakeya_recollision_assets_2026-08-11/unpacked/"
REQUIRED_ASSETS = {
    "eight_atom_no_go_certificate_2026-08-11.json",
    "eight_atom_cycle_search.py",
    "verify_eight_atom_certificate.py",
    "strict_time_c12_worldline_certificate.json",
    "verify_strict_time_c12_worldlines.py",
    "dual_core_residual_audit.py",
    "entropy_prover_residual_target_edge_certificate_2026-08-13.json",
    "weighted_cycle_cutting.py",
    "psitip_cycle_audit.py",
    "dual_core_tertiary_certificate_2026-08-12.json",
    "protein_tertiary_fold_certificate_2026-08-12.json",
    "heterogeneous_tertiary_fold_certificate_2026-08-12.json",
    "proof_route_matrix_2026-08-13.md",
    "SHA256SUMS_2026-08-13.txt",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def need(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def contains_all(path: Path, needles: list[str], label: str) -> None:
    text = path.read_text(encoding="utf-8")
    missing = [x for x in needles if x not in text]
    need(not missing, f"{label}: missing tokens {missing}")


def json_atoms(path: Path) -> set[str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    out: set[str] = set()
    def walk(value):
        if isinstance(value, dict):
            for key, item in value.items():
                out.add(str(key))
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)
        else:
            out.add(str(value))
    walk(data)
    return out


class AuditHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.hrefs: list[str] = []
        self.external_resources: list[tuple[str, str]] = []
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.add(a["id"])
        if "href" in a:
            self.hrefs.append(a["href"])
        if tag in {"script", "img", "iframe", "object", "embed"} and "src" in a:
            self.external_resources.append((tag, a["src"]))
        if tag == "link" and a.get("rel") in {"stylesheet", "preload"}:
            self.external_resources.append((tag, a.get("href", "")))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    args = parser.parse_args()
    draft = Path.cwd().resolve()
    project = args.project.resolve()
    source = args.source.resolve()

    required = ["main.tex", "references.bib", "SUMMARY.zh-CN.md", "paper.html", "main.pdf", "main.log"]
    need(all((draft / name).is_file() for name in required), "required draft/build file absent")
    print("PASS package-presence: manuscript, bibliography, Chinese summary, HTML, PDF, and final log")

    for name, expected in EXPECTED_MATRIX_HASHES.items():
        need(sha256(project / name) == expected, f"frozen hash mismatch: {name}")
    print("PASS frozen-matrix-hashes: 2/2 match supplied SHA-256 values")

    # Every internal-manifest hash row is checked against the selected source root.
    manifest_text = (project / "INTERNAL_EVIDENCE_MANIFEST.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| `([^`]+)` \| `([0-9a-f]{64})` \|", manifest_text, flags=re.M)
    need(rows, "no source hash rows parsed from internal evidence manifest")
    for filename, expected in rows:
        need((source / filename).is_file(), f"internal-manifest source missing: {filename}")
        need(sha256(source / filename) == expected, f"internal-manifest hash mismatch: {filename}")
    print(f"PASS decisive-source-hashes: {len(rows)}/{len(rows)} internal-manifest rows match")

    # Check all entries in the public payload manifest without running computations.
    sha_manifest = source / "SHA256SUMS_2026-08-13.txt"
    sha_rows = []
    for line in sha_manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-f]{64})\s+\*?(.+)", line)
        need(m is not None, f"unparsed SHA manifest line: {line!r}")
        sha_rows.append((m.group(2), m.group(1)))
    need(len(sha_rows) == 38, f"expected 38 payload hashes, got {len(sha_rows)}")
    for filename, expected in sha_rows:
        need((source / filename).is_file(), f"payload missing: {filename}")
        need(sha256(source / filename) == expected, f"payload hash mismatch: {filename}")
    print("PASS source-payload-hashes: 38/38 manifest entries match")

    tex = (draft / "main.tex").read_text(encoding="utf-8")
    bib = (draft / "references.bib").read_text(encoding="utf-8")
    lit = (project / "literature" / "VERIFIED_LITERATURE_MATRIX.md").read_text(encoding="utf-8")
    cite_keys: set[str] = set()
    for match in re.finditer(r"\\cite(?:p|t)?(?:\[[^\]]*\])?(?:\[[^\]]*\])?\{([^}]+)\}", tex):
        cite_keys.update(k.strip() for k in match.group(1).split(","))
    bib_keys = set(re.findall(r"^@[A-Za-z]+\{([^,]+),", bib, flags=re.M))
    need(cite_keys == bib_keys, f"citation/bibliography mismatch: cited={sorted(cite_keys)} bib={sorted(bib_keys)}")
    missing_retained = sorted(k for k in cite_keys if k not in lit)
    need(not missing_retained, f"citation key absent from retained verified rows: {missing_retained}")
    need(len(cite_keys) == 6, f"unexpected citation-key count: {len(cite_keys)}")
    print("PASS citations: 6 cited keys = 6 bibliography keys; every key occurs in a retained verified row")

    log = (draft / "main.log").read_text(encoding="utf-8", errors="replace")
    final_bad = [
        r"Citation .* undefined", r"There were undefined citations",
        r"Reference .* undefined", r"There were undefined references",
        r"Label\(s\) may have changed", r"Rerun to get cross-references right",
    ]
    hits = [pat for pat in final_bad if re.search(pat, log)]
    need(not hits, f"unresolved-reference patterns in final LaTeX log: {hits}")
    print("PASS final-reference-log: no unresolved citation/reference or rerun marker")

    public = {name: (draft / name).read_text(encoding="utf-8") for name in PUBLIC_FILES}
    for name, text in public.items():
        need("/Users/" not in text, f"private absolute path in {name}")
        need(".lingtai" not in text.lower(), f"private framework path in {name}")
        need(not re.search(r"\bem-[0-9a-f]{4,}\b", text, flags=re.I), f"internal worker ID in {name}")
        need("Boltzmann" not in text and "boltzmann" not in text, f"unsupported Boltzmann material in {name}")
    print("PASS public-privacy-and-omission: no private path/worker ID; Boltzmann--Grad layer absent")

    disclosure_requirements = {
        "main.tex": [
            "This manuscript was written and prepared using",
            "https://github.com/Lingtai-AI/lingtai",
            "LingTai AI assisted with evidence organization",
            "Runyuan Wang retains sole authorship, scientific judgment, and responsibility for the final content.",
        ],
        "paper.html": [
            "This manuscript was written and prepared using",
            "https://github.com/Lingtai-AI/lingtai",
            "LingTai AI assisted with evidence organization",
            "Runyuan Wang retains sole authorship, scientific judgment, and responsibility for the final content.",
        ],
        "SUMMARY.zh-CN.md": [
            "本文使用 [LingTai AI（灵台）]",
            "协助证据整理、来源—主张核验、草稿写作、一致性检查和文档构建",
            "Runyuan Wang 保留唯一作者身份、科学判断与最终内容责任。",
        ],
    }
    for filename, tokens in disclosure_requirements.items():
        missing = [token for token in tokens if token not in public[filename]]
        need(not missing, f"LingTai tool-use disclosure absent from {filename}: {missing}")
    print("PASS LingTai-tool-use-disclosure: explicit assistance and sole-human-authorship language agree across TeX, HTML, and Chinese summary")

    # Ban common affirmative overclaims while permitting explicit negations/non-claims.
    combined = "\n".join(public.values())
    prohibited = {
        "affirmative proof": r"\b(?:we|this (?:paper|work|study)) (?:here )?(?:prove|proves|establish|establishes)\b",
        "affirmative disproof": r"\b(?:we|this (?:paper|work|study)) (?:here )?(?:disprove|disproves|refute|refutes) (?:the )?(?:four-dimensional )?Kakeya\b",
        "new bound claim": r"\b(?:we (?:obtain|establish|prove)|this (?:paper|work) (?:obtains|establishes|proves)) (?:a )?new (?:Kakeya )?bound\b",
        "novelty claim": r"\b(?:our|a) novel (?:proof|method|result|bound|theorem)\b",
        # `current-best` may occur only inside an explicit local negation/non-claim.
        "current-best claim": r"(?<!do not describe that number as )\b(?:current best|current-best|best known)\b",
        "exhaustive global search": r"\b(?:we exhaust|exhaustive search of all|exhausts all)\b",
        "reproduction claim": r"\b(?:we|this (?:paper|work)) (?:independently )?(?:reproduced|reran|re-ran)\b",
    }
    bad = {label: re.findall(rx, combined, flags=re.I) for label, rx in prohibited.items() if re.search(rx, combined, flags=re.I)}
    need(not bad, f"forbidden affirmative-claim lint: {bad}")
    required_scope_patterns = {
        "does not prove or disprove": r"does not prove or disprove",
        "does not independently rerun": r"does not independently rerun",
        "PSITIP/GLOP lacks exact rational dual certificate": r"PSITIP/GLOP output does not include an exact rational dual certificate",
        "beam nulls are not no-go theorems": r"beam nulls are not a no-go theorem",
        "manifest/static audit is not reproduction": r"Manifest/hash and static source audit do not constitute independent reproduction",
    }
    for label, pattern in required_scope_patterns.items():
        need(re.search(pattern, tex, flags=re.I) is not None, f"mandatory boundary absent from main.tex: {label}")
    need("Runyuan Wang" in tex and "\\affil" not in tex and "affiliation" not in tex.lower(), "author/affiliation constraint failed")
    print("PASS claim lint: no affirmative proof/disproof/new-bound/novelty/current-best/exhaustiveness/reproduction claim; mandatory boundaries present")

    # Source correspondence: exact values are required in both source material and manuscript.
    c1 = json_atoms(source / "eight_atom_no_go_certificate_2026-08-11.json")
    for atom in ["1725", "313762", "16", "6", "6885", "110160", "8931424", "13", "8"]:
        need(atom in c1, f"C1 source atom absent: {atom}")
    contains_all(draft / "main.tex", ["1,725", "313,762", "16", "6,885", "110,160", "8,931,424", "1{,}000{,}000{,}007", "10^8", "1{,}000{,}000{,}009", "anchors $0$ and $5$"], "C1 manuscript")

    c2_path = source / "strict_time_c12_worldline_certificate.json"
    c2 = json_atoms(c2_path)
    c2_raw = c2_path.read_text(encoding="utf-8")
    for atom in ["201/100", "99/100", "2017/55600", "529/12200", "7/3000", "43/33000"]:
        need(atom in c2, f"C2 source atom absent: {atom}")
    need("1/25" in c2_raw, "C2 source threshold absent: 1/25")
    contains_all(draft / "main.tex", ["(0,1,2,3,201/100,99/100)", "Q_0<Q_5<Q_1<Q_2<Q_4<Q_3", "2017}{55600", "529}{12200", "7}{3000", "43}{33000", "\\delta<1/25"], "C2 manuscript")

    c3_text = (source / "dual_core_residual_audit.py").read_text(encoding="utf-8") + (source / "entropy_prover_residual_target_edge_certificate_2026-08-13.json").read_text(encoding="utf-8")
    for token in ["27", "20", "18", "2", "6", "7", "(1, -1)"]:
        need(token in c3_text, f"C3 source token absent: {token}")
    contains_all(draft / "main.tex", ["27 rows", "rank}(A)=18", "20-18=2", "X_0,\\ C_4,\\ D_5,\\ Y_0,\\ Y_1,\\ Y_2,\\ X_8", "0,1,2,8,9", "atoms 6 and 7", "=(1,-1)"], "C3 manuscript")

    # C4/N quantitative screen and regression values: exact source files and manuscript table/prose.
    source_to_tex_tokens = {
        "protein_tertiary_fold_certificate_2026-08-12.json": [
            ("9635", "9,635"), ("6076559", "6,076,559"), ("32", "width-32"),
            ("826472", "826,472"), ("8", "$8/10$"),
        ],
        "heterogeneous_tertiary_fold_certificate_2026-08-12.json": [
            ("14931", "14,931"), ("3423723", "3,423,723"), ("67348", "67,348"),
            ("1649587", "1,649,587"), ("8", "$8/11$"),
        ],
        "dual_core_tertiary_certificate_2026-08-12.json": [
            ("15936", "15,936"), ("3332936", "3,332,936"), ("108", "108 folds"),
            ("92983", "92,983"), ("561380", "561,380"),
        ],
        "entropy_prover_residual_target_edge_certificate_2026-08-13.json": [
            ("1.1.7", "PSITIP 1.1.7"), ("9.15.6755", "OR-Tools 9.15.6755"),
            ("256", "width 256"), ("800", "800 structures"), ("16", "width-16"),
            ("857064", "857,064"),
        ],
    }
    for filename, pairs in source_to_tex_tokens.items():
        raw = (source / filename).read_text(encoding="utf-8")
        for source_token, tex_token in pairs:
            need(source_token in raw, f"C4/N source token absent from {filename}: {source_token}")
            need(tex_token in tex, f"C4/N manuscript token absent for {filename}: {tex_token}")
    sampling_tokens = ["24/160", "40/640", "16/64"]
    sampling_source = (source / "protein_tertiary_fold_search.py").read_text(encoding="utf-8")
    for token in [
        'deterministic_sample(pairs, 24, "pair:" + salt)',
        'deterministic_sample(triples, 40, "triple:" + salt)',
        'deterministic_sample(full, 16, "full:" + salt)',
        'deterministic_sample(pairs, 16, "anti-pair:" + salt)',
        'deterministic_sample(triples, 24, "anti-triple:" + salt)',
        'deterministic_sample(full, 8, "anti-full:" + salt)',
    ]:
        need(token in sampling_source, f"open-path sampling source token absent: {token}")
    for label, text in [("TeX", tex), ("HTML", public["paper.html"]), ("Chinese summary", public["SUMMARY.zh-CN.md"])]:
        for token in sampling_tokens:
            need(token in text, f"open-path sampling disclosure absent from {label}: {token}")
    print("PASS screen-regression-correspondence: named source certificates match all reported C4/N counts, widths, ranks, versions, and open-path sampling fractions")

    # Every published artifact link is allowlisted and maps to an actual source payload.
    for name in REQUIRED_ASSETS:
        need((source / name).is_file(), f"required linked source absent: {name}")
        need((ASSET_PREFIX + name) in public["paper.html"], f"HTML missing public asset link: {name}")
    tex_linked = set(re.findall(r"\\assetpath\s+([^}\s]+)", tex))
    need(tex_linked == REQUIRED_ASSETS, f"TeX asset link set mismatch: {sorted(tex_linked ^ REQUIRED_ASSETS)}")
    print("PASS exact-correspondence: C1/C2/C3 values anchored in named frozen artifacts; all public asset links allowlisted")

    # Cross-format minimum parity for quantitative conclusions and boundaries.
    parity_patterns = {
        "1,725": r"1,725", "313,762": r"313,762", "110,160": r"110,160", "8,931,424": r"8,931,424",
        "2017/55600": r"2017(?:/|\}\{)55600", "529/12200": r"529(?:/|\}\{)12200",
        "7/3000": r"7(?:/|\}\{)3000", "43/33000": r"43(?:/|\}\{)33000",
        "27": r"\b27\b", "20": r"\b20\b", "rank": r"rank", "18": r"\b18\b", "nullity": r"nullity|零度", "2": r"\b2\b", "PSITIP": r"PSITIP", "GLOP": r"GLOP",
    }
    for filename in ("main.tex", "SUMMARY.zh-CN.md", "paper.html"):
        text = public[filename]
        missing = [label for label, pattern in parity_patterns.items() if not re.search(pattern, text, flags=re.I)]
        need(not missing, f"cross-format parity tokens missing in {filename}: {missing}")
    print("PASS cross-format claim parity: key exact counts/margins/rank and PSITIP/GLOP limits occur in TeX, Chinese summary, and HTML")

    # HTML structural/static checks.
    html = public["paper.html"]
    hp = AuditHTMLParser()
    hp.feed(html)
    need(not hp.external_resources, f"HTML external resources violate self-contained requirement: {hp.external_resources}")
    broken_fragments = [h for h in hp.hrefs if h.startswith("#") and h[1:] not in hp.ids]
    need(not broken_fragments, f"broken HTML fragments: {broken_fragments}")
    local_assets = [h for h in hp.hrefs if not re.match(r"^(?:https?://|mailto:|#)", h)]
    need(set(local_assets) == {ASSET_PREFIX + x for x in REQUIRED_ASSETS}, "HTML local-link set differs from allowlist")
    need("<style>" in html and "<script" not in html.lower(), "HTML CSS/script self-contained check failed")
    print(f"PASS HTML-static: {len(hp.ids)} IDs, {len(hp.hrefs)} links, {len(local_assets)} mapped local assets, zero external resources/broken fragments")

    # A git-diff-check equivalent: no trailing spaces/tabs, conflict markers, or missing final newline.
    whitespace_bad = []
    for name in PUBLIC_FILES + ("validation_checks.py",):
        path = draft / name
        if not path.exists():
            continue
        data = path.read_text(encoding="utf-8")
        for number, line in enumerate(data.splitlines(), 1):
            if line.endswith((" ", "\t")):
                whitespace_bad.append(f"{name}:{number}: trailing whitespace")
            if line.startswith(("<<<<<<<", "=======", ">>>>>>>")):
                whitespace_bad.append(f"{name}:{number}: conflict marker")
        if data and not data.endswith("\n"):
            whitespace_bad.append(f"{name}: missing final newline")
    need(not whitespace_bad, f"draft-text whitespace/conflict check: {whitespace_bad}")
    print("PASS draft-text-diff-check: no trailing whitespace, conflict marker, or missing final newline")

    print("VALIDATION_CHECKS_PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"VALIDATION_CHECKS_FAIL: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
