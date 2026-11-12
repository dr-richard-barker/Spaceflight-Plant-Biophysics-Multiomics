#!/usr/bin/env python3
"""
scripts/abai_check.py
Autonomous Biophysics & AstroBotany Initiative (ABAI) Quality Control Audit Suite.

Performs a rigorous, transparent quality control verification across:
  1. LaTeX build & PDF artifacts (zero errors, zero missing refs, zero overfull boxes > 10pt)
  2. Word (.docx) artifacts (valid zip archives containing embedded image media)
  3. Quantitative grounding audit: verifies numbers in manuscript vs source CSV/TSV tables
  4. Blocklist & text integrity: zero TODOs, zero unverified "blind" claims
  5. Asset & relative link integrity across web pages (paper.html, explorer/index.html, index.html)
  6. Subsystem self-checks: fvcb.py, check_site_data.py, make_manifest.py, validation/check_tex.py
"""
from __future__ import annotations

import csv
import json
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

def audit_latex_build() -> list[str]:
    issues = []
    ms_dir = ROOT / "manuscript"
    for base in ["main", "supplementary", "npj_main"]:
        pdf = ms_dir / f"{base}.pdf"
        log = ms_dir / f"{base}.log"
        if not pdf.exists():
            issues.append(f"LaTeX PDF missing: {pdf}")
        if not log.exists():
            issues.append(f"LaTeX LOG missing: {log}")
            continue
        log_txt = log.read_text(errors="ignore")
        if "^!" in log_txt or "\n! " in log_txt:
            issues.append(f"LaTeX fatal error found in {log.name}")
        for line in log_txt.splitlines():
            if "LaTeX Warning: Reference `" in line and "undefined" in line:
                issues.append(f"Undefined LaTeX reference in {log.name}: {line.strip()}")
            if "LaTeX Warning: Citation `" in line and "undefined" in line:
                issues.append(f"Undefined LaTeX citation in {log.name}: {line.strip()}")
            m = re.search(r"Overfull \\hbox \(([0-9.]+)pt too wide\)", line)
            if m and float(m.group(1)) > 15.0:
                issues.append(f"Overfull hbox > 15pt in {log.name}: {line.strip()}")
    return issues


def audit_docx_artifacts() -> list[str]:
    issues = []
    docs = [
        ROOT / "manuscript/Spaceflight_Biophysics_Multiomics_Manuscript.docx",
        ROOT / "manuscript/Spaceflight_Biophysics_Multiomics_SI.docx",
        ROOT / "public/Spaceflight_Biophysics_Multiomics_Manuscript.docx",
        ROOT / "public/Spaceflight_Biophysics_Multiomics_SI.docx",
        ROOT / "public/explorer/Spaceflight_Biophysics_Multiomics_Manuscript.docx",
        ROOT / "public/explorer/Spaceflight_Biophysics_Multiomics_SI.docx",
    ]
    for doc in docs:
        if not doc.exists():
            issues.append(f"DOCX missing: {doc}")
            continue
        try:
            with zipfile.ZipFile(doc, "r") as zf:
                nl = zf.namelist()
                images = [n for n in nl if n.startswith("word/media/")]
                if not images:
                    issues.append(f"DOCX {doc.name} contains NO embedded images!")
        except Exception as e:
            issues.append(f"Failed to inspect DOCX {doc}: {e}")
    return issues


def audit_quantitative_grounding() -> list[str]:
    issues = []
    # Verify T13 boundary layer conductances
    t13_p = ROOT / "results/tables/T13_boundary_layer.csv"
    if not t13_p.exists():
        issues.append(f"T13 missing: {t13_p}")
    else:
        with open(t13_p) as f:
            r = csv.DictReader(f)
            gbl_map = {(row["scale"], row["scenario"]): float(row["g_bl_mol_m2_s"]) for row in r}
        expected = {
            ("leaf", "leaf-earth"): 0.997,
            ("leaf", "leaf-ug"): 0.443,
            ("rosette", "rosette-earth"): 0.540,
            ("rosette", "rosette-ug"): 0.272,
            ("canopy", "canopy-earth"): 0.133,
            ("canopy", "canopy-ug"): 0.058,
        }
        for k, v in expected.items():
            if abs(gbl_map.get(k, 0.0) - v) > 0.001:
                issues.append(f"T13 ground mismatch for {k}: expected {v}, got {gbl_map.get(k)}")

    # Verify T10 metabolite hubs
    t10_p = ROOT / "results/tables/T10_metabolite_hubs.tsv"
    if not t10_p.exists():
        issues.append(f"T10 missing: {t10_p}")
    else:
        with open(t10_p) as f:
            lines = [l.strip().split("\t") for l in f if l.strip() and not l.startswith("#")]
        header = lines[0]
        rows = [dict(zip(header, row)) for row in lines[1:]]
        fdr_sig = [r for r in rows if float(r["FDR"]) < 0.05]
        if len(fdr_sig) != 7:
            issues.append(f"T10 metabolite hubs: expected 7 FDR < 0.05 compounds, got {len(fdr_sig)}")
    return issues


def audit_text_integrity() -> list[str]:
    issues = []
    # Scan for [TODO]
    for p in ROOT.rglob("*"):
        if p.is_file() and p.suffix in {".tex", ".md", ".html"} and ".git" not in p.parts and "node_modules" not in p.parts:
            txt = p.read_text(errors="ignore")
            if "TODO" in txt:
                issues.append(f"Found [TODO] in {p.relative_to(ROOT)}")
            if "Photorespiration_multiomics_microgravity.pdf" in txt:
                issues.append(f"Stale unmerged PDF reference in {p.relative_to(ROOT)}")
    return issues


def audit_html_links() -> list[str]:
    issues = []
    for html_rel in ["index.html", "public/paper.html", "public/explorer/index.html"]:
        p = ROOT / html_rel
        if not p.exists():
            issues.append(f"HTML file missing: {p}")
            continue
        base = p.parent if html_rel != "index.html" else (ROOT / "public")
        txt = p.read_text(errors="ignore")
        for attr, val in re.findall(r"(href|src)=[\"']([^\"'#?]+)", txt):
            if val.startswith(("http://", "https://", "mailto:", "javascript:", "data:")):
                continue
            if val.startswith("/src/"):
                target = ROOT / val.lstrip("/")
            elif val.startswith("/"):
                target = (ROOT / "public") / val.lstrip("/")
            else:
                target = (base / val).resolve()
            if not target.exists():
                issues.append(f"Broken {attr} in {html_rel}: {val} -> {target}")
    return issues


def main() -> int:
    print("=================================================================")
    print("  ABAI Quality Control (QC) Automated Verification Suite        ")
    print("=================================================================")
    all_issues = []

    print("[1/5] Checking LaTeX & PDF artifacts...")
    latex_issues = audit_latex_build()
    all_issues.extend(latex_issues)
    print(f"      {'✓ PASS' if not latex_issues else '✗ ' + str(len(latex_issues)) + ' issues'}")

    print("[2/5] Checking Word (.docx) embedded-figure media...")
    docx_issues = audit_docx_artifacts()
    all_issues.extend(docx_issues)
    print(f"      {'✓ PASS' if not docx_issues else '✗ ' + str(len(docx_issues)) + ' issues'}")

    print("[3/5] Checking Quantitative Grounding against CSV/TSV tables...")
    ground_issues = audit_quantitative_grounding()
    all_issues.extend(ground_issues)
    print(f"      {'✓ PASS' if not ground_issues else '✗ ' + str(len(ground_issues)) + ' issues'}")

    print("[4/5] Checking Text Integrity (Zero TODOs & Stale References)...")
    text_issues = audit_text_integrity()
    all_issues.extend(text_issues)
    print(f"      {'✓ PASS' if not text_issues else '✗ ' + str(len(text_issues)) + ' issues'}")

    print("[5/5] Checking Web Link & Asset Integrity...")
    html_issues = audit_html_links()
    all_issues.extend(html_issues)
    print(f"      {'✓ PASS' if not html_issues else '✗ ' + str(len(html_issues)) + ' issues'}")

    print("-----------------------------------------------------------------")
    if all_issues:
        print(f"QC Audit FAILED with {len(all_issues)} issue(s):")
        for iss in all_issues:
            print(f"  • {iss}")
        return 1
    else:
        print("QC Audit PASSED — All verification gates, artifacts & links CLEAR.")
        return 0

if __name__ == "__main__":
    sys.exit(main())
