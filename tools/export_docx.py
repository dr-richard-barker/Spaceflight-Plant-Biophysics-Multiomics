#!/usr/bin/env python3
"""
Generate a clean Word (.docx) manuscript from MANUSCRIPT.md using pure Python
standard library (zipfile + XML), and sync compiled PDF/DOCX artifacts into
public/ and public/explorer/.
"""
from __future__ import annotations

import html
import re
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MS_DIR = ROOT / "manuscript"
PUB_DIR = ROOT / "public"
EXP_DIR = PUB_DIR / "explorer"


def _clean_md_text(text: str) -> str:
    # Strip simple LaTeX math formatting into readable Unicode for Word
    replacements = [
        (r"\mu\mathrm{mol\,m^{-2}\,s^{-1}}", "µmol m⁻² s⁻¹"),
        (r"\mu\mathrm{mol\,h^{-1}\,cm^{-2}}", "µmol h⁻¹ cm⁻²"),
        (r"\mathrm{mol\,m^{-2}\,s^{-1}}", "mol m⁻² s⁻¹"),
        (r"\mathrm{cm\,s^{-1}}", "cm s⁻¹"),
        (r"\mathrm{m^2\,s^{-1}}", "m² s⁻¹"),
        (r"\mathrm{g\,mol^{-1}}", "g mol⁻¹"),
        (r"\mathrm{CO_2}", "CO₂"),
        (r"\mathrm{O_2}", "O₂"),
        (r"\mathrm{H_2O}", "H₂O"),
        (r"\mathrm{C_3}", "C₃"),
        (r"\mathrm{C_2}", "C₂"),
        (r"\mathrm{Re}", "Re"),
        (r"\mathrm{Ra}", "Ra"),
        (r"\mathrm{Pr}", "Pr"),
        (r"\mathrm{Nu}", "Nu"),
        (r"\mathrm{St}", "St"),
        (r"\mathrm{Sh}", "Sh"),
        (r"\mathrm{FDR}", "FDR"),
        (r"\mathrm{ppm}", "ppm"),
        (r"\mathrm{mm}", "mm"),
        (r"\mathrm{cm^2}", "cm²"),
        (r"\mathrm{cm}", "cm"),
        (r"\mathrm{mL}", "mL"),
        (r"\mathrm{min}", "min"),
        (r"\mathrm{ms}", "ms"),
        (r"\mathrm{d}", "d"),
        (r"\mathrm{h}", "h"),
        (r"\mathrm{s}", "s"),
        (r"\mu g", "µg"),
        (r"\varphi", "φ"),
        (r"\Gamma^*", "Γ*"),
        (r"\Gamma", "Γ"),
        (r"\Delta", "Δ"),
        (r"\delta", "δ"),
        (r"\rho", "ρ"),
        (r"\propto", "∝"),
        (r"\approx", "≈"),
        (r"\times", "×"),
        (r"\le", "≤"),
        (r"\ge", "≥"),
        (r"\to", "→"),
        (r"\xrightarrow{\;g_{bl}\;}", " →(g_bl)→ "),
        (r"\xrightarrow{\;g_s\;}", " →(g_s)→ "),
        (r"\xrightarrow{\;g_m\;}", " →(g_m)→ "),
        (r"\tfrac{1}{2}", "1/2"),
        (r"\mathrm{erfc}", "erfc"),
        (r"\log_2\mathrm{FC}", "log₂FC"),
        (r"\log_2", "log₂"),
        (r"\text{–}", "–"),
        (r"\,", " "),
        (r"\;", " "),
    ]
    out = text
    for k, v in replacements:
        out = out.replace(k, v)
    out = re.sub(r"\$\$([^$]+)\$\$", r"\1", out)
    out = re.sub(r"\$([^$]+)\$", r"\1", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"\1", out)
    out = re.sub(r"\*([^*]+)\*", r"\1", out)
    out = re.sub(r"`([^`]+)`", r"\1", out)
    return out.strip()


def _w_p(text: str, style: str | None = None, bold: bool = False) -> str:
    esc = html.escape(_clean_md_text(text))
    pPr = f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>' if style else ""
    rPr = "<w:rPr><w:b/></w:rPr>" if bold else ""
    return f"<w:p>{pPr}<w:r>{rPr}<w:t xml:space=\"preserve\">{esc}</w:t></w:r></w:p>"


def _w_table(rows: list[list[str]]) -> str:
    xml = [
        "<w:tbl>",
        "<w:tblPr><w:tblW w:w=\"5000\" w:type=\"pct\"/>"
        "<w:tblBorders>"
        "<w:top w:val=\"single\" w:sz=\"6\" w:space=\"0\" w:color=\"333333\"/>"
        "<w:bottom w:val=\"single\" w:sz=\"6\" w:space=\"0\" w:color=\"333333\"/>"
        "<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"CCCCCC\"/>"
        "</w:tblBorders></w:tblPr>",
    ]
    for r_idx, row in enumerate(rows):
        xml.append("<w:tr>")
        for cell in row:
            c_txt = html.escape(_clean_md_text(cell))
            rPr = "<w:rPr><w:b/><w:sz w:val=\"18\"/></w:rPr>" if r_idx == 0 else "<w:rPr><w:sz w:val=\"18\"/></w:rPr>"
            xml.append(
                f"<w:tc><w:p><w:r>{rPr}<w:t xml:space=\"preserve\">{c_txt}</w:t></w:r></w:p></w:tc>"
            )
        xml.append("</w:tr>")
    xml.append("</w:tbl>")
    return "".join(xml)


def build_docx(md_path: Path, docx_path: Path) -> None:
    lines = md_path.read_text().splitlines()
    body_parts: list[str] = []
    table_buf: list[list[str]] = []

    def flush_table() -> None:
        nonlocal table_buf
        if table_buf:
            body_parts.append(_w_table(table_buf))
            body_parts.append(_w_p(""))
            table_buf = []

    for line in lines:
        s = line.strip()
        if s.startswith("|") and s.endswith("|"):
            cells = [c.strip() for c in s.strip("|").split("|")]
            if all(set(c) <= set(":-") for c in cells):
                continue
            table_buf.append(cells)
            continue
        else:
            flush_table()

        if not s or s == "---":
            continue
        if s.startswith("# "):
            body_parts.append(_w_p(s[2:], style="Title", bold=True))
        elif s.startswith("## "):
            body_parts.append(_w_p(s[3:], style="Heading1", bold=True))
        elif s.startswith("### "):
            body_parts.append(_w_p(s[4:], style="Heading2", bold=True))
        else:
            body_parts.append(_w_p(s))

    flush_table()

    document_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        "<w:body>"
        + "".join(body_parts)
        + "<w:sectPr><w:pgSz w:w=\"12240\" w:h=\"15840\"/>"
        "<w:pgMar w:top=\"1440\" w:right=\"1440\" w:bottom=\"1440\" w:left=\"1440\"/></w:sectPr>"
        "</w:body></w:document>"
    )

    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        '<Override PartName="/word/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>'
        "</Types>"
    )

    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="word/document.xml"/>'
        "</Relationships>"
    )

    word_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rIdStyles" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        "</Relationships>"
    )

    styles_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        '<w:style w:type="paragraph" w:styleId="Title">'
        '<w:name w:val="Title"/><w:rPr><w:b/><w:sz w:val="32"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Heading1">'
        '<w:name w:val="heading 1"/><w:rPr><w:b/><w:sz w:val="26"/></w:rPr></w:style>'
        '<w:style w:type="paragraph" w:styleId="Heading2">'
        '<w:name w:val="heading 2"/><w:rPr><w:b/><w:sz w:val="22"/></w:rPr></w:style>'
        "</w:styles>"
    )

    with zipfile.ZipFile(docx_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types_xml)
        zf.writestr("_rels/.rels", rels_xml)
        zf.writestr("word/document.xml", document_xml)
        zf.writestr("word/_rels/document.xml.rels", word_rels_xml)
        zf.writestr("word/styles.xml", styles_xml)


def main() -> None:
    md_path = MS_DIR / "MANUSCRIPT.md"
    docx_out = MS_DIR / "Spaceflight_Biophysics_Multiomics_Manuscript.docx"
    build_docx(md_path, docx_out)
    print(f"Created {docx_out} ({docx_out.stat().st_size:,} bytes)")

    PUB_DIR.mkdir(parents=True, exist_ok=True)
    EXP_DIR.mkdir(parents=True, exist_ok=True)

    # Copy DOCX to public/ (both unified name and legacy name for backward compatibility)
    shutil.copy2(docx_out, PUB_DIR / "Spaceflight_Biophysics_Multiomics_Manuscript.docx")
    shutil.copy2(docx_out, PUB_DIR / "LunarLeaf-CFD-manuscript.docx")
    shutil.copy2(docx_out, EXP_DIR / "Spaceflight_Biophysics_Multiomics_Manuscript.docx")

    # Copy compiled PDFs to public/ and public/explorer/
    main_pdf = MS_DIR / "main.pdf"
    supp_pdf = MS_DIR / "supplementary.pdf"
    if main_pdf.exists():
        shutil.copy2(main_pdf, PUB_DIR / "Spaceflight_Biophysics_Multiomics_Manuscript.pdf")
        shutil.copy2(main_pdf, EXP_DIR / "Spaceflight_Biophysics_Multiomics_Manuscript.pdf")
        shutil.copy2(main_pdf, EXP_DIR / "Photorespiration_Spaceflight_Manuscript.pdf")
    if supp_pdf.exists():
        shutil.copy2(supp_pdf, PUB_DIR / "Spaceflight_Biophysics_Multiomics_SI.pdf")
        shutil.copy2(supp_pdf, EXP_DIR / "Spaceflight_Biophysics_Multiomics_SI.pdf")
        shutil.copy2(supp_pdf, EXP_DIR / "Photorespiration_Spaceflight_SI.pdf")
    print("Synced compiled PDF and DOCX manuscripts to public/ and public/explorer/")


if __name__ == "__main__":
    main()
