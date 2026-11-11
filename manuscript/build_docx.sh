#!/usr/bin/env bash
# Build Spaceflight_Biophysics_Multiomics_Manuscript.docx and
# Spaceflight_Biophysics_Multiomics_SI.docx from the LaTeX sources with pandoc,
# embedding high-res PNG figures and resolving cross-references.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/figures"
cp "$ROOT_DIR/results/figures"/*.png "$tmp/figures/" 2>/dev/null || true
cp "$ROOT_DIR/results/figures/supplementary"/*.png "$tmp/figures/" 2>/dev/null || true
cp "$ROOT_DIR/public/figures"/*.png "$tmp/figures/" 2>/dev/null || true
cp "$ROOT_DIR/Microgravity_Stagnation_in_Spaceflight_Plants.png" "$tmp/figures/" 2>/dev/null || true
cp "$SCRIPT_DIR/references.bib" "$tmp/"

# Process main.tex
sed -E 's/(\\includegraphics(\[[^]]*\])?\{[^}]*)\.pdf\}/\1.png}/' "$SCRIPT_DIR/main.tex" > "$tmp/main.tex"
# Process supplementary.tex and resolve xr references from main.aux
sed -E 's/(\\includegraphics(\[[^]]*\])?\{[^}]*)\.pdf\}/\1.png}/' "$SCRIPT_DIR/supplementary.tex" > "$tmp/supplementary.tex"
perl -e 'open A, "'"$SCRIPT_DIR"'/main.aux" or die "main.aux missing: run latexmk on main.tex first\n";
         while (<A>) { $n{$1} = $2 if /\\newlabel\{([^}]+)\}\{\{([^}]*)\}/ }
         local $/; $_ = <STDIN>;
         s/\\ref\{([^}]+)\}/exists $n{$1} ? $n{$1} : "\\ref{$1}"/ge; print' \
    < "$tmp/supplementary.tex" > "$tmp/supplementary.resolved" && mv "$tmp/supplementary.resolved" "$tmp/supplementary.tex"

echo "Building main docx with pandoc..."
(cd "$tmp" && /opt/homebrew/bin/pandoc main.tex --citeproc --bibliography=references.bib \
    --resource-path=.:figures -o Spaceflight_Biophysics_Multiomics_Manuscript.docx)

echo "Building SI docx with pandoc..."
(cd "$tmp" && /opt/homebrew/bin/pandoc supplementary.tex --citeproc --bibliography=references.bib \
    --resource-path=.:figures -o Spaceflight_Biophysics_Multiomics_SI.docx)

cp "$tmp/Spaceflight_Biophysics_Multiomics_Manuscript.docx" "$SCRIPT_DIR/"
cp "$tmp/Spaceflight_Biophysics_Multiomics_SI.docx" "$SCRIPT_DIR/"

# Also sync to public/ and public/explorer/
mkdir -p "$ROOT_DIR/public/explorer"
cp "$SCRIPT_DIR/Spaceflight_Biophysics_Multiomics_Manuscript.docx" "$ROOT_DIR/public/"
cp "$SCRIPT_DIR/Spaceflight_Biophysics_Multiomics_SI.docx" "$ROOT_DIR/public/"
cp "$SCRIPT_DIR/Spaceflight_Biophysics_Multiomics_Manuscript.docx" "$ROOT_DIR/public/explorer/"
cp "$SCRIPT_DIR/Spaceflight_Biophysics_Multiomics_SI.docx" "$ROOT_DIR/public/explorer/"

echo "Successfully built embedded-figure DOCX files:"
ls -lh "$SCRIPT_DIR"/*.docx
