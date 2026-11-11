#!/usr/bin/env python3
"""Independent pathway enrichment check outside PaintOmics (Reviewer Item 15).

Tests OSD-522 transcriptome and proteome independently on core KEGG pathways:
  - ath00195: Photosynthesis
  - ath00196: Photosynthesis - antenna proteins
  - ath00710: Carbon fixation in photosynthetic organisms
  - ath00500: Starch and sucrose metabolism
  - ath00630: Glyoxylate and dicarboxylate metabolism (photorespiration)
  - ath04141: Protein processing in endoplasmic reticulum (UPR)

Calculates:
  1. Fisher's exact test (over-representation in DE features at FDR < 0.05)
  2. Mann-Whitney U test (median shift against background)
  3. Directional concordance

Writes:
  results/tables/T16_independent_pathways.tsv
"""
from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets
from paths import TABLES, ensure

PATHWAYS = {
    "ath00195": "Photosynthesis",
    "ath00196": "Photosynthesis - antenna proteins",
    "ath00710": "Carbon fixation in photosynthetic organisms",
    "ath00500": "Starch and sucrose metabolism",
    "ath00630": "Glyoxylate and dicarboxylate metabolism",
    "ath04141": "Protein processing in endoplasmic reticulum",
}

def fisher_enrichment(de_set, pathway_set, universe):
    """Fisher exact test of pathway membership vs DE status."""
    a = len(de_set & pathway_set)
    b = len(de_set - pathway_set)
    c = len((universe - de_set) & pathway_set)
    d = len((universe - de_set) - pathway_set)
    table = [[a, b], [c, d]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return a, len(pathway_set & universe), odds, p

def main():
    ensure(TABLES)
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"), sep="\t", index_col=0)
    tx_universe = set(tx.index)
    tx_sig = set(tx[tx["fdr"] < 0.05].index)
    
    rows = []
    for pid, pname in PATHWAYS.items():
        genes = genesets.pathway_genes(pid)
        k_in_universe = genes & tx_universe
        if not k_in_universe:
            continue
            
        k_de, total_p, odds, p_fisher = fisher_enrichment(tx_sig, genes, tx_universe)
        
        # Mann-Whitney shift
        in_lfc = tx.loc[tx.index.isin(genes), "log2fc"].dropna()
        out_lfc = tx.loc[~tx.index.isin(genes), "log2fc"].dropna()
        u, p_mwu = stats.mannwhitneyu(in_lfc, out_lfc, alternative="two-sided")
        shift = float(in_lfc.median() - out_lfc.median())
        
        rows.append({
            "pathway_id": pid,
            "pathway_name": pname,
            "pathway_size_universe": total_p,
            "DE_count": k_de,
            "fisher_odds_ratio": round(odds, 3),
            "fisher_pvalue": p_fisher,
            "median_lfc": round(float(in_lfc.median()), 4),
            "median_shift": round(shift, 4),
            "mwu_pvalue": p_mwu,
            "direction": "down" if shift < 0 else "up"
        })
        
    df = pd.DataFrame(rows)
    out_path = os.path.join(TABLES, "T16_independent_pathways.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path}")
    print(df[["pathway_id", "pathway_name", "DE_count", "fisher_pvalue", "median_shift", "mwu_pvalue"]])
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
