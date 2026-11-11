#!/usr/bin/env python3
"""Integration control without the predicted metabolite layer (Reviewer Item 12).

Tests whether headline pathway rankings hold using RNA-seq + Proteomics alone.
Computes Fisher combined p-value across transcriptome and proteome p-values
for each pathway present in PaintOmics Job 1.

Writes:
  results/tables/T17_integration_control_nometa.tsv
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import REPO, TABLES, ensure

RAW = os.path.join(REPO, "results", "paintomics_raw", "job1_oK2h66740Q_full.json.gz")

def main():
    ensure(TABLES)
    with gzip.open(RAW, "rt") as fh:
        job = json.load(fh)
        
    rows = []
    for p in job["pathwaysInfo"]:
        name = p["name"]
        db = "KEGG" if p.get("source") == "KEGG" else "MapMan"
        g = p.get("globalOmicPvalues") or {}
        p_gene = g.get("Gene expression")
        p_prot = g.get("Proteomics")
        
        # Combine available p-values (RNA only or RNA + protein)
        valid_ps = [pv for pv in (p_gene, p_prot) if pv is not None and np.isfinite(pv) and pv > 0]
        if valid_ps:
            stat, p_comb = stats.combine_pvalues(valid_ps, method="fisher")
        else:
            p_comb = np.nan
            
        rows.append({
            "db": db,
            "pathway": name,
            "p_gene": p_gene,
            "p_protein": p_prot,
            "p_combined_RNA_prot_only": p_comb,
            "p_combined_with_metabolites": p["combinedSignificancePvalues"]["Fisher"][0]
        })
        
    df = pd.DataFrame(rows).sort_values("p_combined_RNA_prot_only").reset_index(drop=True)
    df["rank_RNA_prot_only"] = df.index + 1
    
    # Calculate FDR
    n = len(df)
    p_s = df["p_combined_RNA_prot_only"].fillna(1.0)
    order = p_s.sort_values().index
    q = (p_s.loc[order] * n / range(1, n + 1))[::-1].cummin()[::-1].clip(upper=1.0)
    df["fdr_bh_nometa"] = q.reindex(df.index)
    
    out_path = os.path.join(TABLES, "T17_integration_control_nometa.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path}")
    
    # Print key checks
    for target in ["Protein processing", "Photosynthesis", "Glyoxylate"]:
        sub = df[df["pathway"].str.contains(target, case=False)]
        for _, r in sub.iterrows():
            print(f"  {r['pathway'][:35]:<36} RNA+prot rank: {r['rank_RNA_prot_only']:>3}/{n}  "
                  f"p={r['p_combined_RNA_prot_only']:.4e}  (with meta: {r['p_combined_with_metabolites']:.4e})")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
