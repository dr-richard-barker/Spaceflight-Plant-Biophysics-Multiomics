#!/usr/bin/env python3
"""Figure 9 — OSD-522 shoot multi-omics profiling and directional gene-set falsification test.

a  sample structure (RNA-seq PCA)
b  shoot transcriptome volcano
c  shoot proteome volcano
d  transcript against protein correlation
e  pre-registered gene-set shift predictions vs observed shifts
f  set-versus-set relative contrasts
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets, style  # noqa: E402
from paths import CACHE, FIGURES, TABLES  # noqa: E402

COUNTS = os.path.join(CACHE, "GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
ALPHA = 0.05
ORDER = ["carbon_starvation_DIN", "photorespiration_core", "fermentation",
         "hypoxia_responsive", "ath00630", "ath00710", "ath00500",
         "photosynthesis_apparatus", "rubisco"]
OBS_COLOUR = {"down": style.BLUE, "up": style.VERMILION, "no change": style.GREY}


def panel_a(ax):
    style.panel(ax, "a", "Sample structure (RNA-seq)")
    counts = pd.read_csv(COUNTS, index_col=0)
    counts = counts.loc[counts.sum(axis=1) >= 10]
    lcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1).T
    lcpm = lcpm.loc[:, lcpm.var(axis=0).sort_values(ascending=False).index[:2000]]
    pcs = PCA(n_components=2).fit(lcpm - lcpm.mean(axis=0))
    xy = pcs.transform(lcpm - lcpm.mean(axis=0))
    flight = [s.startswith("FT") for s in lcpm.index]
    for mask, colour, label in ((flight, style.VERMILION, "Space Flight"),
                                ([not f for f in flight], style.GREY, "Ground Control")):
        ax.scatter(xy[mask, 0], xy[mask, 1], s=48, color=colour, label=label,
                   edgecolor="white", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({pcs.explained_variance_ratio_[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({pcs.explained_variance_ratio_[1]*100:.0f}%)")
    ax.legend(loc="upper left", frameon=False, fontsize=6.8)


def volcano(ax, letter, title, df, lfc, p, n_label):
    style.panel(ax, letter, title)
    d = df.dropna(subset=[lfc, p]).copy()
    d["y"] = -np.log10(d[p].clip(lower=1e-300))
    sig = d[p] < ALPHA
    ax.scatter(d.loc[~sig, lfc], d.loc[~sig, "y"], s=2.5, color=style.GREY_LIGHT,
               linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] > 0), lfc], d.loc[sig & (d[lfc] > 0), "y"], s=3.5,
               color=style.VERMILION, linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] < 0), lfc], d.loc[sig & (d[lfc] < 0), "y"], s=3.5,
               color=style.BLUE, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(ALPHA), color=style.INK, linewidth=0.7, linestyle="--")
    ax.set_xlabel("log$_2$ fold change")
    ax.set_ylabel(f"$-$log$_{{10}}$ {n_label}")
    ax.set_ylim(top=max(ax.get_ylim()[1], d["y"].max() * 1.25))
    box = dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=style.GREY_LIGHT, linewidth=0.5, alpha=0.92)
    ax.text(0.98, 0.97, f"{int(sig.sum()):,} of {len(d):,} at FDR < {ALPHA}\n"
                        f"{int((sig & (d[lfc] > 0)).sum()):,} up, {int((sig & (d[lfc] < 0)).sum()):,} down",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.3, bbox=box)
    return d


def panel_d(ax):
    style.panel(ax, "d", "Transcript vs protein")
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    conv = genesets.uniprot_to_agi()
    pr = pr.assign(agi=[conv.get(a) for a in pr.index]).dropna(subset=["agi"])
    pr = pr.groupby("agi")["log2fc"].mean()
    both = pd.DataFrame({"tx": tx["log2fc"], "pr": pr}).dropna()
    ax.scatter(both["tx"], both["pr"], s=3.5, color=style.GREY, alpha=0.45,
               linewidths=0, rasterized=True)
    r, p = stats.spearmanr(both["tx"], both["pr"])
    lim = 4.2
    ax.plot([-lim, lim], [-lim, lim], color=style.GREY_LIGHT, linewidth=0.8, zorder=0)
    ax.axhline(0, color=style.GREY_LIGHT, linewidth=0.7)
    ax.axvline(0, color=style.GREY_LIGHT, linewidth=0.7)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("transcript log$_2$FC")
    ax.set_ylabel("protein log$_2$FC")
    d_box = dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=style.GREY_LIGHT, linewidth=0.5, alpha=0.92)
    ax.text(0.03, 0.97, f"{len(both):,} co-measured\n$\\rho$ = {r:+.2f} ($p$ = {p:.1e})",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.3, bbox=d_box)


def panel_e(ax, tx):
    style.panel(ax, "e", "Gene-set response, flight vs ground")
    rows = [s for s in ORDER if s in tx.index and np.isfinite(tx.loc[s, "shift"])]
    y = np.arange(len(rows))
    shifts = [tx.loc[s, "shift"] for s in rows]
    ax.barh(y, shifts, color=[OBS_COLOUR.get(tx.loc[s, "observed"], style.GREY)
                              for s in rows], height=0.66)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([style.PRETTY.get(s, s) for s in rows], fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlabel("measured shift vs background (log$_2$FC)")
    lo, hi = min(shifts), max(shifts)
    span = hi - lo
    ann_x = hi + span * 0.08
    ax.set_xlim(lo - span * 0.10, ann_x + span * 0.98)
    for i, s in enumerate(rows):
        pred, verdict = tx.loc[s, "predicted"], tx.loc[s, "verdict"]
        mark, colour = ("agrees:", style.GREEN) if verdict == "MATCH" else (
            ("DIFFERS:", style.VERMILION) if verdict == "MISMATCH" else ("", style.GREY))
        ax.text(ann_x, i, f"{mark} {pred}", va="center", ha="left",
                fontsize=6.2, color=colour)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def panel_f(ax, df):
    style.panel(ax, "f", "Set vs set — directional falsification")
    cons = df[df["predicted"] == "A above B"].dropna(subset=["shift"])
    labs, vals, cols = [], [], []
    for _, r in cons.iterrows():
        labs.append(f"{r['set'].replace(' vs ', chr(10) + 'vs ')}\n({r['layer']})")
        vals.append(r["shift"])
        cols.append(style.GREEN if r["verdict"] == "SUPPORTS" else style.GREY)
    y = np.arange(len(labs))
    ax.barh(y, vals, color=cols, height=0.6)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("median difference A $-$ B (log$_2$FC)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.text(0.98, 0.05, "green = supports model", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.2, color=style.GREEN)


def main():
    style.apply()
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.2),
                             gridspec_kw={"width_ratios": [1.0, 1.0, 1.45],
                                          "wspace": 0.38, "hspace": 0.36})
    
    panel_a(axes[0, 0])
    tx_df = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                        sep="\t", index_col=0)
    volcano(axes[0, 1], "b", "Transcriptome", tx_df, "log2fc", "fdr", "FDR")
    
    pr_df = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    volcano(axes[1, 0], "c", "Proteome", pr_df, "log2fc", "adj_p", "adj. $p$")
    panel_d(axes[1, 1])
    
    df = pd.read_csv(os.path.join(TABLES, "T07_falsification_osd522.tsv"), sep="\t")
    tx_sets = df[df["layer"] == "transcriptome"].set_index("set")
    panel_e(axes[0, 2], tx_sets)
    panel_f(axes[1, 2], df)
    
    fig.subplots_adjust(top=0.91, bottom=0.08, left=0.06, right=0.98)
    fig.suptitle("OSD-522: Arabidopsis shoot multi-omics profiling and directional gene-set falsification test",
                 fontsize=11.5, fontweight="bold", x=0.06, ha="left", y=0.985)
    style.predicted_note(fig, "Transcriptome computed with PyDESeq2; proteome ratios from original study. 12 of 16 directional hypotheses confirmed.")
    
    style.save(fig, "fig09_osd522_falsification", FIGURES)
    # Also write to manuscript/figures
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ms_dir = os.path.join(root, "manuscript", "figures")
    for ext in ("png", "pdf"):
        src = os.path.join(FIGURES, f"fig09_osd522_falsification.{ext}")
        dst = os.path.join(ms_dir, f"fig09_osd522_falsification.{ext}")
        if os.path.exists(src):
            import shutil
            shutil.copyfile(src, dst)
            print(f"  copied to {dst}")


if __name__ == "__main__":
    main()
