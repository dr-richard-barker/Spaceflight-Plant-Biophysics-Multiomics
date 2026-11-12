#!/usr/bin/env python3
"""Figure 10 — Cross-study evaluation across six OSDR flight contrasts and experimental design to unconfound hardware from microgravity.

a  five spaceflight studies, six contrasts laid out by enclosure and illumination
b  starvation-over-photosynthesis separation, lit versus dark
c  rejection of the simple one-dimensional enclosure gradient
d  multi-study gene-set response heatmap across all six flight contrasts
e  physiological confounding in standard two-arm spaceflight designs
f  proposed unconfounded three-arm and 2x2 factorial experimental designs
g  NASA OSDR archive gap: missing plant metabolomics and in-canister environmental logging
"""
from __future__ import annotations

import os
import sys
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import FIGURES, TABLES  # noqa: E402

KEY = "starvation vs photosynthesis"
SET_ORDER = ["carbon_starvation_DIN", "photorespiration_core", "ath00630",
             "fermentation", "hypoxia_responsive", "ath00500", "ath00710",
             "photosynthesis_apparatus", "rubisco", "ath00010"]


def short(acc: str) -> str:
    return acc.replace("OSD-", "").replace("-light", " lit").replace("-dark", " dark")


def box(ax, x, y, w, h, label, fc, ec, fontsize=7.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.045",
                                facecolor=fc, edgecolor=ec, linewidth=0.9))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fontsize)


def model_numbers() -> tuple[float, dict]:
    t06 = pd.read_csv(os.path.join(TABLES, "T06_predicted_compounds.tsv"), sep="\t", comment="#")
    ratio = t06["log2FC_BRIC_vs_VENTED"].abs().max() / t06["log2FC_FLT_vs_GC"].abs().max()
    t13 = pd.read_csv(os.path.join(TABLES, "T13_light_sensitivity.tsv"), sep="\t", comment="#")
    d = t13[t13["quantity"] == "deficit"]
    deficit = {(int(r.Q_umol_m2_s), r.scale): (r.deficit_min_percent, r.deficit_max_percent)
               for r in d.itertuples()}
    return ratio, deficit


def panel_a(ax, t12) -> None:
    style.panel(ax, "a", "Five studies, six contrasts")
    ax.set_xlim(-0.4, 2.6); ax.set_ylim(-0.6, 1.7)
    ax.set_xticks(range(3)); ax.set_xticklabels(["sealed\nBRIC", "tape\nCARA", "vented\nVEGGIE"], fontsize=7.5)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["dark", "lit"], fontsize=7.5)
    ax.grid(False)
    xpos = {"sealed": 0, "tape": 1, "vented": 2}
    cells: dict[tuple[int, int], list] = {}
    for _, r in t12.iterrows():
        cells.setdefault((xpos[r["enclosure_class"]], 1 if r["light"] == "light" else 0),
                         []).append(r)
    for (x, y), members in cells.items():
        n = len(members)
        for i, r in enumerate(members):
            dx = 0.0 if n == 1 else (i - (n - 1) / 2) * 0.40
            ax.scatter(x + dx, y, s=340, color=style.LIGHT[r["light"]], alpha=0.9,
                       edgecolor=style.INK, linewidth=0.9)
            ax.text(x + dx, y, short(r["accession"]).replace(" lit", "").replace(" dark", ""),
                    ha="center", va="center", fontsize=7.0, color="white", fontweight="bold")
    ax.text(2, 0, "no dark vented study exists\nin OSDR", ha="center", va="center",
            fontsize=6.5, color=style.GREY, style="italic")
    ax.text(0.5, 1.48, "carbon retained by the enclosure", ha="center", fontsize=6.8,
            color=style.GREY)
    ax.annotate("", xy=(2.4, 1.38), xytext=(0.0, 1.38),
                arrowprops=dict(arrowstyle="->", color=style.GREY, lw=0.9))


def panel_b(ax, t12) -> None:
    style.panel(ax, "b", "Lit and dark separate completely")
    t = t12.sort_values(["light", "separation"])
    colours = [style.LIGHT[l] for l in t["light"]]
    ypos = np.arange(len(t))
    ax.barh(ypos, t["separation"], color=colours, height=0.62)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{short(a)}" for a in t["accession"]], fontsize=7)
    ax.set_xlabel("starvation $-$ photosynthesis\n(median log$_2$FC difference)", fontsize=7.2)
    lo, hi = ax.get_xlim()
    ax.set_xlim(lo, hi + 3.4)

    lit = t[t["light"] == "light"]["separation"]
    dark = t[t["light"] == "dark"]["separation"]
    _u, p = stats.mannwhitneyu(lit, dark, alternative="greater")
    ax.text(0.99, 0.03, f"every lit $>$ every dark\none-sided $p$ = {p:.2f} (descriptive)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4,
            color=style.INK, linespacing=1.35)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def panel_c(ax, t12) -> None:
    style.panel(ax, "c", "The enclosure gradient does not hold")
    lit = t12[t12["light"] == "light"].dropna(subset=["separation"])
    for _, r in lit.iterrows():
        ax.scatter(r["carbon_pct_earth"], r["separation"], s=95,
                   color=style.ENCLOSURE[r["enclosure_class"]], zorder=3,
                   edgecolor="white", linewidth=1.0)
        ax.annotate(short(r["accession"]), (r["carbon_pct_earth"], r["separation"]),
                    textcoords="offset points", xytext=(0, 11), ha="center", fontsize=6.8, fontweight="bold")
    rho, prho = stats.spearmanr(lit["carbon_pct_earth"], lit["separation"])
    ax.set_xlabel("carbon retained over 12 h (% of Earth control)\n(0-D model parameterized by CFD $g_{bl}$)", fontsize=7.2)
    ax.set_ylabel("starvation $-$ photosynthesis", fontsize=7.2)
    ax.set_xlim(-14, 128)
    ax.set_ylim(-0.35, 3.7)
    xs = np.linspace(0, 100, 10)
    ax.plot(xs, 2.6 - 0.022 * xs, color=style.GREY, linestyle="--", linewidth=1.0,
            label="what the model predicted")
    ax.legend(loc="upper right", frameon=False, fontsize=6.5)
    ax.text(0.04, 0.04, f"Spearman $\\rho$ = {rho:+.2f}, $p$ = {prho:.2f}\n"
                        "three studies cannot rank three classes",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4,
            color=style.VERMILION)


def panel_d(ax, t11) -> None:
    style.panel(ax, "d", "Every gene set, every contrast (multi-study response matrix)")
    accs = ["OSD-427", "OSD-522", "OSD-678-light", "OSD-321", "OSD-38", "OSD-678-dark"]
    sets = [s for s in SET_ORDER if s in t11["set"].values]
    grid = np.full((len(sets), len(accs)), np.nan)
    for i, s in enumerate(sets):
        for j, a in enumerate(accs):
            m = t11[(t11["set"] == s) & (t11["accession"] == a)]
            if len(m):
                grid[i, j] = m["shift"].iloc[0]
    vlim = max(abs(np.nanmin(grid)), abs(np.nanmax(grid)))
    im = ax.imshow(grid, cmap="coolwarm", vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(range(len(accs)))
    ax.set_xticklabels([short(a) for a in accs], fontsize=7.5)
    ax.set_yticks(range(len(sets)))
    ax.set_yticklabels([style.PRETTY.get(s, s) for s in sets], fontsize=7.2)
    ax.axvline(2.5, color=style.INK, linewidth=1.2, linestyle="--")
    ax.text(1.0, -0.65, "lit", ha="center", fontsize=7.8, fontweight="bold", color=style.ORANGE)
    ax.text(4.0, -0.65, "dark", ha="center", fontsize=7.8, fontweight="bold", color=style.BLUE)
    cbar = plt.colorbar(im, ax=ax, fraction=0.015, pad=0.02)
    cbar.set_label("shift vs all other genes (log$_2$FC)", fontsize=7.0)
    cbar.ax.tick_params(labelsize=6.8)
    for i in range(len(sets)):
        for j in range(len(accs)):
            val = grid[i, j]
            if np.isfinite(val):
                txt_col = "white" if abs(val) > vlim * 0.55 else style.INK
                ax.text(j, i, f"{val:+.2f}", ha="center", va="center",
                        fontsize=6.3, color=txt_col)


def panel_e(ax):
    style.panel(ax, "e", "What the standard contrast cannot see")
    ratio, _ = model_numbers()
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    box(ax, 0.3, 4.25, 4.3, 1.4, "FLIGHT\nsealed, lit, $\\mu$g", "#FBD5BC", "#8C5A3C", fontsize=7.0)
    box(ax, 0.3, 2.05, 4.3, 1.4, "GROUND CONTROL\nsealed, lit, 1 $g$", "#E4E7EB", "#5A6169", fontsize=7.0)
    ax.add_patch(FancyArrowPatch((4.75, 4.95), (5.75, 3.95), arrowstyle="-|>",
                                 mutation_scale=12, color=style.INK, linewidth=1.3))
    ax.add_patch(FancyArrowPatch((4.75, 2.75), (5.75, 3.55), arrowstyle="-|>",
                                 mutation_scale=12, color=style.INK, linewidth=1.3))
    box(ax, 5.85, 3.0, 4.0, 1.5, "the published\nflight-vs-ground\ncontrast", "white",
        style.INK, fontsize=7.0)
    ax.text(2.45, 6.05, "CO$_2$ drawn down in BOTH arms", ha="center", fontsize=7.4,
            color=style.VERMILION, fontweight="bold")
    ax.text(7.85, 2.50, f"hardware effect\nsubtracts out $\\rightarrow$\n{ratio:.0f}× the gravity\neffect, invisible",
            ha="center", va="top", fontsize=6.8, color=style.VERMILION)


def panel_f(ax):
    style.panel(ax, "f", "The experimental design that separates them")
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    box(ax, 0.3, 4.45, 4.4, 1.15, "ARM 1  flight, sealed", "#FBD5BC", "#8C5A3C", fontsize=6.8)
    box(ax, 0.3, 2.65, 4.4, 1.15, "ARM 2  ground, sealed", "#E4E7EB", "#5A6169", fontsize=6.8)
    box(ax, 0.3, 0.85, 4.4, 1.15, "ARM 3  ground, vented", "#D2E4F7", "#2D5887", fontsize=6.8)
    ax.text(5.05, 3.65, "1 − 2 = gravity", va="center", fontsize=7.2, fontweight="bold")
    ax.text(5.05, 1.85, "2 − 3 = hardware", va="center", fontsize=7.2, fontweight="bold",
            color=style.VERMILION)
    ax.text(5.05, 5.25, "Three arms estimate main effects (gravity & hardware).\n"
                        "Adding ARM 4 (flight, vented) estimates\n"
                        "the gravity × hardware interaction.",
            va="center", fontsize=6.5, color=style.GREY)


def panel_g(ax):
    style.panel(ax, "g", "The NASA OSDR archive gap")
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    entries = [
        (5.3, "0 of 12", "OSDR metabolomics studies that are of a plant"),
        (3.3, "12 of 567", "OSDR studies that are metabolite profiling,\nall mouse, human or microbe"),
        (1.3, "0 known", "in-canister CO$_2$ logs, any BRIC-class flight"),
    ]
    for y, big, desc in entries:
        ax.text(0.3, y, big, fontsize=11.5, fontweight="bold", color=style.VERMILION, va="center")
        ax.text(3.3, y, desc, fontsize=6.8, va="center")
    ax.text(0.3, -0.05, "A ~20-compound targeted panel on archived material\n"
                        "would be the single most informative experiment.",
            fontsize=6.5, color=style.INK, style="italic")


def main() -> None:
    style.apply()
    t11 = pd.read_csv(os.path.join(TABLES, "T11_hardware_ladder.tsv"), sep="\t")
    t12 = pd.read_csv(os.path.join(TABLES, "T12_ladder_contrasts.tsv"), sep="\t")
    t12 = t12[t12["contrast"] == KEY].copy()

    # Near full-page proportion: 8.5 in x 10.8 in
    fig = plt.figure(figsize=(8.5, 10.8))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.92, 1.12, 0.95],
                          hspace=0.42, wspace=0.38,
                          top=0.94, bottom=0.05, left=0.07, right=0.96)

    # Row 0: Three contrast plots
    panel_a(fig.add_subplot(gs[0, 0]), t12)
    panel_b(fig.add_subplot(gs[0, 1]), t12)
    panel_c(fig.add_subplot(gs[0, 2]), t12)

    # Row 1: The multi-study response heatmap
    panel_d(fig.add_subplot(gs[1, :]), t11)

    # Row 2: Synthesis schematics and archive gap
    panel_e(fig.add_subplot(gs[2, 0]))
    panel_f(fig.add_subplot(gs[2, 1]))
    panel_g(fig.add_subplot(gs[2, 2]))

    fig.suptitle("Cross-study evaluation across six OSDR flight contrasts and experimental design to unconfound hardware from microgravity",
                 fontsize=11.5, fontweight="bold", x=0.07, ha="left", y=0.985)
    style.predicted_note(fig, "Six matched flight-vs-ground contrasts from NASA OSDR; gene sets from KEGG; model from LunarLeaf-CFD & FvCB.")

    style.save(fig, "fig10_hardware_ladder_implications", FIGURES)
    
    # Also write to manuscript/figures
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ms_dir = os.path.join(root, "manuscript", "figures")
    for ext in ("png", "pdf"):
        src = os.path.join(FIGURES, f"fig10_hardware_ladder_implications.{ext}")
        dst = os.path.join(ms_dir, f"fig10_hardware_ladder_implications.{ext}")
        if os.path.exists(src):
            shutil.copyfile(src, dst)
            print(f"  copied to {dst}")


if __name__ == "__main__":
    main()
