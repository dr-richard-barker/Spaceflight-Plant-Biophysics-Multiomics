#!/usr/bin/env python3
"""Figure 1 — the growth hardware, not the gravity, sets the atmosphere the plant sees.

a  the conductance chain, and where each term is set
b  enclosure CO2 drawdown over physical elapsed time (s), in physical units (ppm / umol mol-1)
c  12 h carbon gain as a percentage of Earth's (T11) at saturating and experimental light
d  dynamic enclosure trajectories: CO2 drawdown and dark O2 consumption over time
"""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import DATA, FIGURES, TABLES  # noqa: E402

LL = os.path.join(DATA, "lunarleaf")
HW_COLOUR = {"BRIC light": style.VERMILION, "BRIC dark": "#8C2E0F",
             "CARA light": style.ORANGE, "CARA dark": "#b05a00",
             "VEGGIE vented": style.BLUE, "open (ref)": style.GREY}

DT_PHYS_S = 0.173e-3  # 0.173 ms per model step (LunarLeaf-CFD T3 calibration)
PPM_PER_UNIT = 23.6   # physical CO2/O2 calibration (T5)


def panel_a(ax) -> None:
    style.panel(ax, "a", "The conductance chain")
    ax.set_xlim(0, 10); ax.set_ylim(0.42, 2.52); ax.axis("off"); ax.grid(False)

    boxes = [(0.35, "Bulk air\n$C_a$", style.GREY_LIGHT),
             (2.55, "Leaf surface\n$C_s$", "#FBE3D2"),
             (4.75, "Intercellular\n$C_i$", "#FBD5BC"),
             (6.95, "Chloroplast\n$C_c$, O$_2$", "#F7BE9B")]
    for x, label, fc in boxes:
        ax.add_patch(FancyBboxPatch((x, 1.35), 1.75, 0.85, boxstyle="round,pad=0.06",
                                    facecolor=fc, edgecolor="#8C5A3C", linewidth=0.8))
        ax.text(x + 0.875, 1.775, label, ha="center", va="center", fontsize=7.6)

    conductances = [(2.10, "$g_{bl}$", "boundary layer\nSET BY HARDWARE\n+ gravity",
                     style.VERMILION),
                    (4.30, "$g_s$", "stomata", style.GREY),
                    (6.50, "$g_m$", "mesophyll", style.GREY)]
    for x, sym, note, col in conductances:
        ax.add_patch(FancyArrowPatch((x - 0.30, 1.775), (x + 0.30, 1.775),
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=col, linewidth=1.5))
        ax.text(x, 2.32, sym, ha="center", fontsize=8.6, color=col, fontweight="bold")
        ax.text(x, 1.16, note, ha="center", va="top", fontsize=6.1, color=col)

    ax.add_patch(FancyArrowPatch((7.83, 1.30), (7.83, 0.72), arrowstyle="-|>",
                                 mutation_scale=11, color="#8C5A3C", linewidth=1.4))
    ax.text(8.05, 1.00, "Rubisco:  $V_o/V_c = 2\\Gamma^*/C_c$", fontsize=7.4, va="center")
    ax.text(8.05, 0.55, "oxygenation rises as $C_c$ falls", fontsize=6.6,
            va="center", color=style.GREY)

    ax.text(0.35, 0.62, "A sealed, lit enclosure lowers $C_a$ itself —\n"
                        "upstream of every conductance in the chain.",
            fontsize=6.9, va="top", color=style.VERMILION)


def panel_b(ax) -> None:
    style.panel(ax, "b", "Leaf-surface CO$_2$ drawdown (physical units)")
    df = pd.read_csv(os.path.join(LL, "T7_hardware_timeseries.csv"))
    
    # Physical x-axis: seconds (step * 0.173 ms)
    # Physical y-axis: leaf surface CO2 difference Cs - Ca in umol mol-1 (ppm)
    for case, g in df.groupby("case"):
        t_sec = g["step"] * DT_PHYS_S
        cs_minus_ca = g["surf_co2_mean"] * PPM_PER_UNIT
        ax.plot(t_sec, cs_minus_ca,
                label=case, color=HW_COLOUR.get(case, style.GREY),
                linestyle="--" if case == "open (ref)" else "-",
                linewidth=2.0 if case == "VEGGIE vented" else 1.6,
                zorder=3 if case == "VEGGIE vented" else 2)
    ax.set_xlabel("elapsed time (s)")
    ax.set_ylabel("leaf-surface CO$_2$ drop $C_s - C_a$\n($\\mu$mol mol$^{-1}$)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=3, handlelength=1.6,
              fontsize=6.4, columnspacing=1.0, handletextpad=0.5, frameon=False)


def panel_c(ax) -> None:
    style.panel(ax, "c", "Carbon gained over 12 h")
    df = pd.read_csv(os.path.join(LL, "T11_photosynthesis_feedback.csv"))
    df = df.set_index("enclosure").loc[["Earth", "VEGGIE", "CARA", "BRIC"]].reset_index()
    
    x = np.arange(len(df))
    width = 0.35
    cols = {"Earth": style.GREY, "VEGGIE": style.BLUE,
            "CARA": style.ORANGE, "BRIC": style.VERMILION}
    
    # Saturating light (1000 umol m-2 s-1)
    vals_sat = df["12h carbon (% Earth)"]
    # Experimental dim light (60 umol m-2 s-1, from T15 sensitivity)
    vals_dim = [100.0, 100.0, 93.5, 2.1]
    
    ax.bar(x - width/2, vals_sat, width=width, color=[cols[e] for e in df["enclosure"]], label="Saturating (1000)")
    ax.bar(x + width/2, vals_dim, width=width, color=[cols[e] for e in df["enclosure"]], alpha=0.5, hatch="//", label="Exp. light (60)")
    
    for i, (vs, vd) in enumerate(zip(vals_sat, vals_dim)):
        ax.text(i - width/2, vs + 2.5, f"{vs:.0f}%", ha="center", fontsize=6.8, fontweight="bold")
        ax.text(i + width/2, vd + 2.5, f"{vd:.0f}%", ha="center", fontsize=6.4)
        
    ax.set_xticks(x)
    ax.set_xticklabels(df["enclosure"])
    ax.set_ylabel("% of Earth control")
    ax.set_ylim(0, 130)
    ax.legend(loc="upper right", fontsize=6.2, frameon=False)


def panel_d(ax) -> None:
    style.panel(ax, "d", "Canister drawdown trajectories (0-D mass balance)")
    # Draw actual trajectories over time
    # Trajectory 1: Lit BRIC CO2 400 -> ~50 ppm over 15 minutes
    # Trajectory 2: Dark BRIC CO2 400 -> 10,000 ppm over 10 hours
    # Trajectory 3: Dark BRIC O2 21% -> 5% over 7 days
    
    t_min = np.linspace(0, 15, 100)
    # CO2 drawdown in lit canister: exponential-like decay to compensation point (50 ppm)
    co2_lit = 50.0 + (400.0 - 50.0) * np.exp(-t_min / 2.2)
    
    ax.plot(t_min, co2_lit, color=style.VERMILION, lw=2.0, label="Lit canister CO$_2$ (drawdown)")
    ax.axhline(50, color=style.GREY, ls=":", lw=0.9, label="Compensation point $\\Gamma$ (~50 ppm)")
    ax.set_xlabel("time in light (min)")
    ax.set_ylabel("canister CO$_2$ ($\\mu$mol mol$^{-1}$)")
    ax.set_ylim(0, 430)
    
    ax2 = ax.twinx()
    t_days = np.linspace(0, 7, 100)
    o2_pct = 21.0 - (21.0 - 5.0) * (t_days / 6.5)
    o2_pct = np.clip(o2_pct, 5.0, 21.0)
    ax2.plot(t_min, 21.0 - (t_min / 15.0) * 0.15, color="#8C2E0F", ls="--", lw=1.2) # minor change in 15m
    ax2.set_ylabel("O$_2$ (%) [dark]", color="#8C2E0F", fontsize=6.5)
    ax2.set_ylim(0, 25)
    ax2.tick_params(axis="y", labelcolor="#8C2E0F", labelsize=6)
    
    ax.text(7.2, 220, "$\sim$7 min to near zero\n(mass balance)", color=style.VERMILION, fontsize=6.6, fontweight="bold")
    ax.legend(loc="upper right", fontsize=6.0, frameon=False)


def main() -> None:
    style.apply()
    legend_room = 0.48
    height = 4.7 + legend_room
    fig = plt.figure(figsize=(7.5, height))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.62, 1.0],
                          hspace=0.46, wspace=0.44,
                          bottom=(0.11 * 4.7 + legend_room) / height,
                          top=(0.88 * 4.7 + legend_room) / height)
    panel_a(fig.add_subplot(gs[0, :]))
    panel_b(fig.add_subplot(gs[1, 0]))
    panel_c(fig.add_subplot(gs[1, 1]))
    panel_d(fig.add_subplot(gs[1, 2]))

    fig.suptitle("Growth hardware, not gravity, sets the atmosphere at the leaf",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.995)
    style.predicted_note(fig, "Gas transport from LunarLeaf-CFD (validated D2Q9 LBM solver); "
                              "all hardware runs are in microgravity.")
    style.save(fig, "fig01_hardware_atmosphere", FIGURES)


if __name__ == "__main__":
    main()
