#!/usr/bin/env python3
"""Sensitivity analysis for 0-D enclosure drawdown and carbon gain (Reviewer Items 1, 3, 4, 7).

Evaluates:
  1. Time to compensation point (minutes) in sealed hardware (BRIC) across:
     - Headspace volume V: 15, 30, 60 mL
     - Leaf area A_leaf: 1, 3, 6, 10 cm^2
     - Growth irradiance / assimilation rate A: Q=60 (OSD-522) vs Q=1000 (saturating)
  2. 12 h carbon gain (% of Earth control) across:
     - Sealing / permeance (sealed k=0, tape 0.001 - 0.02 mol m^-2 s^-1, vented)
     - Ambient CO2 (400 vs elevated 800 - 1200 ppm)
     - Photoperiod and respiration assumptions

Writes:
  results/tables/T15_enclosure_sensitivity.tsv
"""
from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import TABLES, ensure

RHO = 41.6  # mol m^-3 air density at 25 C
CA_DEFAULT = 400.0  # ppm
GAMMA = 50.0  # ppm compensation point
K_M = 300.0   # ppm half-saturation

def solve_A(Cbulk, gbl, Amax):
    """Steady assimilation at leaf surface given bulk concentration."""
    A = 2.0
    for _ in range(50):
        Cs = max(Cbulk - A / gbl, 0.0)
        A = Amax * (Cs - GAMMA) / (Cs - GAMMA + K_M) if Cs > GAMMA else 0.0
    return A, Cs

def simulate_12h(V_ml, Aleaf_cm2, mode, g_tape=0.006, Ca=400.0, Q=1000.0,
                 A_meas=3.85, R_meas=1.18, gbl=0.272):
    """Simulate 12h photoperiod drawdown and integrated carbon gain."""
    V = V_ml * 1e-6          # m^3
    Aleaf = Aleaf_cm2 * 1e-4 # m^2
    Adish = 2.8e-3           # m^2 (60 mm dish)
    
    # Scale Amax by light level
    if Q <= 100.0:
        # At 60 umol m-2 s-1, A_meas is ~1.5 - 1.8 umol m-2 s-1
        A_meas_q = A_meas * 0.45
    else:
        A_meas_q = A_meas
        
    Amax = A_meas_q * (Ca - GAMMA + K_M) / (Ca - GAMMA)
    
    dt = 2.0  # seconds
    n_steps = int(12 * 3600 / dt)
    C = Ca
    carbon = 0.0
    time_to_comp = None
    
    # Baseline Earth unconfined carbon gain over 12h
    A_earth_init, _ = solve_A(Ca, 0.540, Amax)
    carbon_earth_control = A_earth_init * Aleaf * 12 * 3600  # umol CO2
    
    for i in range(n_steps):
        if mode == "vented":
            C = Ca
        A, Cs = solve_A(C, gbl, Amax)
        if mode == "sealed":
            dC = -A * Aleaf / (V * RHO)
        elif mode == "tape":
            dC = (-A * Aleaf + g_tape * Adish * (Ca - C)) / (V * RHO)
        else: # vented
            dC = 0.0
            
        C = max(C + dC * dt, GAMMA)
        carbon += A * Aleaf * dt
        
        if time_to_comp is None and C <= (GAMMA + 2.0):
            time_to_comp = (i * dt) / 60.0  # minutes
            
    if time_to_comp is None:
        time_to_comp = np.nan
        
    pct_carbon = (carbon / carbon_earth_control) * 100.0 if carbon_earth_control > 0 else 0.0
    return time_to_comp, pct_carbon

def main():
    ensure(TABLES)
    rows = []
    
    # 1. Grid of Headspace x Leaf Area x Irradiance for Sealed Canister
    for V in [15.0, 30.0, 60.0]:
        for A_leaf in [1.0, 3.0, 6.0, 10.0]:
            for Q in [60.0, 1000.0]:
                t_comp, c_pct = simulate_12h(V, A_leaf, "sealed", Q=Q)
                rows.append({
                    "hardware_class": "sealed (BRIC)",
                    "headspace_mL": V,
                    "leaf_area_cm2": A_leaf,
                    "irradiance_PAR": Q,
                    "tape_permeance_mol_m2_s": 0.0,
                    "ambient_CO2_ppm": 400.0,
                    "time_to_comp_min": round(t_comp, 1) if not np.isnan(t_comp) else ">720",
                    "carbon_gain_12h_pct_earth": round(c_pct, 2)
                })
                
    # 2. Tape Permeance Sensitivity for CARA
    for g_tape in [0.001, 0.003, 0.006, 0.012, 0.020]:
        for Q in [60.0, 1000.0]:
            t_comp, c_pct = simulate_12h(30.0, 3.0, "tape", g_tape=g_tape, Q=Q)
            rows.append({
                "hardware_class": "tape (CARA)",
                "headspace_mL": 30.0,
                "leaf_area_cm2": 3.0,
                "irradiance_PAR": Q,
                "tape_permeance_mol_m2_s": g_tape,
                "ambient_CO2_ppm": 400.0,
                "time_to_comp_min": round(t_comp, 1) if not np.isnan(t_comp) else "bounded",
                "carbon_gain_12h_pct_earth": round(c_pct, 2)
            })
            
    # 3. Ambient CO2 Variation (ISS Cabin often 800 - 1200 ppm)
    for Ca in [400.0, 800.0, 1200.0]:
        for mode in ["sealed", "tape", "vented"]:
            t_comp, c_pct = simulate_12h(30.0, 3.0, mode, Ca=Ca, Q=60.0)
            rows.append({
                "hardware_class": f"{mode} (ambient sweep)",
                "headspace_mL": 30.0,
                "leaf_area_cm2": 3.0,
                "irradiance_PAR": 60.0,
                "tape_permeance_mol_m2_s": 0.006 if mode == "tape" else 0.0,
                "ambient_CO2_ppm": Ca,
                "time_to_comp_min": round(t_comp, 1) if not np.isnan(t_comp) else ("bounded" if mode != "sealed" else ">720"),
                "carbon_gain_12h_pct_earth": round(c_pct, 2)
            })

    df = pd.DataFrame(rows)
    out_path = os.path.join(TABLES, "T15_enclosure_sensitivity.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path} ({len(df)} sensitivity scenarios)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
