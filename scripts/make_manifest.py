#!/usr/bin/env python3
"""Write unified FAIR provenance manifests across results/tables/, results/cfd/tables/,
results/omics/tables/, and results/MANIFEST.tsv.

FAIR data requires every derived table to document its generating script, byte size,
description, and upstream source data. This script fails with a non-zero exit code if
any table in results/cfd/tables/, results/omics/tables/, or results/tables/ is undescribed.
"""
from __future__ import annotations
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CFD_TABLES, OMICS_TABLES, RESULTS, TABLES  # noqa: E402

CFD_DESCRIBED = {
    "T1_measured_gas_exchange.csv": (
        "validation/analyze_data.py",
        "Measured whole-chamber Arabidopsis gas-exchange fluxes (net assimilation 3.85 umol m-2 s-1, dark respiration 1.18)",
        "Chew et al. (2022) Edinburgh DataShare doi:10.7488/ds/3453 + Vernier whole-chamber logs",
    ),
    "T2_model_sweep.csv": (
        "validation/export_cfd.ts",
        "Steady-state LBM CFD output across 8 gravity x canopy-scale scenarios: u_max, Ra, mean/peak dC for H2O, CO2, O2",
        "2D D2Q9/D2Q5 lattice-Boltzmann solver (150k steps; 900k for 0g canopy)",
    ),
    "T3_calibration.csv": (
        "validation/analyze_model.py",
        "Lattice-to-physical unit calibration (dx = 0.288 mm, dt = 0.173 ms, velocity scale, empirical flux anchor)",
        "T1_measured_gas_exchange.csv + T2_model_sweep.csv",
    ),
    "T4_chamber_accumulation.csv": (
        "validation/chamber_sim.ts",
        "Closed-chamber mass conservation time series (relative error 2e-5) and near-leaf vs far-wall concentration gap",
        "2D D2Q9/D2Q5 lattice-Boltzmann sealed-chamber simulation",
    ),
    "T5_physical_prediction.csv": (
        "validation/analyze_model.py",
        "Predicted leaf-surface CO2 drawdown and O2 build-up in ppm across scale x gravity, anchored on measured flux",
        "T1_measured_gas_exchange.csv + T2_model_sweep.csv (steady state)",
    ),
    "T6_forced_airflow.csv": (
        "validation/fan_sweep.ts",
        "Forced-ventilation sweep for an isolated microgravity leaf (pressure outlet + density-matched inlet, steady state)",
        "2D D2Q9/D2Q5 lattice-Boltzmann ventilation channel sweep",
    ),
    "T7_hardware_timeseries.csv": (
        "validation/hardware_sim.ts",
        "Spaceflight hardware (BRIC sealed, CARA micropore tape, VEGGIE vented) enclosure CO2 drift and surface gap vs step",
        "2D D2Q9/D2Q5 lattice-Boltzmann hardware boundary-condition simulations",
    ),
    "T8_enclosure_timescales.csv": (
        "validation/hardware_analysis.py",
        "Analytic sealed-dish (BRIC, 30 mL headspace, 3 cm2 leaf) atmosphere timescales: CO2 exhaustion (min) and O2 hypoxia (d)",
        "0-D mass balance anchored to T1_measured_gas_exchange.csv",
    ),
    "T9_fan_by_scale.csv": (
        "validation/scales.ts",
        "Forced-airflow sweep across leaf, rosette, and canopy scales with Mach number and outlet throughput verification",
        "2D D2Q9/D2Q5 lattice-Boltzmann multi-scale ventilation sweep",
    ),
    "T10_hardware_by_scale.csv": (
        "validation/scales.ts",
        "BRIC, CARA, and VEGGIE surface gradients and dish-mean CO2 across leaf, rosette, and canopy scales",
        "2D D2Q9/D2Q5 lattice-Boltzmann hardware x scale simulations",
    ),
    "T11_photosynthesis_feedback.csv": (
        "validation/feedback_0d.py",
        "0-D closed-loop CO2-limited photosynthesis over a 12 h photoperiod: start/end assimilation and integrated carbon gain (% Earth)",
        "Coupled 0-D enclosure mass balance + rectangular-hyperbola CO2 response",
    ),
    "T12_feedback_spatial.csv": (
        "validation/feedback_test.ts",
        "Spatial closed-loop CO2-limited net assimilation (% of potential) across gravity, scale, and hardware at steady state",
        "2D D2Q9/D2Q5 lattice-Boltzmann solver with local surface CO2 feedback (Ca = 8)",
    ),
    "T13_boundary_layer.csv": (
        "validation/export_cfd.ts",
        "Steady-state boundary-layer conductance g_bl (mol m-2 s-1), diffusive film delta (mm), Sherwood number, and 5s drift",
        "2D D2Q9/D2Q5 lattice-Boltzmann solver (handoff table to scripts/fvcb.py)",
    ),
    "T14_shape_sweep.csv": (
        "validation/export_shape.ts",
        "Leaf-shape sensitivity sweep at steady state: ellipse vs realistic outline (longitudinal) and flat vs curled blade (transverse)",
        "2D D2Q9/D2Q5 lattice-Boltzmann solver + NeuraLeaf cross-section geometry (src/leaf/)",
    ),
    "T15_stationarity.csv": (
        "analysis/transient_stats.py",
        "Snapshot (30k steps) vs steady-state (150k steps) g_bl stationarity and drift verification across 6 scenarios",
        "validation/timeseries.ts time-series exports",
    ),
    "T16_resolution_gci.csv": (
        "analysis/gci.py",
        "Three-grid resolution convergence ladder (dx = 0.288, 0.192, 0.144 mm) and Celik et al. (2008) Grid Convergence Index (GCI)",
        "validation/resolution_ladder.ts multi-grid simulations",
    ),
    "T17_domain_size.csv": (
        "validation/domain_size.ts",
        "Chamber domain-size sensitivity (1x, 1.5x, 2x production box) for leaf and rosette at 1g and 0g",
        "2D D2Q9/D2Q5 lattice-Boltzmann domain-scaling simulations",
    ),
    "T18_surface_distribution.csv": (
        "validation/surface_dist.ts",
        "Spatial distribution of leaf-surface CO2 drawdown (CoV, 10th/50th/90th percentiles, fraction worse than Earth mean)",
        "2D D2Q9/D2Q5 lattice-Boltzmann steady-state surface cell distributions",
    ),
    "T19_dish_hardware.csv": (
        "validation/dish_hardware.ts",
        "Leaf 1 mm above agar inside real BRIC (60 mm round) and CARA (120 mm square) dish cross-sections at 1g and 0g",
        "2D D2Q9/D2Q5 lattice-Boltzmann real dish geometry simulations",
    ),
    "T20_dish_feedback.csv": (
        "validation/dish_feedback.ts",
        "Time series (415 s) for sealed BRIC and taped CARA dishes with CO2-limited photosynthesis and humidity-limited transpiration",
        "2D D2Q9/D2Q5 lattice-Boltzmann coupled transient dish simulations",
    ),
    "T20_dish_feedback_summary.csv": (
        "validation/dish_feedback.ts",
        "Summary metrics for sealed BRIC and taped CARA transient dish feedback simulations at 1g and 0g",
        "T20_dish_feedback.csv",
    ),
    "T20_dish_feedback_constant_transpiration.csv": (
        "validation/dish_feedback.ts",
        "Reference sealed BRIC transient series with constant (unbounded) transpiration showing late-time 1g buoyancy runaway",
        "2D D2Q9/D2Q5 lattice-Boltzmann diagnostic simulation",
    ),
    "T21_dish_buoyancy_diag.csv": (
        "validation/dish_buoyancy_diag.ts",
        "Sealed BRIC 1g buoyancy decomposition (all species vs H2O-only vs CO2+O2-only vs fixed uptake)",
        "2D D2Q9/D2Q5 lattice-Boltzmann buoyancy diagnostic series",
    ),
    "T22_dish_humidity.csv": (
        "validation/dish_humidity.ts",
        "Sealed BRIC dish sensitivity across three humidity saturation thresholds (h2oSatExcess = 0.65, 1.3, 2.6) at 1g and 0g",
        "2D D2Q9/D2Q5 lattice-Boltzmann humidity-limited transpiration sweep",
    ),
    "T23_dish_wet_agar.csv": (
        "validation/dish_wet_agar.ts",
        "Wet-agar boundary condition (agarH2O = 1.3) in BRIC and CARA at 0g (valid) and 1g (negative result: exceeds Boussinesq limit)",
        "2D D2Q9/D2Q5 lattice-Boltzmann BGK wet-agar diagnostic runs",
    ),
    "T24_dish_wet_agar_mrt.csv": (
        "validation/dish_wet_agar_mrt.ts",
        "Wet-agar diagnostic rerun with multiple-relaxation-time (MRT, Lallemand & Luo 2000) collision operator",
        "2D D2Q9/D2Q5 lattice-Boltzmann MRT wet-agar diagnostic runs",
    ),
    "inlet_flux_before.csv": (
        "validation/inlet_flux.ts (INLET_FLUX_BC=equilibrium)",
        "Forced-airflow outlet mass-flux verification under legacy rho=1 inlet (85-98% delivered flow)",
        "2D D2Q9 lattice-Boltzmann inlet/outlet diagnostic",
    ),
    "inlet_flux_after.csv": (
        "validation/inlet_flux.ts (INLET_FLUX_BC=densityMatched)",
        "Forced-airflow outlet mass-flux verification under density-matched inlet (99.6-100.9% delivered flow for Ma <= 0.3)",
        "2D D2Q9 lattice-Boltzmann inlet/outlet diagnostic",
    ),
    "inlet_variants.csv": (
        "validation/inlet_variants.ts",
        "Stability and throughput comparison of 5 inlet/collision variants (equilibrium, densityMatched, Zou-He; BGK, MRT) at high fan speeds",
        "2D D2Q9 lattice-Boltzmann high-speed stability diagnostic",
    ),
}

OMICS_DESCRIBED = {
    "T01_osd522_transcriptome.tsv": (
        "scripts/01_dge_rnaseq.py",
        "OSD-522 shoot transcriptome, flight vs ground: log2FC, p, FDR per AGI locus",
        "NASA OSDR OSD-522 RSEM counts",
    ),
    "T02_osd522_proteome.tsv": (
        "scripts/02_prep_proteomics.py",
        "OSD-522 shoot proteome, flight vs ground: deposited S/G log2 ratio and adjusted p",
        "NASA OSDR OSD-522 shoot SOL + MEM proteomics reports",
    ),
    "T03_osdr_metabolomics_census.tsv": (
        "scripts/osdr_census.py",
        "Every OSDR study with a metabolite-profiling assay, and its organism: none is a plant",
        "OSDR bulk study search, cached by scripts/probe_studies.py",
    ),
    "T04_photorespiration_vs_gravity.csv": (
        "scripts/fvcb.py",
        "FvCB operating point across the CFD gravity/scale sweep: Cc, phi, A, Rp, Rp/A",
        "results/cfd/tables/T13_boundary_layer.csv",
    ),
    "T05_operating_points.tsv": (
        "scripts/03_predict_metabolome.py",
        "The four modelled arms (flight, ground, sealed, vented): g_bl, Ca, Cc, phi, A, Vo, Vc",
        "LunarLeaf-CFD T13 + scripts/fvcb.py",
    ),
    "T06_predicted_compounds.tsv": (
        "scripts/03_predict_metabolome.py",
        "21 PREDICTED metabolite pools with KEGG id, driver flux, tier and both contrasts. MODEL OUTPUT, NEVER MEASURED.",
        "fvcb.py fluxes; KEGG ids from rest.kegg.jp",
    ),
    "T07_falsification_osd522.tsv": (
        "scripts/04_falsification_check.py",
        "Blind prediction vs measurement on OSD-522: gene-set shifts, Mann-Whitney p, verdicts",
        "T01, T02, KEGG gene sets",
    ),
    "T08_paintomics_significant.tsv": (
        "scripts/06_paintomics_tables.py (job oK2h66740Q)",
        "All 32 pathways significant at combined Fisher p < 0.05, with per-omic p-values",
        "results/omics/paintomics_raw/job1_oK2h66740Q_full.json.gz",
    ),
    "T09_paintomics_carbon.tsv": (
        "scripts/06_paintomics_tables.py (job oK2h66740Q)",
        "The carbon pathways the model made claims about, significant or not",
        "results/omics/paintomics_raw/job1_oK2h66740Q_full.json.gz",
    ),
    "T10_metabolite_hubs.tsv": (
        "scripts/06_paintomics_tables.py (job oK2h66740Q)",
        "PaintOmics metabolite hub analysis: predicted compounds ranked by REAL DE genes nearby",
        "results/omics/paintomics_raw/job1_oK2h66740Q_full.json.gz",
    ),
    "T11_hardware_ladder.tsv": (
        "scripts/07_hardware_ladder.py",
        "Six OSDR studies x ten gene sets: median shift, Mann-Whitney p, model prediction",
        "data/study_registry.tsv + GeneLab DGE tables",
    ),
    "T12_ladder_contrasts.tsv": (
        "scripts/07_hardware_ladder.py",
        "Set-vs-set separations per study — the illumination test and the enclosure test",
        "data/study_registry.tsv + GeneLab DGE tables",
    ),
    "T13_light_sensitivity.tsv": (
        "scripts/09_light_sensitivity.py",
        "Table 3 operating points, the lit-canister compensation point, and the microgravity deficit at Q = 1000 and Q = 60",
        "LunarLeaf-CFD T13 + scripts/fvcb.py",
    ),
    "T14_physiospace_scores.csv": (
        "scripts/08_physiospace_projection.R",
        "PhysioSpace stress pattern projection scores across OSDR studies (optional)",
        "data/study_registry.tsv + PlantPhysioSpace",
    ),
    "T15_enclosure_sensitivity.tsv": (
        "scripts/enclosure_sensitivity.py",
        "0-D enclosure drawdown and carbon-gain sensitivity across headspace, leaf area, permeance and light",
        "0-D mass balance model + LunarLeaf-CFD",
    ),
    "T16_independent_pathways.tsv": (
        "scripts/independent_pathway_test.py",
        "Independent pathway enrichment and shift tests outside PaintOmics",
        "T01_osd522_transcriptome.tsv + KEGG",
    ),
    "T17_integration_control_nometa.tsv": (
        "scripts/control_integration_nometa.py",
        "PaintOmics pathway integration control using RNA-seq + Proteomics only (no predicted metabolites)",
        "results/omics/paintomics_raw/job1_oK2h66740Q_full.json.gz",
    ),
}

ALL_DESCRIBED = {**CFD_DESCRIBED, **OMICS_DESCRIBED}


def write_dir_manifest(folder: str, catalog: dict[str, tuple[str, str, str]]) -> list[str]:
    present = sorted(
        f for f in os.listdir(folder)
        if f.endswith((".tsv", ".csv")) and f != "MANIFEST.tsv"
    )
    missing = [f for f in present if f not in catalog]
    out = os.path.join(folder, "MANIFEST.tsv")
    with open(out, "w") as fh:
        fh.write("# Every derived table in this directory, what it holds, and what made it.\n")
        fh.write("# Regenerate with scripts/make_manifest.py, which fails on an undescribed file.\n")
        fh.write("file\tbytes\tgenerating_script\tdescription\tsource_data\n")
        for f in present:
            script, desc, src = catalog.get(f, ("?", "UNDESCRIBED", "?"))
            size = os.path.getsize(os.path.join(folder, f))
            fh.write(f"{f}\t{size}\t{script}\t{desc}\t{src}\n")
    print(f"wrote {out}  ({len(present)} tables)")
    return missing


def main() -> int:
    missing_cfd = write_dir_manifest(CFD_TABLES, CFD_DESCRIBED)
    missing_omics = write_dir_manifest(OMICS_TABLES, OMICS_DESCRIBED)
    missing_all = write_dir_manifest(TABLES, ALL_DESCRIBED)

    # Also write top-level results/MANIFEST.tsv
    top_out = os.path.join(RESULTS, "MANIFEST.tsv")
    with open(os.path.join(TABLES, "MANIFEST.tsv")) as src_fh, open(top_out, "w") as dst_fh:
        dst_fh.write(src_fh.read())
    print(f"wrote {top_out}")

    all_missing = sorted(set(missing_cfd + missing_omics + missing_all))
    if all_missing:
        print(f"  FAIL: undescribed tables: {', '.join(all_missing)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
