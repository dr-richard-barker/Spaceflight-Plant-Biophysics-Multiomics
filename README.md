# Spaceflight-Plant-Biophysics-Multiomics

**From microgravity boundary-layer biophysics to the spaceflight multi-ome: a coupled lattice-Boltzmann and Farquhar–von Caemmerer–Berry model of *Arabidopsis* gas exchange and enclosure carbon starvation**

**Author:** Richard J. Barker ([ORCID 0000-0001-5681-9857](https://orcid.org/0000-0001-5681-9857)), Purdue University, West Lafayette, IN, USA (`dr.richard.barker@gmail.com`)  
**Target Special Issue:** *Progress in Biophysics and Molecular Biology* — Special Issue: [*Biophysics and Molecular Biology for Space Flight*](https://www.sciencedirect.com/special-issue/337891/biophysics-and-molecular-biology-for-space-flight) (SI 337891)  
**Licenses:** Code: **MIT** · Manuscript, figures & derived tables: **CC-BY-4.0** (`LICENSE`)

---

## Live Interactive Suite & Compiled Manuscripts

| Resource | Path / URL | Description |
| :--- | :--- | :--- |
| **Main Manuscript (PDF)** | [`manuscript/main.pdf`](manuscript/main.pdf) | Compiled 23-page Elsevier `elsarticle` manuscript (10 figures, 5 tables) |
| **Supplementary Information (PDF)** | [`manuscript/supplementary.pdf`](manuscript/supplementary.pdf) | Compiled 12-page Supplementary Information (Tables S1–S10, Figures S1–S8) |
| **Word Manuscript (`.docx`)** | [`manuscript/Spaceflight_Biophysics_Multiomics_Manuscript.docx`](manuscript/Spaceflight_Biophysics_Multiomics_Manuscript.docx) | Editable Word version of the integrated manuscript |
| **Markdown Manuscript** | [`manuscript/MANUSCRIPT.md`](manuscript/MANUSCRIPT.md) | Complete GitHub-flavored Markdown edition of the manuscript |
| **Live 2-D LBM CFD Workbench** | [`index.html`](index.html) (`npm run dev`) | In-browser D2Q9/D2Q5 lattice-Boltzmann solver across 26 scenarios ($1\,g \to \text{Mars} \to \text{Moon} \to \mu g$; leaf, rosette, canopy; BRIC, CARA, VEGGIE) |
| **Interactive FvCB & Multi-Omics Explorer** | [`public/explorer/index.html`](public/explorer/index.html) | Browser FvCB $\mathrm{C_3}$ calculator, OSD-522 volcano/scatter browser, PaintOmics 4 viewer, and 6-study OSDR hardware ladder |
| **Web Manuscript Landing Page** | [`public/paper.html`](public/paper.html) | Illustrated HTML overview linking all figures, tables, and live tools |

---

## Scientific Summary

Plants grown in spaceflight exhibit conserved transcriptional and proteomic stress responses whose physical driver has remained unresolved between two candidate transport bottlenecks:
1. **The microgravity leaf boundary layer ($\sim 0.7\text{–}11\,\mathrm{mm}$, $C_a \to C_s$)**: Loss of buoyancy-driven natural convection in microgravity ($\mu g$) collapses external gas transport to molecular diffusion. Our benchmark-validated 2-D D2Q9/D2Q5 lattice-Boltzmann CFD solver shows that removing gravity eliminates convective plumes ($u_{\max}\propto\sqrt{g}\to 0$) and reduces steady-state boundary-layer conductance ($g_{bl}$) from **$0.997 \to 0.443\,\mathrm{mol\,m^{-2}\,s^{-1}}$** for an isolated leaf ($\delta: 0.67 \to 1.50\,\mathrm{mm}$), **$0.540 \to 0.272\,\mathrm{mol\,m^{-2}\,s^{-1}}$** for an overlapping rosette ($\delta: 1.23 \to 2.45\,\mathrm{mm}$), and **$0.133 \to 0.058\,\mathrm{mol\,m^{-2}\,s^{-1}}$** for a dense microgreen canopy ($\delta: 5.02 \to 11.42\,\mathrm{mm}$). Forced ventilation reverses this penalty, requiring **$\approx 2.0\,\mathrm{cm\,s^{-1}}$** for an isolated leaf and **$\approx 6.1\,\mathrm{cm\,s^{-1}}$** for a rosette.
2. **Growth-hardware headspace drawdown ($10^1\text{–}10^2\,\mathrm{mm}$, $C_{\text{cabin}} \to C_a$)**: Coupling the CFD conductances and zero-dimensional (0-D) canister mass balances for sealed **BRIC**, micropore-taped **CARA**, and actively ventilated **VEGGIE** hardware to a Farquhar–von Caemmerer–Berry (FvCB) $\mathrm{C_3}$ model demonstrates that an illuminated sealed BRIC canister exhausts its $\sim 400\,\mathrm{ppm}$ headspace $\mathrm{CO_2}$ to the compensation point within **$\sim 7\,\mathrm{min}$**, collapsing $12\,\mathrm{h}$ carbon gain to **$\approx 1\%$** of an unenclosed Earth control (versus **$\approx 90\%$** in CARA and **$\approx 100\%$** in VEGGIE) and raising Rubisco oxygenation ($\varphi = V_o/(V_o+V_c)$) from **$18.0\%$ to $46.2\%$** at $60\,\mu\mathrm{mol\,m^{-2}\,s^{-1}}$. Because canister drawdown is a mass balance, it occurs symmetrically in both flight and ground controls, cancelling out of the flight-versus-ground differential contrast while shifting both arms into severe carbon starvation (an enclosure effect **$22\times$ larger** than the microgravity boundary-layer effect) and causing the $1\,g$ and $\mu g$ boundary layers to **converge** as net assimilation $A\to 0$.
3. **Empirical falsification across NASA OSDR multi-omics**: After an exhaustive census of all 567 studies in NASA OSDR revealed **0 of 12 metabolomics studies are in plants** (`T03`), we formulated a falsifiable 21-compound flux-propagated metabolite scaffold (`T06`, explicitly labelled as model output) and tested the coupled model against OSD-522 (BRIC-LED) shoot RNA-seq and proteomics (`T01`, `T02`, `T07`) and across six matched flight-versus-ground contrasts from five OSDR studies (`OSD-522`, `OSD-38`, `OSD-321`, `OSD-678`, `OSD-427`; `T11`, `T12`). The model called **8 of 9 transcriptome gene-set shifts** and **4 of 7 in the proteome**, while PaintOmics 4 integration (`T08`, `T09`) revealed that **endoplasmic reticulum unfolded protein response (`ath04141`, $p = 8.6\times 10^{-9}$)** dominates globally and photorespiration ranks in the bottom $8\%$ (**214th of 232 pathways**). Across the six OSDR contrasts, illumination completely separated all lit studies from all dark studies on carbon starvation over photosynthesis ($p = 0.05$).

---

## Repository Structure

```text
Spaceflight-Plant-Biophysics-Multiomics/
├── manuscript/
│   ├── main.tex / main.pdf                    # Elsevier elsarticle main manuscript (23 pp., 10 figs, 5 tables)
│   ├── supplementary.tex / supplementary.pdf  # Supplementary Information (12 pp., Tables S1–S10, Figs S1–S8)
│   ├── MANUSCRIPT.md                          # Complete Markdown manuscript
│   ├── Spaceflight_Biophysics_Multiomics_Manuscript.docx  # Word manuscript
│   ├── references.bib                         # Verified BibTeX bibliography (23 references)
│   └── figures/                               # All vector PDF and high-res PNG figures
├── src/                                       # TypeScript / React 2D D2Q9 + D2Q5 LBM CFD solver & UI
│   ├── solver/lbm/                            # LBMFluid.ts (BGK/MRT + Guo forcing), ScalarField.ts, buoyancy.ts
│   ├── leaf/                                  # Leaf, rosette, canopy, and real BRIC/CARA dish geometries
│   └── ui/                                    # Interactive browser CFD workbench
├── validation/                                # Headless numerical verification & physical sweeps
│   ├── gates.ts / mrt_check.ts                # 4 numerical validation benchmarks (Ghia, Strouhal, de Vahl Davis, erfc)
│   ├── sweep.ts / stationarity_check.ts       # 26-scenario steady-state sweep (150k/900k steps)
│   ├── resolution_ladder.ts / domain_size.ts  # 3-grid Celik (2008) GCI & 1x/1.5x/2x domain-size sweeps
│   ├── dish_hardware.ts / dish_feedback.ts    # Real 60mm BRIC & 120mm CARA dish cross-sections + RH saturation
│   └── feedback_0d.py / analyze_model.py      # 0-D photoperiod carbon gain & physical unit calibration
├── scripts/                                   # Python FvCB biochemical model & NASA OSDR multi-omics pipeline
│   ├── fvcb.py                                # Coupled CFD -> FvCB C3 model (Bernacchi et al. 2001, 2003)
│   ├── 00_fetch_osdr.py .. 08_make_figures.py # End-to-end OSDR RNA-seq (PyDESeq2), proteomics, PaintOmics & figures
│   ├── enclosure_sensitivity.py               # 0-D enclosure sensitivity across V, A_leaf, k_seal, Q, and cabin CO2
│   ├── independent_pathway_test.py            # Metabolite-free Fisher & Mann-Whitney pathway verification
│   ├── check_js_parity.py / check_site_data.py# Automated JS/Python FvCB parity (1e-6) & site-data consistency checks
│   └── make_manifest.py                       # Generates machine-verified MANIFEST.tsv across all 45 tables
├── data/                                      # Primary & cached reference inputs (with SHA-256 provenance)
│   ├── study_registry.tsv                     # Pre-specified inclusion/exclusion audit of 11 OSDR plant studies
│   ├── kegg_cache/ & mapman/                  # Cached KEGG REST pathways & GoMapMan/MapManStore definitions
│   └── paintomics_inputs/                     # Exact inputs submitted to PaintOmics 4 (job oK2h66740Q)
├── results/
│   ├── MANIFEST.tsv                           # Master manifest of all 45 derived tables (29 CFD + 16 Omics)
│   ├── cfd/tables/ & cfd/figures/             # CFD tables T1–T24 and figures F1–F17
│   ├── omics/tables/ & omics/figures/         # Omics/FvCB tables T01–T17 and figures fig01–fig08, figS1
│   └── tables/ & figures/                     # Unified symlinks/copies for seamless script execution
├── public/                                    # Static assets bundled into the GitHub Pages deployment
│   ├── paper.html                             # Unified illustrated web manuscript
│   └── explorer/                              # Interactive FvCB & Multi-Omics web explorer
├── .zenodo.json / CITATION.cff / LICENSE      # Zenodo deposition metadata, citation file & MIT/CC-BY-4.0 licenses
└── package.json / requirements.txt            # Node.js and Python dependency manifests
```

---

## Reproducing Everything From Scratch

### 1. Lattice-Boltzmann CFD Benchmarks, Sweeps & Web Build (Node.js $\ge 20$)

```bash
npm ci
npm run gates        # Run the 4 numerical validation gates (lid cavity, wake, Ra=10^4 cavity, erfc)
npm run mrt-check    # Verify BGK vs MRT collision operator parity across all 4 gates
npm run build        # Typecheck and build the interactive web suite into dist/
```

### 2. Coupled FvCB Model, Multi-Omics Pipeline & Verification Checks (Python $\ge 3.10$)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Verify FvCB biochemical self-checks, JS/Python parity, web data integrity, and table manifests:
python3 scripts/fvcb.py
python3 scripts/check_js_parity.py
python3 scripts/check_site_data.py
python3 scripts/make_manifest.py
```

### 3. Compiling the Elsevier LaTeX Manuscript & Supplementary PDF

```bash
cd manuscript
latexmk -pdf supplementary.tex
latexmk -pdf main.tex
latexmk -pdf supplementary.tex
```

---

## FAIR Compliance & Scientific Integrity Statement

- **Findable & Accessible**: Complete deposition metadata in [`.zenodo.json`](.zenodo.json) and [`CITATION.cff`](CITATION.cff), linking all five primary NASA OSDR study accessions (`OSD-522`, `OSD-38`, `OSD-321`, `OSD-678`, `OSD-427`) and the Edinburgh DataShare gas-exchange dataset (`doi:10.7488/ds/3453`).
- **Interoperable & Reusable**: All 45 derived `.csv`/`.tsv` tables are documented with row counts, column counts, generating scripts, and descriptions in [`results/MANIFEST.tsv`](results/MANIFEST.tsv).
- **Explicit Distinction Between Measured and Predicted Data**: Because NASA OSDR currently contains zero plant metabolomics datasets (`T03_osdr_metabolomics_census.tsv`: 0 of 12 studies), all 21 compound values in `T06_predicted_compounds.tsv` are **quasi-steady-state FvCB flux predictions** and are explicitly labelled as model output across every table, figure, and web interface.
