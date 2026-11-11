#!/usr/bin/env python3
"""Central path configuration for the unified Spaceflight-Plant-Biophysics-Multiomics repository.

Every Python script imports from here so that the CFD results (`results/cfd/`),
multi-omics results (`results/omics/`), combined tables (`results/tables/`),
manuscript (`manuscript/`), and interactive web suite (`public/explorer/`)
share a single source of truth.
"""

from __future__ import annotations

import os

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(SCRIPTS)

DATA = os.path.join(REPO, "data")
CACHE = os.path.join(DATA, "cache")
UPLOAD = os.path.join(DATA, "paintomics_upload")
STUDY_REGISTRY = os.path.join(DATA, "study_registry.tsv")
CFD_SWEEP = os.path.join(REPO, "results", "cfd", "tables", "T13_boundary_layer.csv")

RESULTS = os.path.join(REPO, "results")
CFD_RESULTS = os.path.join(RESULTS, "cfd")
CFD_TABLES = os.path.join(CFD_RESULTS, "tables")
CFD_FIGURES = os.path.join(CFD_RESULTS, "figures")

OMICS_RESULTS = os.path.join(RESULTS, "omics")
OMICS_TABLES = os.path.join(OMICS_RESULTS, "tables")
OMICS_FIGURES = os.path.join(OMICS_RESULTS, "figures")

TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(RESULTS, "figures")
SUPP_FIGURES = os.path.join(FIGURES, "supplementary")
PAINTOMICS_RAW = os.path.join(OMICS_RESULTS, "paintomics_raw")

MANUSCRIPT = os.path.join(REPO, "manuscript")
DOCS = os.path.join(REPO, "public", "explorer")
DOCS_DATA = os.path.join(DOCS, "data")
METHODS = os.path.join(REPO, "methods")

# Direct in-repo path to the validated LunarLeaf-CFD tables.
LUNARLEAF = CFD_TABLES


def ensure(*dirs: str) -> None:
    """Create output directories on demand."""
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def table(name: str) -> str:
    return os.path.join(TABLES, name)


def figure(name: str) -> str:
    return os.path.join(FIGURES, name)
