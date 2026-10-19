"""Checks on the saved results: the numbers in the report must match the files the pipeline wrote."""
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DERIVED = ROOT / "data" / "derived"


def load(name):
    p = RESULTS / name
    if not p.exists():
        pytest.skip(f"{name} not produced yet")
    return json.loads(p.read_text())


def test_qc_keeps_most_markers_and_samples():
    qc = load("qc.json")
    assert qc["samples_kept"] >= 350, qc
    assert qc["markers_kept"] >= 20000, qc
    assert qc["markers_kept"] <= qc["markers_in"]
    assert qc["chromosomes"] == 12, "rice has 12 chromosomes"


def test_structure_finds_several_groups():
    s = load("structure.json")
    assert 2 <= s["k"] <= 8
    assert sum(s["group_sizes"].values()) == load("qc.json")["samples_kept"]
    assert s["variance_explained"][0] > 0.05, "PC1 should carry real structure"


def test_gwas_recovers_a_known_gene():
    g = load("gwas.json")
    checks = g["known_gene_checks"]
    if not checks:
        pytest.skip("no known-gene checks were possible")
    assert any(c["significant"] for c in checks), "no textbook gene recovered; the GWAS is suspect"


def test_gwas_is_not_inflated():
    g = load("gwas.json")
    lam = [t["lambda_gc"] for t in g["traits"] if t["lambda_gc"] == t["lambda_gc"]]
    assert lam, "no traits analysed"
    assert max(lam) < 2.0, f"genomic inflation too high: {max(lam):.2f}"


def test_prediction_beats_the_mean_baseline():
    p = load("prediction.json")
    m = p["mean_r_random"]
    assert m["ridge"] > 0.3, m
    assert m["ridge"] > m.get("ridge (100 markers)", 0), "all markers should beat 100 markers"


def test_across_group_is_harder_than_random():
    p = load("prediction.json")
    if not p.get("mean_r_across_group"):
        pytest.skip("no across-group test")
    assert p["mean_r_across_group"]["ridge"] < p["mean_r_random"]["ridge"], \
        "predicting an unseen genetic group should be harder, not easier"


def test_tables_agree_with_summary():
    if not (RESULTS / "prediction_summary.csv").exists():
        pytest.skip("not produced yet")
    s = pd.read_csv(RESULTS / "prediction_summary.csv")
    folds = pd.read_csv(RESULTS / "prediction_folds.csv")
    row = s.iloc[0]
    sub = folds[(folds.trait == row.trait) & (folds.test == row.test) & (folds.model == row.model)]
    assert abs(sub.r.mean() - row.r) < 1e-6
