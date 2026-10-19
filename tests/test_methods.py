"""Checks on the methods themselves, on data made up for the test."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))


def test_genotype_encoding_counts_the_minor_allele():
    from s02_genotypes import encode
    # marker 1: A common, G rare; marker 2: everything missing for one sample
    alleles = np.array([[["A", "A"], ["C", "C"]],
                        [["A", "G"], ["C", "T"]],
                        [["G", "G"], ["0", "0"]],
                        [["A", "A"], ["T", "T"]]])
    g, minor, major = encode(alleles)
    assert list(g[:, 0]) == [0, 1, 2, 0], g[:, 0]
    assert minor[0] == "G" and major[0] == "A"
    assert g[2, 1] == -1, "missing genotypes must stay flagged"


def test_gwas_finds_a_planted_effect_and_ignores_noise():
    from s04_gwas import gwas
    rng = np.random.default_rng(0)
    n, m = 300, 200
    g = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    cov = rng.normal(size=(n, 2))
    y = 2.5 * g[:, 7] + cov @ np.array([1.0, -0.5]) + rng.normal(scale=1.0, size=n)
    beta, se, p = gwas(y, g, cov)
    assert int(np.argmin(p)) == 7, "the marker with the real effect should win"
    assert p[7] < 1e-10
    assert abs(beta[7] - 2.5) < 0.3
    assert (p[np.arange(m) != 7] < 0.05 / m).sum() == 0, "no false positives at the Bonferroni threshold"


def test_gwas_covariates_remove_a_confounder():
    from s04_gwas import gwas
    rng = np.random.default_rng(1)
    n, m = 400, 150
    group = rng.integers(0, 2, n)
    # markers differ between groups, and so does the trait: a classic structure confounder
    g = rng.binomial(2, np.where(group[:, None] == 1, 0.7, 0.2) * np.ones((n, m)))
    y = 3.0 * group + rng.normal(size=n)
    _, _, p_naive = gwas(y, g.astype(float), np.zeros((n, 1)))
    _, _, p_adj = gwas(y, g.astype(float), group.reshape(-1, 1).astype(float))
    assert (p_naive < 0.05 / m).mean() > 0.3, "without correction, structure should create false hits"
    assert (p_adj < 0.05 / m).mean() < 0.05, "the covariate should remove them"


def test_accuracy_handles_missing_and_flat_predictions():
    from s05_prediction import accuracy
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, np.nan])
    r, rmse = accuracy(y, np.array([1.1, 2.2, 2.9, 4.1, 5.2, 5.8, 99.0]))
    assert r > 0.99 and rmse < 0.3, (r, rmse)
    r2, _ = accuracy(y, np.full(7, 2.0))
    assert np.isnan(r2), "a constant prediction has no correlation"
    r3, _ = accuracy(y, np.full(7, np.nan))
    assert np.isnan(r3), "too few usable points must give no score"
