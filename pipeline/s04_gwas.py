"""Step 4: GWAS for each trait, with population structure as a covariate, checked against genes already known.

Every marker is tested one at a time in a linear model, trait ~ marker + leading principal components. The PCs
absorb the differences between rice groups, which otherwise produce associations across the whole genome.

The point of the step is not to find new genes but to show the analysis behaves: for traits whose genes are
textbook (plant height and sd1, amylose content and Waxy, grain length and GS3), the strongest signal should land
on the right chromosome near the right position.
"""
import numpy as np
import pandas as pd
from scipy import stats

from util import DERIVED, RESULTS, config, log, save_json

N_PC = 5
# Published major genes used as positive controls (MSU/RAP positions on the Nipponbare reference, Mb)
KNOWN = {
    "Plant.height": [("sd1", 1, 38.4)],
    "Flowering.time.at.Arkansas": [("Hd1", 6, 9.3)],
    "Amylose.content": [("Waxy (GBSSI)", 6, 1.77)],
    "Seed.length": [("GS3", 3, 16.73)],
    "Seed.width": [("GW5/qSW5", 5, 5.36)],
    "Seed.volume": [("GS3", 3, 16.73)],
}
WINDOW_MB = 1.0


def gwas(y, g, cov):
    """Association of every marker with y, controlling for cov. Returns beta, se, p."""
    ok = ~np.isnan(y)
    y, g = y[ok], g[ok]
    x0 = np.column_stack([np.ones(ok.sum()), cov[ok]])
    # project trait and markers onto the space orthogonal to the covariates
    q, _ = np.linalg.qr(x0)
    yr = y - q @ (q.T @ y)
    gr = g - q @ (q.T @ g)
    ss = (gr ** 2).sum(axis=0)
    ss[ss == 0] = np.nan
    beta = (gr * yr[:, None]).sum(axis=0) / ss
    resid_ss = ((yr[:, None] - gr * beta) ** 2).sum(axis=0)
    df = len(y) - x0.shape[1] - 1
    se = np.sqrt(resid_ss / df / ss)
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return beta, se, p


def main():
    g = np.load(DERIVED / "genotypes.npy")
    markers = pd.read_csv(DERIVED / "markers.csv")
    struct = pd.read_csv(DERIVED / "structure.csv")
    phen = pd.read_csv(DERIVED / "phenotypes.csv")
    cov = struct[[f"PC{i+1}" for i in range(N_PC)]].to_numpy()
    traits = [c for c in phen.columns if c not in ("fid", "iid")]
    bonferroni = 0.05 / len(markers)
    log(f"{len(traits)} traits, {len(markers)} markers, genome-wide threshold p < {bonferroni:.2e}")

    hits, checks = [], []
    for trait in traits:
        y = phen[trait].to_numpy(dtype=float)
        if np.isfinite(y).sum() < 100:
            continue
        beta, se, p = gwas(y, g, cov)
        res = markers.assign(beta=beta, se=se, p=p)
        sig = res[res.p < bonferroni].sort_values("p")
        lam = float(np.median(stats.chi2.isf(p[np.isfinite(p)], 1)) / stats.chi2.ppf(0.5, 1))
        top = res.loc[res.p.idxmin()]
        hits.append({"trait": trait, "n": int(np.isfinite(y).sum()), "n_significant": int(len(sig)),
                     "lambda_gc": lam, "top_snp": top.snp, "top_chrom": int(top.chrom),
                     "top_pos_mb": float(top.pos) / 1e6, "top_p": float(top.p)})
        sig.head(25).assign(trait=trait).to_csv(RESULTS / f"gwas_{trait.replace('.', '_')}_top.csv", index=False)

        for gene, chrom, pos_mb in KNOWN.get(trait, []):
            near = res[(res.chrom == chrom) & (abs(res.pos / 1e6 - pos_mb) <= WINDOW_MB)]
            best = near.loc[near.p.idxmin()] if len(near) else None
            checks.append({"trait": trait, "gene": gene, "chrom": chrom, "expected_mb": pos_mb,
                           "markers_in_window": int(len(near)),
                           "best_p_in_window": float(best.p) if best is not None else None,
                           "best_pos_mb": float(best.pos) / 1e6 if best is not None else None,
                           "significant": bool(best is not None and best.p < bonferroni),
                           "distance_of_top_hit_mb": float(abs(top.pos / 1e6 - pos_mb)) if int(top.chrom) == chrom else None})
    pd.DataFrame(hits).to_csv(RESULTS / "gwas_summary.csv", index=False)
    pd.DataFrame(checks).to_csv(RESULTS / "gwas_known_genes.csv", index=False)
    save_json(RESULTS / "gwas.json", {"threshold": bonferroni, "n_markers": int(len(markers)),
                                      "n_pcs": N_PC, "traits": hits, "known_gene_checks": checks})
    ok = sum(c["significant"] for c in checks)
    log(f"known-gene checks: {ok}/{len(checks)} recovered a genome-wide significant marker in the window")


if __name__ == "__main__":
    main()
