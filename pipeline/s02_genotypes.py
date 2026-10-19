"""Step 2: read the PLINK genotypes, run quality control and build a numeric genotype matrix.

A .ped file holds two alleles per marker per line. Each marker is turned into the count of its minor allele
(0, 1 or 2), which is what both GWAS and genomic prediction need. Quality control follows the usual thresholds:
markers with a low call rate or rare minor allele, and samples with many missing calls, are dropped.
"""
import gzip
import zipfile

import numpy as np
import pandas as pd

from util import DERIVED, RAW, RESULTS, log, save_json

SNP_CALL_RATE = 0.95
SAMPLE_CALL_RATE = 0.90
MAF = 0.05


def find(ext):
    hits = sorted(RAW.rglob(f"*{ext}"))
    if not hits:
        raise FileNotFoundError(f"no {ext} file under {RAW}")
    return max(hits, key=lambda p: p.stat().st_size)


def read_map(path):
    m = pd.read_csv(path, sep=r"\s+", header=None, names=["chrom", "snp", "cm", "pos"], dtype={"snp": str})
    return m


def read_ped(path, n_snps):
    """-> (sample table, allele array of shape (n_samples, n_snps, 2))."""
    samples, alleles = [], []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            samples.append({"fid": f[0], "iid": f[1], "sex": f[4], "pheno": f[5]})
            g = np.array(f[6:], dtype="U1")
            if g.size != 2 * n_snps:
                raise ValueError(f"{samples[-1]['iid']}: {g.size} allele fields, expected {2 * n_snps}")
            alleles.append(g.reshape(n_snps, 2))
    return pd.DataFrame(samples), np.stack(alleles)


def encode(alleles):
    """Count of the minor allele per genotype; -1 marks missing."""
    n_samples, n_snps, _ = alleles.shape
    out = np.full((n_samples, n_snps), -1, dtype=np.int8)
    minor = np.empty(n_snps, dtype="U1")
    major = np.empty(n_snps, dtype="U1")
    missing_codes = {"0", "N", "-", "."}
    for j in range(n_snps):
        col = alleles[:, j, :]
        called = ~np.isin(col, list(missing_codes))
        both = called.all(axis=1)
        vals, counts = np.unique(col[both], return_counts=True)
        if len(vals) == 0:
            continue
        order = np.argsort(-counts)
        major[j] = vals[order[0]]
        minor[j] = vals[order[1]] if len(vals) > 1 else vals[order[0]]
        out[both, j] = (col[both] == minor[j]).sum(axis=1)
    return out, minor, major


def read_phenotypes(samples):
    """PLINK phenotype table -> one row per kept sample, numeric traits only."""
    path = RAW / "phenotypes.txt"
    df = pd.read_csv(path, sep="\t" if "\t" in path.read_text(errors="ignore")[:5000] else r"\s+",
                     na_values=["NA", "", "-", "."])
    cols = list(df.columns)
    log(f"phenotype file: {df.shape[0]} rows, {df.shape[1]} columns; first columns {cols[:3]}")
    id_cols = cols[:2] if df[cols[1]].dtype == object or df[cols[1]].astype(str).str.isalpha().any() else cols[:1]
    df = df.rename(columns={id_cols[0]: "fid", **({id_cols[1]: "iid"} if len(id_cols) > 1 else {})})
    if "iid" not in df.columns:
        df["iid"] = df["fid"]
    traits = [c for c in df.columns if c not in ("fid", "iid")]
    for c in traits:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    keep = [c for c in traits if df[c].notna().sum() >= 100]
    log(f"{len(keep)} of {len(traits)} traits have at least 100 measured varieties")
    out = samples[["fid", "iid"]].merge(df[["fid", "iid"] + keep], on=["fid", "iid"], how="left")
    if out[keep].notna().sum().sum() == 0:  # ids may differ in the first column only
        out = samples[["fid", "iid"]].merge(df[["iid"] + keep], on="iid", how="left")
    log(f"phenotypes matched for {int(out[keep].notna().any(axis=1).sum())} of {len(samples)} genotyped varieties")
    return out


def main():
    mp = read_map(find(".map"))
    log(f"{len(mp)} markers in the map file, chromosomes {sorted(mp.chrom.unique())}")
    samples, alleles = read_ped(find(".ped"), len(mp))
    log(f"{len(samples)} samples read")
    geno, minor, major = encode(alleles)
    del alleles

    missing = geno < 0
    snp_call = 1 - missing.mean(axis=0)
    sample_call = 1 - missing.mean(axis=1)
    g = np.where(missing, np.nan, geno).astype(np.float32)
    freq = np.nanmean(g, axis=0) / 2
    maf = np.minimum(freq, 1 - freq)

    keep_s = sample_call >= SAMPLE_CALL_RATE
    keep_m = (snp_call >= SNP_CALL_RATE) & (maf >= MAF)
    log(f"QC: {keep_s.sum()}/{len(keep_s)} samples and {keep_m.sum()}/{len(keep_m)} markers pass "
        f"(call rate >= {SNP_CALL_RATE}, MAF >= {MAF})")

    g = g[np.ix_(keep_s, keep_m)]
    col_mean = np.nanmean(g, axis=0)
    idx = np.where(np.isnan(g))
    g[idx] = np.take(col_mean, idx[1])  # mean imputation, the usual choice for genomic prediction

    np.save(DERIVED / "genotypes.npy", g.astype(np.float32))
    kept_samples = samples[keep_s].assign(call_rate=sample_call[keep_s]).reset_index(drop=True)
    kept_samples.to_csv(DERIVED / "samples.csv", index=False)
    read_phenotypes(kept_samples).to_csv(DERIVED / "phenotypes.csv", index=False)
    mp[keep_m].assign(maf=maf[keep_m], call_rate=snp_call[keep_m],
                      minor=minor[keep_m], major=major[keep_m]).to_csv(DERIVED / "markers.csv", index=False)
    save_json(RESULTS / "qc.json", {
        "samples_in": int(len(keep_s)), "samples_kept": int(keep_s.sum()),
        "markers_in": int(len(keep_m)), "markers_kept": int(keep_m.sum()),
        "thresholds": {"snp_call_rate": SNP_CALL_RATE, "sample_call_rate": SAMPLE_CALL_RATE, "maf": MAF},
        "missing_rate_before_qc": float(missing.mean()),
        "chromosomes": int(mp[keep_m].chrom.nunique()),
    })
    log(f"genotype matrix: {g.shape[0]} samples x {g.shape[1]} markers")


if __name__ == "__main__":
    main()
