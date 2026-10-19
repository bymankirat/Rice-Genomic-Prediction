"""Run the rice genomic prediction analysis.

    python run_all.py            # everything (needs internet for step 1)
    python run_all.py --from 2   # from the files already downloaded
"""
import argparse
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "pipeline"))

STEPS = [
    ("s01_fetch", "download the Rice Diversity Panel genotypes and phenotypes"),
    ("s02_genotypes", "read PLINK files, quality control, genotype matrix and traits"),
    ("s03_structure", "population structure: PCA and genetic groups"),
    ("s04_gwas", "GWAS per trait with structure correction, checked on known genes"),
    ("s05_prediction", "genomic prediction: within and across genetic groups"),
    ("s06_export", "write docs/results.js"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", type=int, default=1)
    args = ap.parse_args()
    t0 = time.time()
    for i, (mod, what) in enumerate(STEPS, 1):
        if i < args.start:
            continue
        print(f"\n=== step {i}: {what} ===", flush=True)
        importlib.import_module(mod).main()
    print(f"\nfinished in {(time.time() - t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
