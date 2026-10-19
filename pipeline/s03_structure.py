"""Step 3: population structure.

Rice is divided into a few genetic groups (indica, aus, the japonicas, aromatic). They differ in both markers and
traits, so structure has to be measured first: it is the main confounder in GWAS, and it decides how honest a
prediction test is. Principal components on the marker matrix summarise the structure; k-means on the leading
components gives discrete groups for the across-group prediction test in step 5.
"""
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from util import DERIVED, RESULTS, config, log, save_json

N_PC = 10
K_RANGE = range(2, 9)


def main():
    g = np.load(DERIVED / "genotypes.npy")
    samples = pd.read_csv(DERIVED / "samples.csv")
    z = (g - g.mean(axis=0)) / np.where(g.std(axis=0) == 0, 1, g.std(axis=0))

    pca = PCA(n_components=N_PC, random_state=config.SEED).fit(z)
    pcs = pca.transform(z)
    var = pca.explained_variance_ratio_
    log("variance explained by the first PCs: " + ", ".join(f"PC{i+1} {v:.1%}" for i, v in enumerate(var[:5])))

    scores = {}
    for k in K_RANGE:
        lab = KMeans(n_clusters=k, n_init=20, random_state=config.SEED).fit_predict(pcs[:, :5])
        scores[k] = float(silhouette_score(pcs[:, :5], lab))
    best_k = max(scores, key=scores.get)
    labels = KMeans(n_clusters=best_k, n_init=50, random_state=config.SEED).fit_predict(pcs[:, :5])
    log(f"silhouette by k: { {k: round(v, 3) for k, v in scores.items()} } -> k = {best_k}")

    out = samples.copy()
    for i in range(N_PC):
        out[f"PC{i+1}"] = pcs[:, i]
    out["group"] = [f"G{l+1}" for l in labels]
    out.to_csv(DERIVED / "structure.csv", index=False)
    sizes = out.group.value_counts().sort_index()
    log("group sizes:\n" + sizes.to_string())
    save_json(RESULTS / "structure.json", {
        "n_pcs": N_PC, "variance_explained": [float(v) for v in var],
        "silhouette": scores, "k": int(best_k),
        "group_sizes": {str(k): int(v) for k, v in sizes.items()},
    })


if __name__ == "__main__":
    main()
