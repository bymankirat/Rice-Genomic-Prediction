"""Step 5: genomic prediction - can a model trained on some varieties predict the traits of others?

This is the question a breeding programme actually asks. Three tests, in increasing difficulty:

  1. random 5-fold cross-validation, repeated 5 times: the easy case, where training and test varieties are mixed;
  2. leave-one-group-out: train on all genetic groups but one and predict the group that was held out, the honest
     test of whether a model moves between rice types;
  3. a marker-free baseline (predict the training mean) and a 100-marker baseline, so the gain from using all
     markers is visible rather than assumed.

Models: ridge regression on all markers (equivalent to GBLUP with a marker kernel) and gradient boosting, which
can pick up non-additive effects. Accuracy is the correlation between predicted and observed values, the measure
breeders use.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

from util import DERIVED, RESULTS, config, log, save_json

REPEATS, FOLDS = 5, 5
ALPHAS = np.logspace(1, 5, 17)
MIN_GROUP = 20


def accuracy(y_true, y_pred):
    ok = np.isfinite(y_true) & np.isfinite(y_pred)
    if ok.sum() < 5 or np.std(y_pred[ok]) == 0:
        return np.nan, np.nan
    r = float(np.corrcoef(y_true[ok], y_pred[ok])[0, 1])
    rmse = float(np.sqrt(np.mean((y_true[ok] - y_pred[ok]) ** 2)))
    return r, rmse


def fit_predict(model, xtr, ytr, xte):
    if model == "ridge":
        m = RidgeCV(alphas=ALPHAS)
    elif model == "gbm":
        m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, max_depth=3,
                                          random_state=config.SEED, early_stopping=True)
    elif model == "mean":
        return np.full(len(xte), float(np.mean(ytr)))
    m.fit(xtr, ytr)
    return m.predict(xte)


def main():
    g = np.load(DERIVED / "genotypes.npy")
    struct = pd.read_csv(DERIVED / "structure.csv")
    phen = pd.read_csv(DERIVED / "phenotypes.csv")
    traits = [c for c in phen.columns if c not in ("fid", "iid")]
    groups = struct.group.to_numpy()
    # markers with the strongest variance make the small-marker baseline; chosen inside each training fold below
    rows = []

    for trait in traits:
        y = phen[trait].to_numpy(dtype=float)
        ok = np.isfinite(y)
        if ok.sum() < 150:
            continue
        gi, yi, grp = g[ok], y[ok], groups[ok]

        for rep in range(REPEATS):
            kf = KFold(FOLDS, shuffle=True, random_state=config.SEED + rep)
            for fold, (tr, te) in enumerate(kf.split(gi)):
                mu, sd = gi[tr].mean(axis=0), gi[tr].std(axis=0)
                sd[sd == 0] = 1
                xtr, xte = (gi[tr] - mu) / sd, (gi[te] - mu) / sd
                for model in ("mean", "ridge", "gbm"):
                    r, rmse = accuracy(yi[te], fit_predict(model, xtr, yi[tr], xte))
                    rows.append({"trait": trait, "test": "random 5-fold", "model": model, "rep": rep,
                                 "fold": fold, "n_train": len(tr), "n_test": len(te), "r": r, "rmse": rmse})
            # a small-marker baseline, using the most variable markers chosen on the training data only
            kf = KFold(FOLDS, shuffle=True, random_state=config.SEED + rep)
            for fold, (tr, te) in enumerate(kf.split(gi)):
                keep = np.argsort(-gi[tr].std(axis=0))[:100]
                mu, sd = gi[tr][:, keep].mean(axis=0), gi[tr][:, keep].std(axis=0)
                sd[sd == 0] = 1
                r, rmse = accuracy(yi[te], fit_predict("ridge", (gi[tr][:, keep] - mu) / sd, yi[tr],
                                                       (gi[te][:, keep] - mu) / sd))
                rows.append({"trait": trait, "test": "random 5-fold", "model": "ridge (100 markers)", "rep": rep,
                             "fold": fold, "n_train": len(tr), "n_test": len(te), "r": r, "rmse": rmse})

        for grp_name in sorted(set(grp)):
            te = grp == grp_name
            tr = ~te
            if te.sum() < MIN_GROUP or tr.sum() < 100:
                continue
            mu, sd = gi[tr].mean(axis=0), gi[tr].std(axis=0)
            sd[sd == 0] = 1
            for model in ("mean", "ridge", "gbm"):
                r, rmse = accuracy(yi[te], fit_predict(model, (gi[tr] - mu) / sd, yi[tr], (gi[te] - mu) / sd))
                rows.append({"trait": trait, "test": "leave-one-group-out", "model": model, "group": grp_name,
                             "n_train": int(tr.sum()), "n_test": int(te.sum()), "r": r, "rmse": rmse})
        log(f"  {trait}: done ({ok.sum()} varieties)")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "prediction_folds.csv", index=False)
    summary = (df.groupby(["trait", "test", "model"])
               .agg(r=("r", "mean"), r_sd=("r", "std"), rmse=("rmse", "mean"), n=("r", "size"))
               .reset_index())
    summary.to_csv(RESULTS / "prediction_summary.csv", index=False)
    piv = summary[summary.test == "random 5-fold"].pivot(index="trait", columns="model", values="r")
    log("\nmean accuracy (r), random 5-fold:\n" + piv.round(3).to_string())
    across = summary[summary.test == "leave-one-group-out"].pivot(index="trait", columns="model", values="r")
    if len(across):
        log("\nmean accuracy (r), leave-one-group-out:\n" + across.round(3).to_string())
    save_json(RESULTS / "prediction.json", {
        "repeats": REPEATS, "folds": FOLDS,
        "summary": summary.replace({np.nan: None}).to_dict("records"),
        "mean_r_random": {m: float(piv[m].mean()) for m in piv.columns},
        "mean_r_across_group": {m: float(across[m].mean()) for m in across.columns} if len(across) else {},
    })


if __name__ == "__main__":
    main()
