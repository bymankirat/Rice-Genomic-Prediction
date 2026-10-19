"""Step 6: write docs/results.js for the report page."""
import json

import numpy as np
import pandas as pd

from util import DERIVED, RESULTS, ROOT, load_json, log

DOCS = ROOT / "docs"


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if o is None or (isinstance(o, float) and np.isnan(o)) else float(o)
    if o is pd.NA:
        return None
    return o


def csv(path, **kw):
    return clean(pd.read_csv(path, **kw).to_dict("records")) if path.exists() else []


def js(path):
    return clean(load_json(path)) if path.exists() else {}


def main():
    DOCS.mkdir(exist_ok=True)
    struct = pd.read_csv(DERIVED / "structure.csv")
    pcs = struct[["iid", "group", "PC1", "PC2", "PC3"]].round(3)
    data = {
        "generated": pd.Timestamp.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
        "provenance": js(RESULTS / "provenance.json"),
        "qc": js(RESULTS / "qc.json"),
        "structure": js(RESULTS / "structure.json"),
        "pcs": clean(pcs.to_dict("records")),
        "gwas": js(RESULTS / "gwas.json"),
        "prediction": js(RESULTS / "prediction.json"),
        "prediction_summary": csv(RESULTS / "prediction_summary.csv"),
        "markers_per_chrom": clean(pd.read_csv(DERIVED / "markers.csv").groupby("chrom").size()
                                   .rename("markers").reset_index().to_dict("records")),
    }
    out = "window.R = " + json.dumps(data, separators=(",", ":")) + ";\n"
    (DOCS / "results.js").write_text(out)
    log(f"wrote docs/results.js ({len(out)/1e3:.0f} kB)")


if __name__ == "__main__":
    main()
