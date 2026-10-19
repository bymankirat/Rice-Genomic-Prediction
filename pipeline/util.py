"""Shared helpers: paths, logging, HTTP with retries and fallbacks."""
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import config  # noqa: E402

RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
RESULTS = ROOT / "results"
for _d in (RAW, DERIVED, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)


def log(*parts):
    print(time.strftime("%H:%M:%S"), *parts, flush=True)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    def conv(x):
        import numpy as np
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.floating):
            return None if np.isnan(x) else float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x))
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=conv)


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def get(url, tries=5, **kw):
    import requests
    kw.setdefault("timeout", 300)
    for i in range(tries):
        try:
            r = requests.get(url, headers={"User-Agent": "rice-genomic-prediction (github.com/bymankirat)"}, **kw)
            if r.status_code < 400:
                return r
            log(f"  HTTP {r.status_code}: {url}")
            if r.status_code == 404:
                return None
        except requests.RequestException as e:
            log(f"  {e.__class__.__name__}: {url}")
        time.sleep(2 ** i)
    return None


def download_any(urls, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest, "cache"
    for url in urls:
        r = get(url, stream=True)
        if r is None:
            continue
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
        tmp.rename(dest)
        return dest, url
    raise RuntimeError(f"all sources failed for {dest.name}")
