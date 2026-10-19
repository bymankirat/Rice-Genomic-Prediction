"""Step 1: download the Rice Diversity Panel 1 genotypes (PLINK) and phenotypes, with fingerprints.

The host serves a plain directory index, and the file names on it have changed over the years, so the index is
parsed and the right files picked from it rather than hard-coding a name that may 404. Candidate names are kept as
a fallback in case the index is ever unavailable.
"""
import re
import zipfile

from util import RAW, RESULTS, config, download_any, get, log, md5, save_json


def index_links():
    """-> absolute URLs listed on the data page, in the order they appear."""
    for page in config.DATA_PAGES:
        r = get(page, tries=2)
        if r is None:
            continue
        (RESULTS / "source_index.html").write_text(r.text)
        base = page if page.endswith("/") else page + "/"
        hrefs = re.findall(r'href=["\']([^"\']+)["\']', r.text, re.I)
        out = []
        for h in hrefs:
            if h.startswith("?") or h.startswith("#") or h.lower().startswith("javascript"):
                continue
            out.append(h if h.startswith("http") else base + h.lstrip("./"))
        log(f"index {page}: {len(out)} links")
        return out, page
    return [], None


def pick(links, must, exts):
    """Links whose name contains every word in `must` and ends with one of `exts`, longest match first."""
    hits = [u for u in links
            if u.lower().endswith(tuple(exts)) and all(w.lower() in u.lower() for w in must)]
    return sorted(set(hits), key=len)


def main():
    links, page = index_links()
    geno_urls = pick(links, ["genotype"], [".zip"]) + config.GENOTYPE_URLS
    phen_urls = pick(links, ["phenotype"], [".txt", ".csv"]) + config.PHENOTYPE_URLS
    log("genotype candidates: " + "; ".join(geno_urls[:5]))
    log("phenotype candidates: " + "; ".join(phen_urls[:5]))

    g, gsrc = download_any(geno_urls, RAW / "genotypes.zip")
    p, psrc = download_any(phen_urls, RAW / "phenotypes.txt")
    with zipfile.ZipFile(g) as z:
        names = z.namelist()
        z.extractall(RAW / "genotypes")
    log("genotype archive contains: " + ", ".join(names))

    with open(p, errors="ignore") as fh:
        header = fh.readline().strip()
    prov = {"index": page,
            "genotypes": {"url": gsrc, "md5": md5(g), "bytes": g.stat().st_size, "files": names},
            "phenotypes": {"url": psrc, "md5": md5(p), "bytes": p.stat().st_size},
            "phenotype_header": header.split("\t") if "\t" in header else header.split()}
    log(f"phenotype columns ({len(prov['phenotype_header'])}): {prov['phenotype_header'][:6]} …")
    save_json(RESULTS / "provenance.json", prov)


if __name__ == "__main__":
    main()
