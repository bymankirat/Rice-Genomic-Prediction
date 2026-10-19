"""Step 1: download the Rice Diversity Panel 1 genotypes (PLINK) and phenotypes, with fingerprints."""
import zipfile

from util import RAW, RESULTS, config, download_any, get, log, md5, save_json


def main():
    prov = {}
    for page in config.DATA_PAGES:  # keep a copy of the index page so the file list is on record
        r = get(page, tries=2)
        if r is not None:
            (RESULTS / "source_index.html").write_text(r.text)
            prov["index"] = page
            break
    g, gsrc = download_any(config.GENOTYPE_URLS, RAW / "genotypes.zip")
    p, psrc = download_any(config.PHENOTYPE_URLS, RAW / "phenotypes.txt")
    with zipfile.ZipFile(g) as z:
        names = z.namelist()
        z.extractall(RAW / "genotypes")
    log("genotype archive:", names)
    prov.update({"genotypes": {"url": gsrc, "md5": md5(g), "files": names},
                 "phenotypes": {"url": psrc, "md5": md5(p)}})
    with open(p) as fh:
        head = [next(fh) for _ in range(3)]
    prov["phenotype_header"] = head[0].strip().split("\t")
    log("phenotype columns:", prov["phenotype_header"])
    save_json(RESULTS / "provenance.json", prov)


if __name__ == "__main__":
    main()
