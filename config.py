"""Settings for the rice genomic prediction pipeline (Rice Diversity Panel 1, 44K SNP array)."""

# Zhao et al. 2011, Nature Communications 2:467 - 413 diverse rice accessions genotyped at ~44,000 SNPs,
# phenotyped for 34 traits. Several mirrors are tried in turn.
DATA_PAGES = ["http://www.ricediversity.org/data/sets/44kgwas/", "https://www.ricediversity.org/data/sets/44kgwas/"]
GENOTYPE_URLS = [
    "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Genotypes_PLINK.zip",
    "https://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Genotypes_PLINK.zip",
]
PHENOTYPE_URLS = [
    "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
    "https://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
]
SEED = 20261019
