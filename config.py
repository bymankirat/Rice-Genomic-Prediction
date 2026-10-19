"""Settings for the rice genomic prediction pipeline (Rice Diversity Panel 1, 44K SNP array)."""

# Zhao et al. 2011, Nature Communications 2:467 - 413 diverse rice accessions genotyped at ~44,000 SNPs,
# phenotyped for 34 traits. Several mirrors are tried in turn.
DATA_PAGES = ["http://ricediversity.org/data/sets/44kgwas/", "http://www.ricediversity.org/data/sets/44kgwas/"]
# Fallbacks, used only if the index page cannot be read. The site has served several names over the years.
_B = "http://ricediversity.org/data/sets/44kgwas/"
GENOTYPE_URLS = [_B + n for n in (
    "RiceDiversity.44K.MSU6.Genotypes_PLINK.zip",
    "RiceDiversity_44K_Genotypes_PLINK.zip",
    "RiceDiversity.44K.MSU6.Genotypes.PLINK.zip",
)]
PHENOTYPE_URLS = [_B + n for n in (
    "RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
    "RiceDiversity.44K.MSU6.Phenotypes.txt",
)]
SEED = 20261019
