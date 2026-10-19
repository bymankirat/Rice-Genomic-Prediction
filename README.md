# Can you predict a rice variety's traits from its DNA?

Plant breeding is slow because a cross has to be grown out before anyone knows whether it was worth making.
**Genomic prediction** tries to shorten that: train a model on varieties whose DNA and traits are both known,
then predict the traits of new varieties from DNA alone. This project tests how well that works, and where it
breaks, on a public rice panel.

**Data:** the Rice Diversity Panel 1 (Zhao et al. 2011, *Nature Communications* 2:467): **413 rice varieties from
82 countries**, genotyped at ~44,000 SNPs and measured for 34 traits.

## Questions

| Question | How it is answered |
|---|---|
| How is the panel structured? | PCA on the marker matrix; genetic groups from k-means, chosen by silhouette |
| Does the association analysis work at all? | GWAS per trait with structure correction, then checked against genes that are already textbook: *sd1* (plant height), *Waxy* (amylose), *GS3* (grain length) |
| How accurately can traits be predicted? | Ridge regression on all markers and gradient boosting, 5×5-fold cross-validation |
| Does a model transfer between rice types? | Leave-one-group-out: train on every genetic group but one, predict the group never seen |
| Is the gain real? | Two baselines in the same folds: predicting the training mean, and ridge on 100 markers |

## Pipeline

```bash
pip install -r requirements.txt
python run_all.py        # downloads the panel, then QC → structure → GWAS → prediction → report
python -m pytest -q      # method tests on simulated data, plus checks on the saved results
```

The GitHub workflow runs the whole thing and commits the results back here, so every number on the report page
comes from a file the pipeline wrote.

## Honesty notes

- **Structure is the confounder.** Rice groups differ in both markers and traits, so an uncorrected GWAS produces
  associations all over the genome. The tests include a simulation where exactly that happens, and show the
  covariates remove it.
- **Random cross-validation flatters a model.** When training and test varieties come from the same groups, part
  of the accuracy is relatedness, not biology. The leave-one-group-out test is the one that matters for a breeder
  working with a new population.
- **Markers are not causes.** A marker associated with a trait is usually just close to whatever causes it.

Mankirat Singh · [Portfolio](https://bymankirat.github.io/myportfolio/) · [GitHub](https://github.com/bymankirat)
