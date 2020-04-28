

# Introduction

*[SingleR](https://bioconductor.org/packages/3.11/SingleR)* is an automatic annotation method for single-cell RNA sequencing (scRNAseq) data [@aran2019reference].
Given a reference dataset of samples (single-cell or bulk) with known labels, it labels new cells from a test dataset based on similarity to the reference set.
Specifically, for each test cell:

1. We compute the Spearman correlation between its expression profile and that of each reference sample.
   This is done across the union of marker genes identified between all pairs of labels.
2. We define the per-label score as a fixed quantile (by default, 0.8) of the distribution of correlations.
3. We repeat this for all labels and we take the label with the highest score as the annotation for this cell.
4. We optionally perform a fine-tuning step:
  - The reference dataset is subsetted to only include labels with scores close to the maximum.
  - Scores are recomputed using only marker genes for the subset of labels.
  - This is iterated until one label remains.

Automatic annotation provides a convenient way of transferring biological knowledge across datasets.
In this manner, the burden of manually interpreting clusters and defining marker genes only has to be done once, for the reference dataset, and this knowledge can be propagated to new datasets in an automated manner.

