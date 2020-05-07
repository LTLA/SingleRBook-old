# Advanced options



## Improving efficiency

Advanced users can split the `SingleR()` workflow into two separate training and classification steps.
This means that training (e.g., marker detection, assembling of nearest-neighbor indices) only needs to be performed once.
The resulting data structures can then be re-used across multiple classifications with different test datasets, provided the test feature set is identical to or a superset of the features in the training set.
For example:


```r
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')

library(scater)
hESCs <- logNormCounts(hESCs)

library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
common <- intersect(rownames(hESCs), rownames(hpca.se))
trained <- trainSingleR(hpca.se[common,], labels=hpca.se$label.main)
pred.hesc2 <- classifySingleR(hESCs[common,], trained)
table(pred.hesc$labels, pred.hesc2$labels)
```

```
##                       
##                        Astrocyte Chondrocytes Embryonic_stem_cells iPS_cells
##   Astrocyte                   32            0                    0         0
##   Chondrocytes                 0            1                    0         0
##   Embryonic_stem_cells         0            0                  127         0
##   iPS_cells                    0            0                    0       195
##   Neuroepithelial_cell         0            0                    0         0
##   Neurons                      0            0                    0         0
##   Smooth_muscle_cells          0            0                    0         0
##                       
##                        Neuroepithelial_cell Neurons Smooth_muscle_cells
##   Astrocyte                               0       0                   0
##   Chondrocytes                            0       0                   0
##   Embryonic_stem_cells                    0       0                   0
##   iPS_cells                               0       0                   0
##   Neuroepithelial_cell                 1030       0                   0
##   Neurons                                 0     325                   0
##   Smooth_muscle_cells                     0       0                   5
```

Other efficiency improvements are possible through several arguments:

- Switching to an approximate algorithm for the nearest neighbor search in `trainSingleR()` via the `BNPARAM=` argument from the *[BiocNeighbors](https://bioconductor.org/packages/3.12/BiocNeighbors)* package.
- Parallelizing the fine-tuning step in `classifySingleR()` with the `BPPARAM=` argument from the *[BiocParallel](https://bioconductor.org/packages/3.12/BiocParallel)* package.

These arguments can also be specified in the `SingleR()` command.
