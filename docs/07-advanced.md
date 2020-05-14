---
bibliography: ../ref.bib
---

# Advanced options



## Preconstructed indices

Advanced users can split the `SingleR()` workflow into two separate training and classification steps.
This means that training (e.g., marker detection, assembling of nearest-neighbor indices) only needs to be performed once.
The resulting data structures can then be re-used across multiple classifications with different test datasets, 
provided the gene annotation in the test dataset is identical to or a superset of the genes in the training set.
To illustrate, we will consider the DICE reference dataset [@diceRef].


```r
library(SingleR)
dice <- DatabaseImmuneCellExpressionData(ensembl=TRUE)
dice
```

```
## class: SummarizedExperiment 
## dim: 29914 1561 
## metadata(0):
## assays(1): logcounts
## rownames(29914): ENSG00000121410 ENSG00000268895 ... ENSG00000159840
##   ENSG00000074755
## rowData names(0):
## colnames(1561): TPM_1 TPM_2 ... TPM_101 TPM_102
## colData names(3): label.main label.fine label.ont
```

```r
table(dice$label.fine)
```

```
## 
##                   B cells, naive                 Monocytes, CD14+ 
##                              106                              106 
##                 Monocytes, CD16+                         NK cells 
##                              105                              105 
##       T cells, CD4+, memory TREG             T cells, CD4+, naive 
##                              104                              103 
##        T cells, CD4+, naive TREG T cells, CD4+, naive, stimulated 
##                              104                              102 
##               T cells, CD4+, TFH               T cells, CD4+, Th1 
##                              104                              104 
##            T cells, CD4+, Th1_17              T cells, CD4+, Th17 
##                              104                              104 
##               T cells, CD4+, Th2             T cells, CD8+, naive 
##                              104                              104 
## T cells, CD8+, naive, stimulated 
##                              102
```

Let's say we want to use the DICE reference to annotate the PBMC dataset from Chapter \@ref(introduction).


```r
library(TENxPBMCData)
sce <- TENxPBMCData("pbmc3k")
```



We use the `trainSingleR()` function to do all the necessary calculations 
that are independent of the test dataset.
(Well, almost; see comments below about `common`.)
This yields a list of various components that contains all identified marker genes
and precomputed rank indices to be used in the score calculation.


```r
common <- intersect(rownames(sce), rownames(dice))
trained <- trainSingleR(dice[common,], labels=dice$label.fine)
```

It is then a simple matter to use the `trained` object to annotate our dataset of interest
through the `classifySingleR()` function.
As we can see, this yields exactly the same result as applying `SingleR()` directly;
however, the advantage here is that `trained` can be re-used for multiple `classifySingleR()` calls - 
possibly on different datasets - without having to repeat unnecessary steps when the reference is unchanged.


```r
pred <- classifySingleR(sce, trained, assay.type=1)
table(pred$labels)
```

```
## 
##             B cells, naive           Monocytes, CD14+ 
##                        344                        516 
##           Monocytes, CD16+                   NK cells 
##                        186                        320 
## T cells, CD4+, memory TREG       T cells, CD4+, naive 
##                        148                        111 
##  T cells, CD4+, naive TREG         T cells, CD4+, TFH 
##                         28                        461 
##         T cells, CD4+, Th1      T cells, CD4+, Th1_17 
##                        213                         61 
##        T cells, CD4+, Th17         T cells, CD4+, Th2 
##                         57                         41 
##       T cells, CD8+, naive 
##                        214
```

```r
# Comparing to the direct approach.
direct <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.fine)
identical(pred$labels, direct$labels)
```

```
## [1] TRUE
```



The big caveat is that the universe of genes in the test dataset cannot be greater than that in the reference.
This is the reason behind the intersection to `common` genes and the subsequent subsetting of `dice`.
Practical use of preconstructed indices is best combined with some prior information about the gene-level annotation;
for example, we might know that we always use a particular version of the Ensembl gene models,
so we would filter out any genes in the reference dataset that are not in the Ensembl annotation.

## Parallelization

Parallelization is an obvious approach to increasing annotation throughput.
This is done using the framework in the *[BiocParallel](https://bioconductor.org/packages/3.12/BiocParallel)* package, 
which provides several options for parallelization depending on the available hardware.
On POSIX-compliant systems (i.e., Linux and MacOS), the simplest method is to use forking 
by passing `MulticoreParam()` to the `BPPARAM=` argument:


```r
library(BiocParallel)
pred2a <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.fine,
    BPPARAM=MulticoreParam(8)) # 8 CPUs.
identical(pred$labels, pred2a$labels) 
```

```
## [1] TRUE
```



Alternatively, one can use separate processes with `SnowParam()`, 
which is slower but can be used on all systems including Windows.


```r
pred2b <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.fine,
    BPPARAM=SnowParam(8))
identical(pred$labels, pred2b$labels) 
```

```
## [1] TRUE
```



When working on a cluster, the `BatchtoolsParam()` function allows `SingleR()` to 
seamlessly interface with various job schedulers like SLURM, LSF and so on.
This permits heavy-duty parallelization across hundreds of CPUs for highly intensive jobs,
though often some configuration is required - 
see the [vignette](https://bioconductor.org/packages/3.12/BiocParallel/vignettes/BiocParallel_BatchtoolsParam.pdf) for more details.

## Approximate algorithms

It is possible to sacrifice accuracy to squeeze more speed out of *[SingleR](https://bioconductor.org/packages/3.12/SingleR)*.
The most obvious approach is to simply turn off the fine-tuning with `fine.tune=FALSE`,
which avoids the time-consuming fine-tuning iterations.
When the reference labels are well-separated, this is probably an acceptable trade-off.


```r
pred3a <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.main,
    fine.tune=FALSE)
table(pred3a$labels)
```

```
## 
##       B cells     Monocytes      NK cells T cells, CD4+ T cells, CD8+ 
##           348           705           357           950           340
```

Another approximation is based on the fact that the score calculation is done using a nearest-neighbors search.
By default, this is an exact seach but we can switch to an approximate algorithm via the `BNPARAM=` argument.
In the example below, we use the [Annoy algorithm](https://github.com/spotify/annoy) 
via the *[BiocNeighbors](https://bioconductor.org/packages/3.12/BiocNeighbors)* framework, which yields mostly similar results.
(Note, though, that the Annoy method does involve a considerable amount of overhead,
so for small jobs it will actually be slower than the exact search.)


```r
library(BiocNeighbors)
pred3b <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.main,
    BNPARAM=AnnoyParam())
table(pred3a$labels, pred3b$labels)
```

```
##                
##                 B cells Monocytes NK cells T cells, CD4+ T cells, CD8+
##   B cells           347         1        0             0             0
##   Monocytes           0       705        0             0             0
##   NK cells            0         0      341            11             5
##   T cells, CD4+       0         0        4           900            46
##   T cells, CD8+       0         0        2            59           279
```



## Session information {-}


```r
sessionInfo()
```

```
## R version 4.0.0 Patched (2020-05-01 r78341)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS
## 
## Matrix products: default
## BLAS:   /home/luna/Software/R/R-4-0-branch-dev/lib/libRblas.so
## LAPACK: /home/luna/Software/R/R-4-0-branch-dev/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BiocNeighbors_1.7.0         BiocParallel_1.23.0        
##  [3] TENxPBMCData_1.7.0          HDF5Array_1.17.0           
##  [5] rhdf5_2.33.0                ensembldb_2.13.1           
##  [7] AnnotationFilter_1.13.0     GenomicFeatures_1.41.0     
##  [9] AnnotationDbi_1.51.0        SingleR_1.3.2              
## [11] scater_1.17.0               ggplot2_3.3.0              
## [13] scRNAseq_2.3.0              SingleCellExperiment_1.11.1
## [15] SummarizedExperiment_1.19.2 DelayedArray_0.15.1        
## [17] matrixStats_0.56.0          Biobase_2.49.0             
## [19] GenomicRanges_1.41.1        GenomeInfoDb_1.25.0        
## [21] IRanges_2.23.4              S4Vectors_0.27.5           
## [23] BiocGenerics_0.35.2         BiocStyle_2.17.0           
## 
## loaded via a namespace (and not attached):
##  [1] ggbeeswarm_0.6.0              colorspace_1.4-1             
##  [3] ellipsis_0.3.0                XVector_0.29.0               
##  [5] bit64_0.9-7                   interactiveDisplayBase_1.27.0
##  [7] codetools_0.2-16              knitr_1.28                   
##  [9] Rsamtools_2.5.0               dbplyr_1.4.3                 
## [11] shiny_1.4.0.2                 BiocManager_1.30.10          
## [13] compiler_4.0.0                httr_1.4.1                   
## [15] assertthat_0.2.1              Matrix_1.2-18                
## [17] fastmap_1.0.1                 lazyeval_0.2.2               
## [19] later_1.0.0                   BiocSingular_1.5.0           
## [21] htmltools_0.4.0               prettyunits_1.1.1            
## [23] tools_4.0.0                   rsvd_1.0.3                   
## [25] gtable_0.3.0                  glue_1.4.0                   
## [27] GenomeInfoDbData_1.2.3        dplyr_0.8.5                  
## [29] rappdirs_0.3.1                Rcpp_1.0.4.6                 
## [31] vctrs_0.2.4                   Biostrings_2.57.0            
## [33] ExperimentHub_1.15.0          rtracklayer_1.49.1           
## [35] DelayedMatrixStats_1.11.0     xfun_0.13                    
## [37] stringr_1.4.0                 mime_0.9                     
## [39] lifecycle_0.2.0               irlba_2.3.3                  
## [41] XML_3.99-0.3                  AnnotationHub_2.21.0         
## [43] zlibbioc_1.35.0               scales_1.1.0                 
## [45] hms_0.5.3                     promises_1.1.0               
## [47] ProtGenerics_1.21.0           yaml_2.2.1                   
## [49] curl_4.3                      memoise_1.1.0                
## [51] gridExtra_2.3                 biomaRt_2.45.0               
## [53] stringi_1.4.6                 RSQLite_2.2.0                
## [55] BiocVersion_3.12.0            rlang_0.4.6                  
## [57] pkgconfig_2.0.3               bitops_1.0-6                 
## [59] evaluate_0.14                 lattice_0.20-41              
## [61] purrr_0.3.4                   Rhdf5lib_1.11.0              
## [63] GenomicAlignments_1.25.0      bit_1.1-15.2                 
## [65] tidyselect_1.0.0              magrittr_1.5                 
## [67] bookdown_0.18                 R6_2.4.1                     
## [69] snow_0.4-3                    DBI_1.1.0                    
## [71] pillar_1.4.4                  withr_2.2.0                  
## [73] RCurl_1.98-1.2                tibble_3.0.1                 
## [75] crayon_1.3.4                  BiocFileCache_1.13.0         
## [77] rmarkdown_2.1                 viridis_0.5.1                
## [79] progress_1.2.2                grid_4.0.0                   
## [81] blob_1.2.1                    digest_0.6.25                
## [83] xtable_1.8-4                  httpuv_1.5.2                 
## [85] openssl_1.4.1                 munsell_0.5.0                
## [87] beeswarm_0.2.3                viridisLite_0.3.0            
## [89] vipor_0.4.5                   askpass_1.1
```
