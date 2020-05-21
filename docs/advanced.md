---
bibliography: ref.bib
---

# Advanced options

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("aaron-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.aaron-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.aaron-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>

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
(Almost; see comments below about `common`.)
This yields a list of various components that contains all identified marker genes
and precomputed rank indices to be used in the score calculation.
We can also turn on aggregation with `aggr.ref=TRUE` as described in Chapter \@ref(using-single-cell-references).


```r
common <- intersect(rownames(sce), rownames(dice))

set.seed(2000)
trained <- trainSingleR(dice[common,], labels=dice$label.fine, aggr.ref=TRUE)
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
##                        346                        514 
##           Monocytes, CD16+                   NK cells 
##                        187                        325 
## T cells, CD4+, memory TREG       T cells, CD4+, naive 
##                        163                         97 
##  T cells, CD4+, naive TREG         T cells, CD4+, TFH 
##                         41                        456 
##         T cells, CD4+, Th1      T cells, CD4+, Th1_17 
##                        180                         73 
##        T cells, CD4+, Th17         T cells, CD4+, Th2 
##                         60                         46 
##       T cells, CD8+, naive 
##                        212
```

```r
# Comparing to the direct approach.
set.seed(2000)
direct <- SingleR(sce, ref=dice, labels=dice$label.fine,
    assay.type.test=1, aggr.ref=TRUE)
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
```

Alternatively, one can use separate processes with `SnowParam()`, 
which is slower but can be used on all systems - including Windows, our old nemesis.


```r
pred2b <- SingleR(sce, ref=dice, assay.type.test=1, labels=dice$label.fine,
    BPPARAM=SnowParam(8))
identical(pred$labels, pred2b$labels) 
```

```
## [1] FALSE
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
pred3a <- SingleR(sce, ref=dice, assay.type.test=1, 
    labels=dice$label.main, fine.tune=FALSE)
table(pred3a$labels)
```

```
## 
##       B cells     Monocytes      NK cells T cells, CD4+ T cells, CD8+ 
##           348           705           357           950           340
```

Another approximation is based on the fact that the initial score calculation is done using a nearest-neighbors search.
By default, this is an exact seach but we can switch to an approximate algorithm via the `BNPARAM=` argument.
In the example below, we use the [Annoy algorithm](https://github.com/spotify/annoy) 
via the *[BiocNeighbors](https://bioconductor.org/packages/3.12/BiocNeighbors)* framework, which yields mostly similar results.
(Note, though, that the Annoy method does involve a considerable amount of overhead,
so for small jobs it will actually be slower than the exact search.)


```r
library(BiocNeighbors)
pred3b <- SingleR(sce, ref=dice, assay.type.test=1, 
    labels=dice$label.main, fine.tune=FALSE, # for comparison with pred3a.
    BNPARAM=AnnoyParam())
table(pred3a$labels, pred3b$labels)
```

```
##                
##                 B cells Monocytes NK cells T cells, CD4+ T cells, CD8+
##   B cells           348         0        0             0             0
##   Monocytes           0       705        0             0             0
##   NK cells            0         0      357             0             0
##   T cells, CD4+       0         0        0           950             0
##   T cells, CD8+       0         0        0             0           340
```



## Session information {-}

<button class="aaron-collapse">View session info</button>
<div class="aaron-content">
```
R version 4.0.0 Patched (2020-05-01 r78341)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS:   /home/luna/Software/R/R-4-0-branch-dev/lib/libRblas.so
LAPACK: /home/luna/Software/R/R-4-0-branch-dev/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] BiocNeighbors_1.7.0         BiocParallel_1.23.0        
 [3] TENxPBMCData_1.7.0          HDF5Array_1.17.0           
 [5] rhdf5_2.33.0                SingleCellExperiment_1.11.1
 [7] ensembldb_2.13.1            AnnotationFilter_1.13.0    
 [9] GenomicFeatures_1.41.0      AnnotationDbi_1.51.0       
[11] SingleR_1.3.4               SummarizedExperiment_1.19.4
[13] DelayedArray_0.15.1         matrixStats_0.56.0         
[15] Biobase_2.49.0              GenomicRanges_1.41.1       
[17] GenomeInfoDb_1.25.0         IRanges_2.23.4             
[19] S4Vectors_0.27.6            BiocGenerics_0.35.2        
[21] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] ProtGenerics_1.21.0           bitops_1.0-6                 
 [3] bit64_0.9-7                   progress_1.2.2               
 [5] httr_1.4.1                    tools_4.0.0                  
 [7] R6_2.4.1                      irlba_2.3.3                  
 [9] lazyeval_0.2.2                DBI_1.1.0                    
[11] prettyunits_1.1.1             tidyselect_1.1.0             
[13] processx_3.4.2                bit_1.1-15.2                 
[15] curl_4.3                      compiler_4.0.0               
[17] graph_1.67.0                  rtracklayer_1.49.1           
[19] bookdown_0.19                 askpass_1.1                  
[21] callr_3.4.3                   rappdirs_0.3.1               
[23] Rsamtools_2.5.0               stringr_1.4.0                
[25] digest_0.6.25                 rmarkdown_2.1                
[27] XVector_0.29.0                pkgconfig_2.0.3              
[29] htmltools_0.4.0               dbplyr_1.4.3                 
[31] fastmap_1.0.1                 rlang_0.4.6                  
[33] RSQLite_2.2.0                 shiny_1.4.0.2                
[35] DelayedMatrixStats_1.11.0     dplyr_0.8.5                  
[37] RCurl_1.98-1.2                magrittr_1.5                 
[39] BiocSingular_1.5.0            GenomeInfoDbData_1.2.3       
[41] Matrix_1.2-18                 Rhdf5lib_1.11.0              
[43] Rcpp_1.0.4.6                  lifecycle_0.2.0              
[45] stringi_1.4.6                 yaml_2.2.1                   
[47] zlibbioc_1.35.0               BiocFileCache_1.13.0         
[49] AnnotationHub_2.21.0          grid_4.0.0                   
[51] blob_1.2.1                    promises_1.1.0               
[53] ExperimentHub_1.15.0          crayon_1.3.4                 
[55] lattice_0.20-41               Biostrings_2.57.0            
[57] hms_0.5.3                     CodeDepends_0.6.5            
[59] knitr_1.28                    ps_1.3.3                     
[61] pillar_1.4.4                  codetools_0.2-16             
[63] biomaRt_2.45.0                XML_3.99-0.3                 
[65] glue_1.4.1                    BiocVersion_3.12.0           
[67] evaluate_0.14                 BiocManager_1.30.10          
[69] vctrs_0.3.0                   httpuv_1.5.2                 
[71] openssl_1.4.1                 purrr_0.3.4                  
[73] assertthat_0.2.1              xfun_0.13                    
[75] rsvd_1.0.3                    mime_0.9                     
[77] xtable_1.8-4                  later_1.0.0                  
[79] tibble_3.0.1                  GenomicAlignments_1.25.0     
[81] memoise_1.1.0                 ellipsis_0.3.1               
[83] interactiveDisplayBase_1.27.0
```
</div>
