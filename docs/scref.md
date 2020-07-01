---
bibliography: ref.bib
---

# Controlling marker detection {#more-markers}

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

## Overview

One of the most important steps in *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* (beyond the choice of reference, of course)
is the derivation of the marker genes used in the score calculation.
We have already introduced the classic approach in the previous chapter,
but it is similarly straightforward to perform marker detection with conventional statistical tests.
In particular, we identify top-ranked markers based on pairwise Wilcoxon rank sum tests or $t$-tests between labels;
this allows us to account for the variability across cells to choose genes that are robustly upregulated in each label.

The availability of variance-aware marker detection methods is most relevant for reference datasets
that contain a reasonable number (i.e., at least two) of replicate samples for each label.
An obvious use case is that of single-cell datasets that are used as a reference to annotate other single-cell datasets.
It is also possible for users to supply their own custom marker lists to `SingleR()`,
facilitating incorporation of prior biological knowledge into the annotation process.
We will demonstrate these capabilities below in this chapter.

## Annotation with test-based marker detection

To demonstrate, we will use two human pancreas scRNA-seq datasets from the *[scRNAseq](https://bioconductor.org/packages/3.12/scRNAseq)* package.
The aim is to use one pre-labelled dataset to annotate the other unlabelled dataset.
First, we set up the @muraro2016singlecell dataset to be our reference,
computing log-normalized expression values as discussed in Section \@ref(choices-of-assay-data).


```r
library(scRNAseq)
sceM <- MuraroPancreasData()

# Removing unlabelled cells or cells without a clear label.
sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 

library(scater)
sceM <- logNormCounts(sceM)
sceM
```

```
## class: SingleCellExperiment 
## dim: 19059 2122 
## metadata(0):
## assays(2): counts logcounts
## rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(2): symbol chr
## colnames(2122): D28-1_1 D28-1_2 ... D30-8_93 D30-8_94
## colData names(4): label donor plate sizeFactor
## reducedDimNames(0):
## altExpNames(1): ERCC
```

```r
# Seeing the available labels in this dataset.
table(sceM$label)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         219         812         448         193         245          21 
##     epsilon mesenchymal          pp 
##           3          80         101
```

We then set up our test dataset from @grun2016denovo, applying some basic quality control as discusssed
[here](https://osca.bioconductor.org/grun-human-pancreas-cel-seq2.html#quality-control-8) 
and in Section \@ref(interaction-with-quality-control).
We also compute the log-transformed values here, not because it is strictly necessary
but so that we don't have to keep on typing `assay.type.test=1` in later calls to `SingleR()`.


```r
sceG <- GrunPancreasData()

sceG <- addPerCellQC(sceG)
qc <- quickPerCellQC(colData(sceG), 
    percent_subsets="altexps_ERCC_percent",
    batch=sceG$donor,
    subset=sceG$donor %in% c("D17", "D7", "D2"))
sceG <- sceG[,!qc$discard]

sceG <- logNormCounts(sceG)
sceG
```

```
## class: SingleCellExperiment 
## dim: 20064 1064 
## metadata(0):
## assays(2): counts logcounts
## rownames(20064): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(2): symbol chr
## colnames(1064): D2ex_1 D2ex_2 ... D17TGFB_94 D17TGFB_95
## colData names(9): donor sample ... total sizeFactor
## reducedDimNames(0):
## altExpNames(1): ERCC
```

We run `SingleR()` as described previously but with a marker detection mode that considers the variance of expression across cells.
Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels.
This is slower but more appropriate for single-cell data compared to the default marker detection algorithm,
as the latter may fail for low-coverage data where the median for each label is often zero.


```r
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
table(pred.grun$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         277         203         181          50         306           5 
##     epsilon mesenchymal          pp 
##           1          22          19
```



By default, the function will take the top `de.n` (default: 10) genes from each pairwise comparison between labels.
A larger number of markers increases the robustness of the annotation by ensuring that relevant genes are not omitted,
especially if the reference dataset has study-specific effects that cause uninteresting genes to dominate the top set. 
However, this comes at the cost of increasing noise and computational time.


```r
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="wilcox", de.n=50)
table(pred.grun$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         275         203         177          55         307           5 
##     epsilon mesenchymal          pp 
##           1          23          18
```

## Defining custom markers

The marker detection in `SingleR()` is built on top of the testing framework in *[scran](https://bioconductor.org/packages/3.12/scran)*,
so most options in `?pairwiseWilcox` and friends can be applied via the `de.args=` option.
For example, we could use the $t$-test and test against a log-fold change threshold with `de.args=list(lfc=1)`.


```r
library(SingleR)
pred.grun2 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="t", de.args=list(lfc=1))
table(pred.grun2$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         285         200         177          54         296           5 
##     epsilon mesenchymal          pp 
##           5          24          18
```

However, users can also construct their own marker lists with any DE testing machinery.
For example, we can perform pairwise binomial tests to identify genes that are differentially detected 
(i.e., have differences in the proportion of cells with non-zero counts) between labels in the reference Muraro dataset.
We then take the top 10 marker genes from each pairwise comparison,
obtaining a list of lists of character vectors containing the identities of the markers for that comparison.


```r
library(scran)
out <- pairwiseBinom(counts(sceM), sceM$label, direction="up")
markers <- getTopMarkers(out$statistics, out$pairs, n=10)

# Upregulated in acinar compared to alpha:
markers$acinar$alpha
```

```
##  [1] "KCNQ1__chr11"  "FAM129A__chr1" "KLK1__chr19"   "NTN4__chr12"  
##  [5] "RASEF__chr9"   "CTRL__chr16"   "LGALS2__chr22" "NUPR1__chr16" 
##  [9] "LGALS3__chr14" "NR5A2__chr1"
```

```r
# Upregulated in alpha compared to acinar:
markers$alpha$acinar
```

```
##  [1] "SLC38A4__chr12" "ARX__chrX"      "CRYBA2__chr2"   "FSTL5__chr4"   
##  [5] "GNG2__chr14"    "NOL4__chr18"    "IRX2__chr5"     "KCNMB2__chr3"  
##  [9] "CFC1__chr2"     "KCNJ6__chr21"
```

Once we have this list of lists, we supply it to `SingleR()` via the `genes=` argument,
which causes the function to bypass the internal marker detection to use the supplied gene sets instead.
The most obvious benefit of this approach is that the user can achieve greater control of the markers,
allowing integration of prior biological knowledge to obtain more relevant genes and a more robust annotation.


```r
pred.grun2b <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=markers)
table(pred.grun2b$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         276         202         175          54         302           5 
##     epsilon mesenchymal          pp 
##           2          25          23
```

In some cases, markers may only be available for specific labels rather than for pairwise comparisons between labels.
This is accommodated by supplying a named list of character vectors to `genes`.
Note that this is likely to be less powerful than the list-of-lists approach as information about pairwise differences is discarded.


```r
# Creating label-specific markers.
label.markers <- lapply(markers, unlist)
label.markers <- lapply(label.markers, unique)
str(label.markers)
```

```
## List of 9
##  $ acinar     : chr [1:40] "KCNQ1__chr11" "FAM129A__chr1" "KLK1__chr19" "NTN4__chr12" ...
##  $ alpha      : chr [1:41] "SLC38A4__chr12" "ARX__chrX" "CRYBA2__chr2" "FSTL5__chr4" ...
##  $ beta       : chr [1:47] "ELAVL4__chr1" "PRUNE2__chr9" "NMNAT2__chr1" "PLCB4__chr20" ...
##  $ delta      : chr [1:44] "NOL4__chr18" "CABP7__chr22" "UNC80__chr2" "HEPACAM2__chr7" ...
##  $ duct       : chr [1:50] "ADCY5__chr3" "PDE3A__chr12" "SLC3A1__chr2" "BICC1__chr10" ...
##  $ endothelial: chr [1:26] "GPR4__chr19" "TMEM204__chr16" "GPR116__chr6" "CYYR1__chr21" ...
##  $ epsilon    : chr [1:14] "BHMT__chr5" "JPH3__chr16" "SERPINA10__chr14" "UGT2B4__chr4" ...
##  $ mesenchymal: chr [1:34] "TNFAIP6__chr2" "THBS2__chr6" "CDH11__chr16" "SRPX2__chrX" ...
##  $ pp         : chr [1:44] "SERTM1__chr13" "ETV1__chr7" "ARX__chrX" "ELAVL4__chr1" ...
```

```r
pred.grun2c <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=label.markers)
table(pred.grun2c$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         262         204         169          59         317           6 
##     epsilon mesenchymal          pp 
##           2          24          21
```

## Pseudo-bulk aggregation

Single-cell reference datasets provide a like-for-like comparison to our test single-cell datasets, 
yielding a more accurate classification of the cells in the latter (hopefully).
However, there are frequently many more samples in single-cell references compared to bulk references, 
increasing the computational work involved in classification.
We overcome this by aggregating cells into one "pseudo-bulk" sample per label (e.g., by averaging across log-expression values) 
and using that as the reference profile, which allows us to achieve the same efficiency as the use of bulk references.

The obvious cost of this approach is that we discard potentially useful information 
about the distribution of cells within each label.
Cells that belong to a heterogeneous population may not be correctly assigned if they are far from the population center.
To preserve some of this information, we perform $k$-means clustering within each label
to create pseudo-bulk samples that are representative of a particular region of the expression space (i.e., vector quantization).
We create $\sqrt{N}$ clusters given a label with $N$ cells, 
which provides a reasonable compromise between reducing computational work and preserving the label's internal distribution.

To enable this aggregation, we simply set `aggr.ref=TRUE` in the `SingleR()` call.
This uses the `aggregateReference()` function to perform $k$-means clustering _within_ each label
(typically after principal components analysis on the log-expression matrix, for greater speed)
and average expression values for each within-label cluster.
Note that marker detection is still performed on the unaggregated data 
so as to make full use of the distribution of expression values across cells.


```r
set.seed(100) # for the k-means step.
pred.grun3 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="wilcox", aggr.ref=TRUE)
table(pred.grun3$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         277         202         184          47         306           5 
##     epsilon mesenchymal          pp 
##           1          22          20
```

Obviously, the aggregation itself requires computational work so setting `aggr.ref=TRUE` in `SingleR()` itself may not improve speed.
Rather, the real power of this approach lies in pre-aggregating the reference dataset
so that it can be repeatedly applied to quickly annotate multiple test datasets.
This approach is discussed in more detail in Chapter \@ref(advanced-options).

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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] scran_1.17.4                SingleR_1.3.6              
 [3] scater_1.17.3               ggplot2_3.3.2              
 [5] scRNAseq_2.3.6              SingleCellExperiment_1.11.6
 [7] SummarizedExperiment_1.19.5 DelayedArray_0.15.6        
 [9] matrixStats_0.56.0          Matrix_1.2-18              
[11] Biobase_2.49.0              GenomicRanges_1.41.5       
[13] GenomeInfoDb_1.25.5         IRanges_2.23.10            
[15] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[17] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] bitops_1.0-6                  bit64_0.9-7                  
 [3] httr_1.4.1                    tools_4.0.0                  
 [5] R6_2.4.1                      irlba_2.3.3                  
 [7] vipor_0.4.5                   DBI_1.1.0                    
 [9] colorspace_1.4-1              withr_2.2.0                  
[11] gridExtra_2.3                 tidyselect_1.1.0             
[13] processx_3.4.2                bit_1.1-15.2                 
[15] curl_4.3                      compiler_4.0.0               
[17] graph_1.67.1                  BiocNeighbors_1.7.0          
[19] bookdown_0.20                 scales_1.1.1                 
[21] callr_3.4.3                   rappdirs_0.3.1               
[23] stringr_1.4.0                 digest_0.6.25                
[25] rmarkdown_2.3                 XVector_0.29.3               
[27] pkgconfig_2.0.3               htmltools_0.5.0              
[29] limma_3.45.7                  dbplyr_1.4.4                 
[31] fastmap_1.0.1                 rlang_0.4.6                  
[33] RSQLite_2.2.0                 shiny_1.5.0                  
[35] DelayedMatrixStats_1.11.1     generics_0.0.2               
[37] BiocParallel_1.23.0           dplyr_1.0.0                  
[39] RCurl_1.98-1.2                magrittr_1.5                 
[41] BiocSingular_1.5.0            GenomeInfoDbData_1.2.3       
[43] scuttle_0.99.10               Rcpp_1.0.4.6                 
[45] ggbeeswarm_0.6.0              munsell_0.5.0                
[47] viridis_0.5.1                 lifecycle_0.2.0              
[49] edgeR_3.31.4                  stringi_1.4.6                
[51] yaml_2.2.1                    zlibbioc_1.35.0              
[53] BiocFileCache_1.13.0          AnnotationHub_2.21.1         
[55] grid_4.0.0                    blob_1.2.1                   
[57] dqrng_0.2.1                   promises_1.1.1               
[59] ExperimentHub_1.15.0          crayon_1.3.4                 
[61] lattice_0.20-41               locfit_1.5-9.4               
[63] CodeDepends_0.6.5             knitr_1.29                   
[65] ps_1.3.3                      pillar_1.4.4                 
[67] igraph_1.2.5                  codetools_0.2-16             
[69] XML_3.99-0.3                  glue_1.4.1                   
[71] BiocVersion_3.12.0            evaluate_0.14                
[73] BiocManager_1.30.10           vctrs_0.3.1                  
[75] httpuv_1.5.4                  gtable_0.3.0                 
[77] purrr_0.3.4                   assertthat_0.2.1             
[79] xfun_0.15                     rsvd_1.0.3                   
[81] mime_0.9                      xtable_1.8-4                 
[83] later_1.1.0.1                 viridisLite_0.3.0            
[85] tibble_3.0.1                  AnnotationDbi_1.51.1         
[87] beeswarm_0.2.3                memoise_1.1.0                
[89] statmod_1.4.34                ellipsis_0.3.1               
[91] interactiveDisplayBase_1.27.5
```
</div>
