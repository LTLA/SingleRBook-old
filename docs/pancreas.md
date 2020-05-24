---
bibliography: ref.bib
---

# (PART) Workflows {-}

# Cross-annotating pancreas

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

## Loading the data

We load the @muraro2016singlecell dataset as our reference, 
removing unlabelled cells or cells without a clear label.
We also need to compute log-expression values for use in `SingleR()`.


```r
library(scRNAseq)
sceM <- MuraroPancreasData()

sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 

library(scater)
sceM <- logNormCounts(sceM)
```


```r
# Examining the distribution of labels in this reference.
table(sceM$label)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         219         812         448         193         245          21 
##     epsilon mesenchymal          pp 
##           3          80         101
```

We load the @grun2016denovo dataset as our test,
applying some basic quality control to remove low-quality cells in some of the batches
(see [here](https://osca.bioconductor.org/grun-human-pancreas-cel-seq2.html#quality-control-8) for details).
Technically speaking, this does not need log-expression values but we compute them anyway for convenience.


```r
sceG <- GrunPancreasData()

sceG <- addPerCellQC(sceG)
qc <- quickPerCellQC(colData(sceG), 
    percent_subsets="altexps_ERCC_percent",
    batch=sceG$donor,
    subset=sceG$donor %in% c("D17", "D7", "D2"))
sceG <- sceG[,!qc$discard]

sceG <- logNormCounts(sceG)
```


```r
ncol(sceG)
```

```
## [1] 1064
```

## Applying the annotation

We apply `SingleR()` with Wilcoxon rank sum test-based marker detection to annotate the Grun dataset.


```r
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
```

We examine the distribution of predicted labels:


```r
table(pred.grun$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         277         203         181          50         306           5 
##     epsilon mesenchymal          pp 
##           1          22          19
```

## Diagnostics

## Session information {-}


```r
prettySessionInfo()
```

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
 [1] SingleR_1.3.4               scater_1.17.0              
 [3] ggplot2_3.3.0               scRNAseq_2.3.0             
 [5] SingleCellExperiment_1.11.1 SummarizedExperiment_1.19.4
 [7] DelayedArray_0.15.1         matrixStats_0.56.0         
 [9] Biobase_2.49.0              GenomicRanges_1.41.1       
[11] GenomeInfoDb_1.25.0         IRanges_2.23.4             
[13] S4Vectors_0.27.6            BiocGenerics_0.35.2        
[15] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] bitops_1.0-6                  bit64_0.9-7                  
 [3] httr_1.4.1                    tools_4.0.0                  
 [5] R6_2.4.1                      irlba_2.3.3                  
 [7] vipor_0.4.5                   DBI_1.1.0                    
 [9] colorspace_1.4-1              withr_2.2.0                  
[11] gridExtra_2.3                 tidyselect_1.1.0             
[13] processx_3.4.2                bit_1.1-15.2                 
[15] curl_4.3                      compiler_4.0.0               
[17] graph_1.67.0                  BiocNeighbors_1.7.0          
[19] bookdown_0.19                 scales_1.1.1                 
[21] callr_3.4.3                   rappdirs_0.3.1               
[23] stringr_1.4.0                 digest_0.6.25                
[25] rmarkdown_2.1                 XVector_0.29.0               
[27] pkgconfig_2.0.3               htmltools_0.4.0              
[29] limma_3.45.0                  dbplyr_1.4.3                 
[31] fastmap_1.0.1                 rlang_0.4.6                  
[33] RSQLite_2.2.0                 shiny_1.4.0.2                
[35] DelayedMatrixStats_1.11.0     BiocParallel_1.23.0          
[37] dplyr_0.8.5                   RCurl_1.98-1.2               
[39] magrittr_1.5                  BiocSingular_1.5.0           
[41] scuttle_0.99.6                GenomeInfoDbData_1.2.3       
[43] Matrix_1.2-18                 Rcpp_1.0.4.6                 
[45] ggbeeswarm_0.6.0              munsell_0.5.0                
[47] viridis_0.5.1                 lifecycle_0.2.0              
[49] edgeR_3.31.0                  stringi_1.4.6                
[51] yaml_2.2.1                    zlibbioc_1.35.0              
[53] BiocFileCache_1.13.0          AnnotationHub_2.21.0         
[55] grid_4.0.0                    blob_1.2.1                   
[57] dqrng_0.2.1                   promises_1.1.0               
[59] ExperimentHub_1.15.0          crayon_1.3.4                 
[61] lattice_0.20-41               locfit_1.5-9.4               
[63] CodeDepends_0.6.5             knitr_1.28                   
[65] ps_1.3.3                      pillar_1.4.4                 
[67] igraph_1.2.5                  codetools_0.2-16             
[69] XML_3.99-0.3                  glue_1.4.1                   
[71] BiocVersion_3.12.0            evaluate_0.14                
[73] scran_1.17.0                  BiocManager_1.30.10          
[75] vctrs_0.3.0                   httpuv_1.5.2                 
[77] gtable_0.3.0                  purrr_0.3.4                  
[79] assertthat_0.2.1              xfun_0.13                    
[81] rsvd_1.0.3                    mime_0.9                     
[83] xtable_1.8-4                  later_1.0.0                  
[85] viridisLite_0.3.0             tibble_3.0.1                 
[87] AnnotationDbi_1.51.0          beeswarm_0.2.3               
[89] memoise_1.1.0                 statmod_1.4.34               
[91] ellipsis_0.3.1                interactiveDisplayBase_1.27.0
```
</div>
