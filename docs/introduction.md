---
output:
  html_document
bibliography: ref.bib
---

# (PART) Basic usage {-}

# Introduction

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

## Motivation

The Bioconductor package *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* implements an automatic annotation method 
for single-cell RNA sequencing (scRNA-seq) data [@aran2019reference].
Given a reference dataset of samples (single-cell or bulk) with known labels, 
it assigns those labels to new cells from a test dataset based on similarities in their expression profiles.
This provides a convenient way of transferring biological knowledge across datasets,
allowing users to leverage the domain expertise implicit in the creation of each reference.
The most common application of *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* involves predicting cell type (or "state", or "kind") in a new dataset,
a process that is facilitated by the availability of built-in references and compatibility with user-supplied datasets.
In this manner, the burden of manually interpreting clusters and defining marker genes only has to be done once, for the reference dataset, and this knowledge can be propagated to new datasets in an automated manner.

## Method description

*[SingleR](https://bioconductor.org/packages/3.12/SingleR)* can be considered a robust variant of nearest-neighbors classification,
with some tweaks to improve resolution for closely related labels.
For each test cell:

1. We compute the Spearman correlation between its expression profile and that of each reference sample.
The use of Spearman's correlation provides a measure of robustness to batch effects across datasets.
The calculation only uses the union of marker genes identified by pairwise comparisons between labels in the reference data,
so as to improve resolution of separation between labels.
2. We define the per-label score as a fixed quantile (by default, 0.8) of the correlations across all samples with that label.
This accounts for differences in the number of reference samples for each label, 
which interferes with simpler flavors of nearest neighbor classification;
it also avoids penalizing classifications to heterogeneous labels by only requiring a good match to a minority of samples.
3. We repeat the score calculation for all labels in the reference dataset.
The label with the highest score is used as *[SingleR](https://bioconductor.org/packages/3.12/SingleR)*'s prediction for this cell.
4. We optionally perform a fine-tuning step to improve resolution between closely related labels.
The reference dataset is subsetted to only include labels with scores close to the maximum;
scores are recomputed using only marker genes for the subset of labels, thus focusing on the most relevant features;
and this process is iterated until only one label remains.

## Quick start

We will demonstrate the use of `SingleR()` on a classic 10X Genomics dataset [@zheng2017massively]
with the built-in Human Primary Cell Atlas dataset [@hpcaRef] as the reference.


```r
# Loading test data.
library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")

# Loading reference data with Ensembl annotations.
library(SingleR)
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

# Performing predictions.
predictions <- SingleR(test=new.data, assay.type.test=1, 
    ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
```

```
## 
##           B_cell              CMP               DC              GMP 
##              606                8                1                2 
##         Monocyte          NK_cell        Platelets Pre-B_cell_CD34- 
##             1164              217                3               46 
##          T_cells 
##             2293
```

And that's it, really.

## Where to get help

Questions on the general use of *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* should be posted to 
the [Bioconductor support site](https://support.bioconductor.org).
Please send requests for general assistance and advice to the
support site rather than to the individual authors.
Bug reports or feature requests should be made to the [GitHub repository](https://github.com/LTLA/SingleR);
well-considered suggestions for improvements are always welcome.

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
 [1] ensembldb_2.13.1            AnnotationFilter_1.13.0    
 [3] GenomicFeatures_1.41.0      AnnotationDbi_1.51.0       
 [5] SingleR_1.3.5               TENxPBMCData_1.7.0         
 [7] HDF5Array_1.17.0            rhdf5_2.33.3               
 [9] SingleCellExperiment_1.11.4 SummarizedExperiment_1.19.5
[11] DelayedArray_0.15.1         matrixStats_0.56.0         
[13] Biobase_2.49.0              GenomicRanges_1.41.5       
[15] GenomeInfoDb_1.25.1         IRanges_2.23.9             
[17] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[19] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] ProtGenerics_1.21.0           bitops_1.0-6                 
 [3] bit64_0.9-7                   progress_1.2.2               
 [5] httr_1.4.1                    tools_4.0.0                  
 [7] R6_2.4.1                      irlba_2.3.3                  
 [9] lazyeval_0.2.2                DBI_1.1.0                    
[11] rhdf5filters_1.1.0            prettyunits_1.1.1            
[13] tidyselect_1.1.0              processx_3.4.2               
[15] bit_1.1-15.2                  curl_4.3                     
[17] compiler_4.0.0                graph_1.67.1                 
[19] BiocNeighbors_1.7.0           rtracklayer_1.49.3           
[21] bookdown_0.19                 askpass_1.1                  
[23] callr_3.4.3                   rappdirs_0.3.1               
[25] Rsamtools_2.5.1               stringr_1.4.0                
[27] digest_0.6.25                 rmarkdown_2.2                
[29] XVector_0.29.2                pkgconfig_2.0.3              
[31] htmltools_0.4.0               dbplyr_1.4.4                 
[33] fastmap_1.0.1                 rlang_0.4.6                  
[35] RSQLite_2.2.0                 shiny_1.4.0.2                
[37] DelayedMatrixStats_1.11.0     generics_0.0.2               
[39] BiocParallel_1.23.0           dplyr_1.0.0                  
[41] RCurl_1.98-1.2                magrittr_1.5                 
[43] BiocSingular_1.5.0            GenomeInfoDbData_1.2.3       
[45] Matrix_1.2-18                 Rcpp_1.0.4.6                 
[47] Rhdf5lib_1.11.2               lifecycle_0.2.0              
[49] stringi_1.4.6                 yaml_2.2.1                   
[51] zlibbioc_1.35.0               BiocFileCache_1.13.0         
[53] AnnotationHub_2.21.0          grid_4.0.0                   
[55] blob_1.2.1                    promises_1.1.1               
[57] ExperimentHub_1.15.0          crayon_1.3.4                 
[59] lattice_0.20-41               beachmat_2.5.0               
[61] Biostrings_2.57.2             hms_0.5.3                    
[63] CodeDepends_0.6.5             knitr_1.28                   
[65] ps_1.3.3                      pillar_1.4.4                 
[67] biomaRt_2.45.0                codetools_0.2-16             
[69] XML_3.99-0.3                  glue_1.4.1                   
[71] BiocVersion_3.12.0            evaluate_0.14                
[73] BiocManager_1.30.10           vctrs_0.3.1                  
[75] httpuv_1.5.4                  openssl_1.4.1                
[77] purrr_0.3.4                   assertthat_0.2.1             
[79] xfun_0.14                     rsvd_1.0.3                   
[81] mime_0.9                      xtable_1.8-4                 
[83] later_1.1.0.1                 tibble_3.0.1                 
[85] GenomicAlignments_1.25.3      memoise_1.1.0                
[87] ellipsis_0.3.1                interactiveDisplayBase_1.27.5
```
</div>
