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
a process that is facilitated by the availability of curated references and compatibility with user-supplied datasets.
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

We will demonstrate the use of `SingleR()` on a well-known 10X Genomics dataset [@zheng2017massively]
with the Human Primary Cell Atlas dataset [@hpcaRef] as the reference.


```r
# Loading test data.
library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")

# Loading reference data with Ensembl annotations.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

# Performing predictions.
library(SingleR)
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
 [1] SingleR_1.3.6               ensembldb_2.13.1           
 [3] AnnotationFilter_1.13.0     GenomicFeatures_1.41.0     
 [5] AnnotationDbi_1.51.1        celldex_0.99.1             
 [7] TENxPBMCData_1.7.0          HDF5Array_1.17.3           
 [9] rhdf5_2.33.4                SingleCellExperiment_1.11.6
[11] SummarizedExperiment_1.19.5 DelayedArray_0.15.6        
[13] matrixStats_0.56.0          Matrix_1.2-18              
[15] Biobase_2.49.0              GenomicRanges_1.41.5       
[17] GenomeInfoDb_1.25.5         IRanges_2.23.10            
[19] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[21] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] ProtGenerics_1.21.0           bitops_1.0-6                 
 [3] bit64_0.9-7                   progress_1.2.2               
 [5] httr_1.4.1                    tools_4.0.0                  
 [7] irlba_2.3.3                   R6_2.4.1                     
 [9] lazyeval_0.2.2                DBI_1.1.0                    
[11] rhdf5filters_1.1.1            tidyselect_1.1.0             
[13] prettyunits_1.1.1             processx_3.4.2               
[15] bit_1.1-15.2                  curl_4.3                     
[17] compiler_4.0.0                graph_1.67.1                 
[19] BiocNeighbors_1.7.0           rtracklayer_1.49.3           
[21] bookdown_0.20                 askpass_1.1                  
[23] callr_3.4.3                   rappdirs_0.3.1               
[25] Rsamtools_2.5.3               stringr_1.4.0                
[27] digest_0.6.25                 rmarkdown_2.3                
[29] XVector_0.29.3                pkgconfig_2.0.3              
[31] htmltools_0.5.0               dbplyr_1.4.4                 
[33] fastmap_1.0.1                 rlang_0.4.6                  
[35] RSQLite_2.2.0                 shiny_1.5.0                  
[37] DelayedMatrixStats_1.11.1     generics_0.0.2               
[39] BiocParallel_1.23.0           dplyr_1.0.0                  
[41] BiocSingular_1.5.0            RCurl_1.98-1.2               
[43] magrittr_1.5                  GenomeInfoDbData_1.2.3       
[45] Rcpp_1.0.4.6                  Rhdf5lib_1.11.2              
[47] lifecycle_0.2.0               stringi_1.4.6                
[49] yaml_2.2.1                    zlibbioc_1.35.0              
[51] BiocFileCache_1.13.0          AnnotationHub_2.21.1         
[53] grid_4.0.0                    blob_1.2.1                   
[55] promises_1.1.1                ExperimentHub_1.15.0         
[57] crayon_1.3.4                  lattice_0.20-41              
[59] beachmat_2.5.0                Biostrings_2.57.2            
[61] hms_0.5.3                     CodeDepends_0.6.5            
[63] knitr_1.29                    ps_1.3.3                     
[65] pillar_1.4.4                  codetools_0.2-16             
[67] biomaRt_2.45.1                XML_3.99-0.3                 
[69] glue_1.4.1                    BiocVersion_3.12.0           
[71] evaluate_0.14                 BiocManager_1.30.10          
[73] vctrs_0.3.1                   httpuv_1.5.4                 
[75] openssl_1.4.2                 purrr_0.3.4                  
[77] assertthat_0.2.1              xfun_0.15                    
[79] rsvd_1.0.3                    mime_0.9                     
[81] xtable_1.8-4                  later_1.1.0.1                
[83] tibble_3.0.1                  GenomicAlignments_1.25.3     
[85] memoise_1.1.0                 ellipsis_0.3.1               
[87] interactiveDisplayBase_1.27.5
```
</div>
