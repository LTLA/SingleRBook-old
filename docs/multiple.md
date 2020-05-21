---
bibliography: ref.bib
---

# (PART) Advanced usage {-}

# Using multiple references

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

In some cases, we may wish to use multiple references for annotation of a test dataset.
This yields a more comprehensive set of cell types that are not covered by any individual reference, 
especially when differences in the resolution are considered.
However, it is not trivial due to the presence of batch effects across references
(from differences in technology, experimental protocol or the biological system)
as well as differences in the annotation vocabulary between investigators.

Several strategies are available to combine inferences from multiple references:

- using reference-specific labels in a combined reference
- using harmonized labels in a combined reference
- combining scores across multiple references

This chapter discusses the various strengths and weaknesses of each strategy
and provides some practical demonstrations of each.
Here, we will use the HPCA and BlueprintEncode datasets as our references
and (yet another) PBMC dataset as the test.


```r
library(TENxPBMCData)
pbmc <- TENxPBMCData("pbmc8k")

library(SingleR)
hpca <- HumanPrimaryCellAtlasData(ensembl=TRUE)
bpe <- BlueprintEncodeData(ensembl=TRUE)
```

## Using reference-specific labels

In this strategy, each label is defined in the context of its reference dataset.
This means that a label - say, "B cell" - in reference dataset X is 
considered to be different from a "B cell" label in reference dataset Y.
Use of reference-specific labels is most appropriate if there are relevant biological differences between the references;
for example, if one reference is concerned with healthy tissue while the other reference considers diseased tissue,
it can be helpful to distinguish between the same cell type in different biological contexts.

We can easily implement this approach by combining the expression matrices together 
and pasting the reference name onto the corresponding character vector of labels. 
This modification ensures that the downstream `SingleR()` call
will treat each label-reference combination as a distinct entity.


```r
hpca2 <- hpca
hpca2$label.main <- paste0("HPCA.", hpca2$label.main)

bpe2 <- bpe
bpe2$label.main <- paste0("BPE.", bpe2$label.main)

shared <- intersect(rownames(hpca2), rownames(bpe2))
combined <- cbind(hpca2[shared,], bpe2[shared,])
```

It is then straightforward to perform annotation with the usual methods.


```r
com.res1 <- SingleR(pbmc, ref=combined, labels=combined$label.main, assay.type.test=1)
table(com.res1$labels)
```

```
## 
##      BPE.B-cells BPE.CD4+ T-cells BPE.CD8+ T-cells          BPE.HSC 
##             1178             1708             2656               20 
##    BPE.Monocytes     BPE.NK cells  HPCA.HSC_-G-CSF   HPCA.Platelets 
##             2349              460                1                7 
##     HPCA.T_cells 
##                2
```

However, this strategy identifies markers by directly comparing expression values across references,
meaning that the marker set is likely to contain genes responsible for uninteresting batch effects. 
This will increase noise during the calculation of the score in each reference, 
possibly leading to a loss of precision and a greater risk of technical variation dominating the classification results.
It also complicates interpretation as the cell type is always qualified by its reference of origin in the results.

## Using harmonized labels

This strategy also involves combining the reference datasets into a single matrix 
but the labels are now harmonized so that the same cell type is given the same label across references.
This allows feature selection methods to identify robust sets of label-specific markers 
that are more likely to generalize to other datasets.
It also simplifies interpretation as there is no need to worry about the reference of origin for each label.

Many of the *[SingleR](https://bioconductor.org/packages/3.12/SingleR)* datasets already have their labels 
mapped to the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl).
This provides a standard vocabulary to refer to the same cell type across multiple references.
To simplify interpretation, we set `cell.ont="nonna"` to remove all samples that could not be mapped to the ontology.


```r
hpca.ont <- HumanPrimaryCellAtlasData(ensembl=TRUE, cell.ont="nonna")
bpe.ont <- BlueprintEncodeData(ensembl=TRUE, cell.ont="nonna")

shared <- intersect(rownames(hpca.ont), rownames(bpe.ont))
combined.ont <- cbind(hpca.ont[shared,], bpe.ont[shared,])
combined.ont$source <- rep(c("HPCA", "BPE"), c(ncol(hpca.ont), ncol(bpe.ont)))

# Showing the top 10 most frequent terms:
tab <- table(combined.ont$label.ont, combined.ont$source)
head(tab[order(rowSums(tab), decreasing=TRUE),])
```

```
##             
##              BPE HPCA
##   CL:0000235  18   83
##   CL:0000576  16   57
##   CL:0000451   1   57
##   CL:0000134   0   55
##   CL:0000775  23   21
##   CL:0002618   0   42
```

We perform annotation with `SingleR()` as previously described.
We set the `block=` argument in `de.args=` so as to indicate that DE genes should only be performed _within_ each reference,
and then the statistics merged _across_ references to identify the top markers.
This ensures that we do not directly compare expression values across references,
which reduces the susceptibility of the marker detection to batch effects.


```r
# TODO: add blocking mode for default marker detection.
com.res2 <- SingleR(pbmc, ref=combined.ont, labels=combined.ont$label.ont, 
    assay.type.test=1)
tab <- table(com.res2$labels)
head(sort(tab, decreasing=TRUE))
```

```
## 
## CL:0000624 CL:0000576 CL:0000788 CL:0000815 CL:0000623 CL:0000625 
##       2571       2200        899        681        646        426
```

The most obvious problem with this approach is that it assumes that harmonized labels are available.
This is not always the case due to the use of different (often author-specific!) naming schemes,
which requires at best semi-automated mapping of the author-derived labels to the standard vocabulary.
The mapping process also runs the risk of discarding relevant information about the biological status
(e.g., activation status, disease condition) if there is no obvious counterpart for that state in the ontology.

## Comparing scores across references

The final strategy - and the default approach implemented in `SingleR()` -
involves performing classification separately within each reference, 
and then collating the results to choose the label with the highest score across references. 
This is a relatively expedient approach that avoids the need for explicit harmonization 
while also reducing exposure to reference-specific batch effects.

Use of this approach simply involves passing multiple objects to the `ref=` and `label=` argument in `SingleR()`.
This instructs the function to annotate the test dataset with each reference individually;
it then collects the best labels for each cell across all references and 
identifies the overall best-scoring label as the final prediction for that cell.
The second step requires a recomputation of scores across a subset of relevant markers
to ensure that these scores are comparable across references.


```r
com.res3 <- SingleR(test = pbmc, assay.type.test=1,
    ref = list(BPE=bpe, HPCA=hpca), 
    labels = list(bpe$label.main, hpca$label.main))
table(com.res3$labels)
```

```
## 
##           B_cell          B-cells     CD4+ T-cells     CD8+ T-cells 
##               14             1170             1450             2936 
##              GMP              HSC         Monocyte        Monocytes 
##                1               22              753             1560 
##         NK cells          NK_cell        Platelets Pre-B_cell_CD34- 
##              372               10                9               16 
##          T_cells 
##               68
```

The main appeal of this approach lies in the fact that it is based on the results of annotation with individual references.
This avoids batch effects from comparing expression values across references;
it reduces the need for any coordination in the label scheme between references;
and simultaneously provides the per-reference annotations in the results.
The last feature is particularly useful as it allows for more detailed diagnostics, troubleshooting and further analysis.


```r
head(com.res3$orig.results$BPE$labels)
```

```
## [1] "B-cells"      "Monocytes"    "CD8+ T-cells" "CD8+ T-cells" "Monocytes"   
## [6] "Monocytes"
```

```r
head(com.res3$orig.results$HPCA$labels)
```

```
## [1] "B_cell"   "Monocyte" "T_cells"  "T_cells"  "Monocyte" "Monocyte"
```

The main downside is that it is somewhat suboptimal if there are many reference-specific labels, 
as markers are not identified with the aim of distinguishing a label in one reference from another label in another reference.
The lack of harmonization in the labels also complicates interpretation of the results,
though this can be addressed in the same manner as described above (i.e., replacing `label.main` with `label.ont`).

## Session info {-}

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
 [5] SingleR_1.3.4               TENxPBMCData_1.7.0         
 [7] HDF5Array_1.17.0            rhdf5_2.33.0               
 [9] SingleCellExperiment_1.11.1 SummarizedExperiment_1.19.4
[11] DelayedArray_0.15.1         matrixStats_0.56.0         
[13] Biobase_2.49.0              GenomicRanges_1.41.1       
[15] GenomeInfoDb_1.25.0         IRanges_2.23.4             
[17] S4Vectors_0.27.6            BiocGenerics_0.35.2        
[19] BiocStyle_2.17.0            rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] ProtGenerics_1.21.0           bitops_1.0-6                 
 [3] bit64_0.9-7                   progress_1.2.2               
 [5] httr_1.4.1                    tools_4.0.0                  
 [7] R6_2.4.1                      irlba_2.3.3                  
 [9] lazyeval_0.2.2                DBI_1.1.0                    
[11] prettyunits_1.1.1             tidyselect_1.1.0             
[13] processx_3.4.2                bit_1.1-15.2                 
[15] curl_4.3                      compiler_4.0.0               
[17] graph_1.67.0                  BiocNeighbors_1.7.0          
[19] rtracklayer_1.49.1            bookdown_0.19                
[21] askpass_1.1                   callr_3.4.3                  
[23] rappdirs_0.3.1                Rsamtools_2.5.0              
[25] stringr_1.4.0                 digest_0.6.25                
[27] rmarkdown_2.1                 XVector_0.29.0               
[29] pkgconfig_2.0.3               htmltools_0.4.0              
[31] dbplyr_1.4.3                  fastmap_1.0.1                
[33] rlang_0.4.6                   RSQLite_2.2.0                
[35] shiny_1.4.0.2                 DelayedMatrixStats_1.11.0    
[37] BiocParallel_1.23.0           dplyr_0.8.5                  
[39] RCurl_1.98-1.2                magrittr_1.5                 
[41] BiocSingular_1.5.0            GenomeInfoDbData_1.2.3       
[43] Matrix_1.2-18                 Rcpp_1.0.4.6                 
[45] Rhdf5lib_1.11.0               lifecycle_0.2.0              
[47] stringi_1.4.6                 yaml_2.2.1                   
[49] zlibbioc_1.35.0               BiocFileCache_1.13.0         
[51] AnnotationHub_2.21.0          grid_4.0.0                   
[53] blob_1.2.1                    promises_1.1.0               
[55] ExperimentHub_1.15.0          crayon_1.3.4                 
[57] lattice_0.20-41               beachmat_2.5.0               
[59] Biostrings_2.57.0             hms_0.5.3                    
[61] CodeDepends_0.6.5             knitr_1.28                   
[63] ps_1.3.3                      pillar_1.4.4                 
[65] codetools_0.2-16              biomaRt_2.45.0               
[67] XML_3.99-0.3                  glue_1.4.1                   
[69] BiocVersion_3.12.0            evaluate_0.14                
[71] BiocManager_1.30.10           vctrs_0.3.0                  
[73] httpuv_1.5.2                  openssl_1.4.1                
[75] purrr_0.3.4                   assertthat_0.2.1             
[77] xfun_0.13                     rsvd_1.0.3                   
[79] mime_0.9                      xtable_1.8-4                 
[81] later_1.0.0                   tibble_3.0.1                 
[83] GenomicAlignments_1.25.0      memoise_1.1.0                
[85] ellipsis_0.3.1                interactiveDisplayBase_1.27.0
```
</div>
