---
bibliography: ref.bib
---

# Using the built-in references

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

*[SingleR](https://bioconductor.org/packages/3.12/SingleR)* detects markers in a pairwise manner between labels in the reference dataset.
Specifically, for each label of interest, it performs pairwise comparisons to every other label in the reference
and identifies the genes that are upregulated in the label of interest for each comparison.
The initial score calculation is then performed on the union of marker genes across all comparisons for all label.
This approach ensures that the selected subset of features will contain genes that distinguish each label from any other label.
(In contrast, other approaches that treat the "other" labels as a single group do not offer this guarantee;
see [here](https://osca.bioconductor.org/marker-detection.html#standard-application) for a discussion.)
It also allows the fine-tuning step to aggressively improve resolution by only using marker genes 
from comparisons where both labels have scores close to the maximum.

The original ("classic") marker detection algorithm used in @aran2019reference identified marker genes 
based on their log-fold changes in each pairwise comparison.
Specifically, it used the genes with the largest positive differences in the per-label median log-expression values between labels.
The number of genes taken from each pairwise comparison was defined as $500 (\frac{2}{3})^{\log_{2}(n)}$,
where $n$ is the number of unique labels in the reference;
this aimed to reduce the number of genes (and thus the computational time) as the number of labels and pairwise comparisons increased.
This classic mode is primarily intended for reference datasets that have little or no replication,
a description that covers most of the built-in references 
and precludes more complicated marker detection procedures (Chapter \@ref(using-single-cell-references)).

## Annotation with default marker detection 

For demonstration purposes, we will use the @grun2016denovo haematopoietic stem cell (HSC)
dataset from the *[scRNAseq](https://bioconductor.org/packages/3.12/scRNAseq)* package.
The `GrunHSCData()` function conveniently returns a `SingleCellExperiment` 
object containing the count matrix for this dataset.


```r
library(scRNAseq)
sce <- GrunHSCData(ensembl=TRUE)
sce
```

```
## class: SingleCellExperiment 
## dim: 21817 1915 
## metadata(0):
## assays(1): counts
## rownames(21817): ENSMUSG00000109644 ENSMUSG00000007777 ...
##   ENSMUSG00000055670 ENSMUSG00000039068
## rowData names(3): symbol chr originalName
## colnames(1915): JC4_349_HSC_FE_S13_ JC4_350_HSC_FE_S13_ ...
##   JC48P6_1203_HSC_FE_S8_ JC48P6_1204_HSC_FE_S8_
## colData names(2): sample protocol
## reducedDimNames(0):
## altExpNames(0):
```

Our plan is to annotate each cell with the built-in ImmGen reference dataset [@ImmGenRef].
Calling the `ImmGenData()` function returns a `SummarizedExperiment` object 
containing a matrix of log-expression values with sample-level labels.
We set `ensembl=TRUE` to match the reference's gene annotation with that in the `sce` object
(the default behavior is to use the gene symbol).


```r
library(SingleR)
immgen <- ImmGenData(ensembl=TRUE)
immgen
```

```
## class: SummarizedExperiment 
## dim: 21352 830 
## metadata(0):
## assays(1): logcounts
## rownames(21352): ENSMUSG00000079681 ENSMUSG00000066372 ...
##   ENSMUSG00000034640 ENSMUSG00000036940
## rowData names(0):
## colnames(830):
##   GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL
##   GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL ...
##   GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL
##   GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
## colData names(3): label.main label.fine label.ont
```

Technically speaking, each built-in dataset actually has three sets of labels that primarily differ in their resolution.
For the purposes of this demonstration, we will use the "fine" labels in the `label.fine` metadata field,
which represents the highest resolution of annotation available in ImmGen.


```r
head(immgen$label.fine)
```

```
## [1] "Macrophages (MF.11C-11B+)" "Macrophages (MF.11C-11B+)"
## [3] "Macrophages (MF.11C-11B+)" "Macrophages (MF.ALV)"     
## [5] "Macrophages (MF.ALV)"      "Macrophages (MF.ALV)"
```

Annotation is then a simple matter of calling `SingleR()` on our test (Grun) dataset and the reference (ImmGen) dataset,
leaving the default of `de.method="classic"` to use the original marker detection scheme.
This applies the algorithm described in Section \@ref(method-description),
returning a `DataFrame` where each row contains prediction results for a single cell in the `sce` object.
Labels are provided before fine-tuning (`first.labels`), after fine-tuning (`labels`) and after pruning (`pruned.labels`);
some of the other fields are discussed in more detail in Chapter \@ref(annotation-diagnostics).


```r
# See 'Choices of assay data' for 'assay.type.test=' explanation.
pred <- SingleR(test = sce, ref = immgen, 
    labels = immgen$label.fine, assay.type.test=1)
colnames(pred)
```

```
## [1] "scores"        "first.labels"  "tuning.scores" "labels"       
## [5] "pruned.labels"
```

## Interaction with quality control

Upon examining the distribution of assigned labels, we see that many of them are related to stem cells.
However, there are quite a large number of more differentiated labels mixed in,
which is not what we expect from a sorted population of HSCs.


```r
head(sort(table(pred$labels), decreasing=TRUE))
```

```
## 
##   Stem cells (SC.MEP) Neutrophils (GN.ARTH)      Macrophages (MF) 
##                   362                   314                   165 
##  Stem cells (SC.STSL)    B cells (proB.FrA) Stem cells (SC.LT34F) 
##                   142                   121                   103
```



This is probably because - despite what its name might suggest -
the dataset obtained by `GrunHSCData()` actually contains more than HSCs.
If we restrict our analysis to the sorted HSCs (obviously) and remove one low-quality batch
(see [the analysis here](https://osca.bioconductor.org/merged-hcsc.html#quality-control-12) for the rationale)
we can see that the distribution of cell type labels is more similar to what we might expect.
Low-quality cells - typically those with low library sizes and low numbers of detected genes -
lack information for accurate label assignment and often need to be removed to enable interpretation of the results.


```r
actual.hsc <- pred$labels[sce$protocol=="sorted hematopoietic stem cells" & sce$sample!="JC4"]
head(sort(table(actual.hsc), decreasing=TRUE))
```

```
## actual.hsc
##        Stem cells (SC.STSL)       Stem cells (SC.LT34F) 
##                         109                          98 
##       Stem cells (SC.ST34F) Stem cells (SC.CD150-CD48-) 
##                          37                          15 
##          Stem cells (LTHSC)            Stem cells (MLP) 
##                          12                           8
```



Filtering the annotation results in the above manner is possible because `SingleR()` operates independently on each cell.
The annotation is orthogonal to any decisions about the relative quality of the cells in the test dataset;
the same results will be obtained regardless of whether the annotation is performed before or after quality control.
This is logisitically convenient as it means that the annotation does not have to be repeated 
if the quality control scheme (or any other downstream step, like clustering) changes throughout the lifetime of the analysis.

## Choices of assay data

For the reference dataset, the assay matrix _must_ contain log-transformed normalized expression values.
This is because the default marker detection scheme computes log-fold changes by subtracting the medians,
which makes little sense unless the input expression values are already log-transformed.
For alternative schemes, this requirement may be relaxed (e.g., Wilcoxon rank sum tests do not require transformation);
similarly, if pre-defined markers are supplied, no transformation or normalization is necessary.

For the test data, the assay data need not be log-transformed or even (scale) normalized.
This is because `SingleR()` computes Spearman correlations within each cell, 
which is unaffected by monotonic transformations like cell-specific scaling or log-transformation.
It is perfectly satisfactory to provide the raw counts for the test dataset to `SingleR()`,
which is the reason for setting `assay.type.test=1` in our previous `SingleR()` call for the Grun dataset.

The exception to this rule occurs when comparing data from full-length technologies to the built-in references.
The built-in references are constructed to be comparable to unique molecular identifier (UMI) protocols,
where the expression values are less sensitive to differences in gene length.
Thus, when comparing Smart-seq2 test datasets to the built-in references,
better performance can often be achieved by processing the test counts to transcripts-per-million values.

We demonstrate below using another HSC dataset that was generated using the Smart-seq2 protocol [@nestorowa2016singlecell].
Again, we see that most of the predicted labels are related to stem cells, which is comforting.


```r
sce.nest <- NestorowaHSCData()

# Getting the exonic gene lengths.
library(AnnotationHub)
mm.db <- AnnotationHub()[["AH73905"]]
mm.exons <- exonsBy(mm.db, by="gene")
mm.exons <- reduce(mm.exons)
mm.len <- sum(width(mm.exons))

# Computing the TPMs with a simple scaling by gene length.
library(scater)
keep <- intersect(names(mm.len), rownames(sce.nest))
tpm.nest <- calculateTPM(sce.nest[keep,], lengths=mm.len[keep])

# Performing the assignment.
pred <- SingleR(test = tpm.nest, ref = immgen, labels = immgen$label.fine)
head(sort(table(pred$labels), decreasing=TRUE), 10)
```

```
## 
##         Stem cells (SC.MEP)       Stem cells (SC.ST34F) 
##                         409                         357 
##      Stem cells (SC.MPP34F)      Stem cells (SC.CMP.DR) 
##                         329                         298 
##            Stem cells (MLP)            Stem cells (GMP) 
##                         167                         102 
##        Stem cells (SC.STSL)         Stem cells (SC.MDP) 
##                          71                          66 
## Stem cells (SC.CD150-CD48-)       Stem cells (SC.LT34F) 
##                          55                          37
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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] scater_1.17.3               ggplot2_3.3.1              
 [3] AnnotationHub_2.21.0        BiocFileCache_1.13.0       
 [5] dbplyr_1.4.4                SingleR_1.3.5              
 [7] ensembldb_2.13.1            AnnotationFilter_1.13.0    
 [9] GenomicFeatures_1.41.0      AnnotationDbi_1.51.0       
[11] scRNAseq_2.3.4              SingleCellExperiment_1.11.4
[13] SummarizedExperiment_1.19.5 DelayedArray_0.15.1        
[15] matrixStats_0.56.0          Biobase_2.49.0             
[17] GenomicRanges_1.41.5        GenomeInfoDb_1.25.1        
[19] IRanges_2.23.9              S4Vectors_0.27.12          
[21] BiocGenerics_0.35.4         BiocStyle_2.17.0           
[23] rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0              colorspace_1.4-1             
 [3] ellipsis_0.3.1                scuttle_0.99.9               
 [5] XVector_0.29.2                BiocNeighbors_1.7.0          
 [7] bit64_0.9-7                   interactiveDisplayBase_1.27.5
 [9] codetools_0.2-16              knitr_1.28                   
[11] Rsamtools_2.5.1               graph_1.67.1                 
[13] shiny_1.4.0.2                 BiocManager_1.30.10          
[15] compiler_4.0.0                httr_1.4.1                   
[17] assertthat_0.2.1              Matrix_1.2-18                
[19] fastmap_1.0.1                 lazyeval_0.2.2               
[21] later_1.1.0.1                 BiocSingular_1.5.0           
[23] htmltools_0.4.0               prettyunits_1.1.1            
[25] tools_4.0.0                   rsvd_1.0.3                   
[27] gtable_0.3.0                  glue_1.4.1                   
[29] GenomeInfoDbData_1.2.3        dplyr_1.0.0                  
[31] rappdirs_0.3.1                Rcpp_1.0.4.6                 
[33] vctrs_0.3.1                   Biostrings_2.57.2            
[35] ExperimentHub_1.15.0          rtracklayer_1.49.3           
[37] DelayedMatrixStats_1.11.0     xfun_0.14                    
[39] stringr_1.4.0                 ps_1.3.3                     
[41] mime_0.9                      lifecycle_0.2.0              
[43] irlba_2.3.3                   XML_3.99-0.3                 
[45] zlibbioc_1.35.0               scales_1.1.1                 
[47] hms_0.5.3                     promises_1.1.1               
[49] ProtGenerics_1.21.0           yaml_2.2.1                   
[51] curl_4.3                      memoise_1.1.0                
[53] gridExtra_2.3                 biomaRt_2.45.0               
[55] stringi_1.4.6                 RSQLite_2.2.0                
[57] BiocVersion_3.12.0            BiocParallel_1.23.0          
[59] rlang_0.4.6                   pkgconfig_2.0.3              
[61] bitops_1.0-6                  evaluate_0.14                
[63] lattice_0.20-41               purrr_0.3.4                  
[65] GenomicAlignments_1.25.3      CodeDepends_0.6.5            
[67] bit_1.1-15.2                  processx_3.4.2               
[69] tidyselect_1.1.0              magrittr_1.5                 
[71] bookdown_0.19                 R6_2.4.1                     
[73] generics_0.0.2                DBI_1.1.0                    
[75] pillar_1.4.4                  withr_2.2.0                  
[77] RCurl_1.98-1.2                tibble_3.0.1                 
[79] crayon_1.3.4                  rmarkdown_2.2                
[81] viridis_0.5.1                 progress_1.2.2               
[83] grid_4.0.0                    blob_1.2.1                   
[85] callr_3.4.3                   digest_0.6.25                
[87] xtable_1.8-4                  httpuv_1.5.4                 
[89] openssl_1.4.1                 munsell_0.5.0                
[91] beeswarm_0.2.3                viridisLite_0.3.0            
[93] vipor_0.4.5                   askpass_1.1                  
```
</div>
