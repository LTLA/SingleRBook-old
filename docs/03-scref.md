# Using single-cell references



## Running the algorithm

Here, we will use two human pancreas datasets from the *[scRNAseq](https://bioconductor.org/packages/3.12/scRNAseq)* package.
The aim is to use one pre-labelled dataset to annotate the other unlabelled dataset.
First, we set up the @muraro2016singlecell dataset to be our reference.


```r
library(scRNAseq)
sceM <- MuraroPancreasData()

# One should normally do cell-based quality control at this point, but for
# brevity's sake, we will just remove the unlabelled libraries here.
sceM <- sceM[,!is.na(sceM$label)]

library(scater)
sceM <- logNormCounts(sceM)
```

We then set up our test dataset from @grun2016denovo.
To speed up this demonstration, we will subset to the first 100 cells.


```r
sceG <- GrunPancreasData()
sceG <- sceG[,colSums(counts(sceG)) > 0] # Remove libraries with no counts.
sceG <- logNormCounts(sceG) 
sceG <- sceG[,1:100]
```

We then run `SingleR()` as described previously but with a marker detection mode that considers the variance of expression across cells.
Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels.
This is slower but more appropriate for single-cell data compared to the default marker detection algorithm (which may fail for low-coverage data where the median is frequently zero).


```r
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
table(pred.grun$labels)
```

```
## 
## acinar   beta  delta   duct 
##     53      4      2     41
```

## Defining custom markers

Users can also construct their own marker lists with any DE testing machinery.
For example, we can perform pairwise $t$-tests using methods from *[scran](https://bioconductor.org/packages/3.12/scran)* and obtain the top 10 marker genes from each pairwise comparison.


```r
library(scran)
out <- pairwiseTTests(logcounts(sceM), sceM$label, direction="up")
markers <- getTopMarkers(out$statistics, out$pairs, n=10)
```

We then supply these genes to `SingleR()` directly via the `genes=` argument.
A more focused gene set also allows annotation to be performed more quickly compared to the default approach.


```r
pred.grun2 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=markers)
table(pred.grun2$labels)
```

```
## 
##  acinar    beta   delta    duct      pp unclear 
##      59       4       1      34       1       1
```

In some cases, markers may only be available for specific labels rather than for pairwise comparisons between labels.
This is accommodated by supplying a named list of character vectors to `genes`.
Note that this is likely to be less powerful than the list-of-lists approach as information about pairwise differences is discarded.


```r
label.markers <- lapply(markers, unlist, recursive=FALSE)
pred.grun3 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=label.markers)
table(pred.grun$labels, pred.grun3$labels)
```

```
##         
##          acinar beta delta duct pp
##   acinar     51    0     0    2  0
##   beta        0    4     0    0  0
##   delta       0    0     1    0  1
##   duct        2    0     0   39  0
```

## Pseudo-bulk aggregation

Single-cell reference datasets provide a like-for-like comparison to our test datasets, yielding a more accurate classification of the cells in the latter (hopefully).
However, there are frequently many more samples in single-cell references compared to bulk references, increasing the computational work involved in classification.
We avoid this by aggregating cells into one "pseudo-bulk" sample per label (e.g., by averaging across log-expression values) and using those as the reference, which allows us to achieve the same efficiency as the use of bulk references.

The obvious cost of this approach is that we discard potentially useful information about the distribution of cells within each label.
Cells that belong to a heterogeneous population may not be correctly assigned if they are far from the population center.
We attempt to preserve some of this information by using $k$-means clustering within each cell to create pseudo-bulk samples that are representative of a particular region of the expression space (i.e., vector quantization).
We create $\sqrt{N}$ clusters given a label with $N$ cells, which provides a reasonable compromise between reducing computational work and preserving the label's internal distribution.

This aggregation approach is implemented in the `aggregateReferences` function, which is shown in action below for the @muraro2016singlecell dataset.
The function returns a `SummarizedExperiment` object containing the pseudo-bulk expression profiles and the corresponding labels.


```r
set.seed(100) # for the k-means step.
aggr <- aggregateReference(sceM, labels=sceM$label)
aggr
```

```
## class: SummarizedExperiment 
## dim: 19059 116 
## metadata(0):
## assays(1): logcounts
## rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(0):
## colnames(116): alpha.1 alpha.2 ... mesenchymal.8 epsilon.1
## colData names(1): label
```

The resulting `SummarizedExperiment` can then be used as a reference in `SingleR()`.


```r
pred.aggr <- SingleR(sceG, aggr, labels=aggr$label)
table(pred.aggr$labels)
```

```
## 
## acinar   beta  delta   duct 
##     53      4      1     42
```



