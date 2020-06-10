# Exploiting the cell ontology

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

As previously discussed in Section \@ref(using-harmonized-labels),
*[SingleR](https://bioconductor.org/packages/3.12/SingleR)* maps the labels in its references to the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl).
The most obvious advantage of doing this is to provide a standardized vocabulary with which to describe cell types,
thus facilitating integrated analyses with multiple references.
However, another useful feature of the Cell Ontology is its hierarchical organization of terms,
allowing us to adjust cell type annotations to the desired resolution.
This represents a more dynamic alternative to the static `label.main` and `label.fine` options in each reference.

## Basic manipulation

We use the *[ontoProc](https://bioconductor.org/packages/3.12/ontoProc)* package to load in the Cell Ontology.
This produces an `ontology_index` object (from the *[ontologyIndex](https://CRAN.R-project.org/package=ontologyIndex)* package)
that we can query for various pieces of information.


```r
# TODO: wrap in utility function.
library(ontoProc)
bfc <- BiocFileCache::BiocFileCache(ask=FALSE)
path <- BiocFileCache::bfcrpath(bfc, "http://purl.obolibrary.org/obo/cl.obo")
cl <- get_ontology(path, extract_tags="everything")
cl
```

```
## Ontology with 2236 terms
## 
## format-version: 1.2
## data-version: cl/2020-05-21/cl-simple.owl
## ontology: cl/cl-simple
## 
## Properties:
## 	id: character
## 	name: character
## 	parents: list
## 	children: list
## 	ancestors: list
## 	obsolete: logical
## 	alt_id: list
## 	comment: character
## 	consider: list
## 	created_by: character
## 	creation_date: character
## 	data-version: list
## 	def: character
## 	develops_from: list
## 	format-version: list
## 	holds_over_chain: list
## 	is_a: list
## 	is_transitive: list
## 	namespace: list
## 	ontology: list
## 	property_value: list
## 	remark: list
## 	replaced_by: list
## 	subset: list
## 	subsetdef: list
## 	synonym: list
## 	synonymtypedef: list
## 	transitive_over: list
## 	xref: list
## Roots:
## 	CL:0000000 - cell
## 	BFO:0000050 - NA
## 	develops_from - develops_from
```

The most immediate use of this object lies in mapping ontology terms to their plain-English descriptions.
We can use this to translate annotations produced by `SingleR()` from the `label.ont` labels into a more interpretable form.
We demonstrate this approach using *[SingleR](https://bioconductor.org/packages/3.12/SingleR)*'s collection of mouse RNA-seq references [@aran2019reference].


```r
head(cl$name) # short name
```

```
##                         CL:0000000                         CL:0000001 
##                             "cell"            "primary cultured cell" 
##                         CL:0000002                         CL:0000003 
## "obsolete immortal cell line cell"                      "native cell" 
##                         CL:0000004                         CL:0000005 
##        "obsolete cell by organism"  "fibroblast neural crest derived"
```

```r
head(cl$def) # longer definition
```

```
##                                                                                                                                                                                                                     CL:0000000 
##                                   "\"A material entity of anatomical origin (part of or deriving from an organism) that has as its parts a maximally connected cell compartment surrounded by a plasma membrane.\" [CARO:mah]" 
##                                                                                                                                                                                                                     CL:0000001 
##                                                                 "\"A cultured cell that is freshly isolated from a organismal source, or derives in culture from such a cell prior to the culture being passaged.\" [ReO:mhb]" 
##                                                                                                                                                                                                                     CL:0000002 
##             "\"OBSOLETE: A cell line cell that is expected to be capable of an unlimited number of divisions, and is thus able to support indefinite growth/propagation in vitro as part of a immortal cell line.\" [ReO:mhb]" 
##                                                                                                                                                                                                                     CL:0000003 
## "\"A cell that is found in a natural setting, which includes multicellular organism cells 'in vivo' (i.e. part of an organism), and unicellular organisms 'in environment' (i.e. part of a natural environment).\" [CARO:mah]" 
##                                                                                                                                                                                                                     CL:0000004 
##                                                                                                                            "\"OBSOLETE: A classification of cells by the organisms within which they are contained.\" [FB:ma]" 
##                                                                                                                                                                                                                     CL:0000005 
##                                                                                                                           "\"Any fibroblast that is deriived from the neural crest.\" [https://orcid.org/0000-0001-5208-3432]"
```

```r
library(SingleR)
ref <- MouseRNAseqData(cell.ont="nonna")
translated <- cl$name[ref$label.ont]
head(translated)
```

```
## CL:0000136 CL:0000136 CL:0000136 CL:0000136 CL:0000136 CL:0000136 
## "fat cell" "fat cell" "fat cell" "fat cell" "fat cell" "fat cell"
```

Another interesting application involves examining the relationship between different terms.
The ontology itself is a directed acyclic graph, so we can can convert it into `graph` object
for advanced queries using the *[igraph](https://CRAN.R-project.org/package=igraph)* package.
Each edge represents an "is a" relationship where each vertex represents a specialized case of the concept of the parent node.


```r
# TODO: wrap in utility function.
parents <- cl$parents
self <- rep(names(parents), lengths(parents))

library(igraph)
g <- make_graph(rbind(unlist(parents), self))
g
```

```
## IGRAPH 5434583 DN-- 2234 3141 -- 
## + attr: name (v/c)
## + edges from 5434583 (vertex names):
##  [1] CL:0000010->CL:0000001 CL:0000000->CL:0000003 CL:0000057->CL:0000005
##  [4] CL:0000101->CL:0000006 CL:0000197->CL:0000006 CL:0002321->CL:0000007
##  [7] CL:0000333->CL:0000008 CL:0000578->CL:0000010 CL:0000333->CL:0000011
## [10] CL:0000034->CL:0000014 CL:0000039->CL:0000014 CL:0000586->CL:0000015
## [13] CL:0000014->CL:0000016 CL:0000015->CL:0000016 CL:0000015->CL:0000017
## [16] CL:0000015->CL:0000018 CL:0000413->CL:0000018 CL:0000408->CL:0000019
## [19] CL:0000015->CL:0000020 CL:0000586->CL:0000021 CL:0000014->CL:0000022
## [22] CL:0000021->CL:0000022 CL:0000021->CL:0000023 CL:0000021->CL:0000024
## + ... omitted several edges
```

One query involves identifying all descendents of a particular term of interest.
This can be useful when searching for a cell type in the presence of variable annotation resolution;
for example, a search for "epithelial cell" can be configured to pick up all child terms 
such as "endothelial cell" and "ependymal cell".


```r
term <- "CL:0000624"
cl$name[term]
```

```
##                        CL:0000624 
## "CD4-positive, alpha-beta T cell"
```

```r
all.kids <- names(subcomponent(g, term))
head(cl$name[all.kids])
```

```
##                                                       CL:0000624 
##                                "CD4-positive, alpha-beta T cell" 
##                                                       CL:0000492 
##                                     "CD4-positive helper T cell" 
##                                                       CL:0001051 
## "CD4-positive, CXCR3-negative, CCR6-negative, alpha-beta T cell" 
##                                                       CL:0000791 
##                                       "mature alpha-beta T cell" 
##                                                       CL:0000792 
##      "CD4-positive, CD25-positive, alpha-beta regulatory T cell" 
##                                                       CL:0000793 
##                "CD4-positive, alpha-beta intraepithelial T cell"
```

Alternatively, we might be interested in the last common ancestor (LCA) for a set of terms.
This is the furthest term - or, in some cases, multiple terms - from the root of the ontology
that is also an ancestor of all of the terms of interest.
We will use this LCA concept in the next section to adjust resolution across multiple references. 


```r
terms <- c("CL:0000624", "CL:0000785", "CL:0000623")
cl$name[terms]
```

```
##                        CL:0000624                        CL:0000785 
## "CD4-positive, alpha-beta T cell"                   "mature B cell" 
##                        CL:0000623 
##             "natural killer cell"
```

```r
# TODO: god, put this in a function somewhere.
all.ancestors <- lapply(terms, subcomponent, graph=g, mode="in")
all.ancestors <- lapply(all.ancestors, names)
common.ancestors <- Reduce(intersect, all.ancestors)

ancestors.of.ancestors <- lapply(common.ancestors, subcomponent, graph=g, mode="in")
ancestors.of.ancestors <- lapply(ancestors.of.ancestors, names)
ancestors.of.ancestors <- mapply(setdiff, ancestors.of.ancestors, common.ancestors) 

latest.common.ancestors <- setdiff(common.ancestors, unlist(ancestors.of.ancestors))
cl$name[latest.common.ancestors]
```

```
##   CL:0000542 
## "lymphocyte"
```

## Adjusting resolution

We can use the ontology graph to adjust the resolution of the reference labels by rolling up overly-specific terms to their LCA.
The `findCommonAncestors()` utility takes a set of terms and returns a list of potential LCAs for various subsets of those terms.
Users can inspect this list to identify LCAs at the desired resolution and then map their descendent terms to those LCAs.


```r
findCommonAncestors <- function(..., g, remove.self=TRUE, names=NULL) {
    terms <- list(...)
    if (is.null(names(terms))) {
        names(terms) <- sprintf("set%i", seq_along(terms))
    }

    all.terms <- unique(unlist(terms))
    all.ancestors <- lapply(all.terms, subcomponent, graph=g, mode="in")
    all.ancestors <- lapply(all.ancestors, names)
    by.ancestor <- split(
        rep(all.terms, lengths(all.ancestors)),
        unlist(all.ancestors)
    )

    # Removing ancestor nodes with the same count as its children.
    available <- names(by.ancestor)
    for (i in available) {
        if (!i %in% names(by.ancestor)) {
            next
        }

        counts <- lengths(by.ancestor)
        cur.ancestors <- subcomponent(g, i, mode="in")
        cur.ancestors <- setdiff(names(cur.ancestors), i)
        drop <- cur.ancestors[counts[i]==counts[cur.ancestors]]
        by.ancestor <- by.ancestor[!names(by.ancestor) %in% drop]
    }

    if (remove.self) {
        by.ancestor <- by.ancestor[lengths(by.ancestor) > 1L]
    }
    by.ancestor <- by.ancestor[order(lengths(by.ancestor))] # most specific terms first.

    # Decorating the output.
    for (i in names(by.ancestor)) {
        current <- by.ancestor[[i]]
        df <- DataFrame(row.names=current)

        curout <- list()
        if (!is.null(names)) {
            curout$name <- unname(names[i])
            df$name <- names[current]
        }

        presence <- list()
        for (b in names(terms)) {
            presence[[b]] <- current %in% terms[[b]]
        }
        df <- cbind(df, do.call(DataFrame, presence))

        curout$descendents <- df
        by.ancestor[[i]] <- curout
    }

    by.ancestor
}

lca <- findCommonAncestors(ref$label.ont, g=g, names=cl$name)
head(lca)
```

```
## $`CL:0000081`
## $`CL:0000081`$name
## [1] "blood cell"
## 
## $`CL:0000081`$descendents
## DataFrame with 2 rows and 2 columns
##                   name      set1
##            <character> <logical>
## CL:0000232 erythrocyte      TRUE
## CL:0000094 granulocyte      TRUE
## 
## 
## $`CL:0000126`
## $`CL:0000126`$name
## [1] "macroglial cell"
## 
## $`CL:0000126`$descendents
## DataFrame with 2 rows and 2 columns
##                       name      set1
##                <character> <logical>
## CL:0000127       astrocyte      TRUE
## CL:0000128 oligodendrocyte      TRUE
## 
## 
## $`CL:0000393`
## $`CL:0000393`$name
## [1] "electrically responsive cell"
## 
## $`CL:0000393`$descendents
## DataFrame with 2 rows and 2 columns
##                           name      set1
##                    <character> <logical>
## CL:0000540              neuron      TRUE
## CL:0000746 cardiac muscle cell      TRUE
## 
## 
## $`CL:0002320`
## $`CL:0002320`$name
## [1] "connective tissue cell"
## 
## $`CL:0002320`$descendents
## DataFrame with 2 rows and 2 columns
##                   name      set1
##            <character> <logical>
## CL:0000136    fat cell      TRUE
## CL:0000057  fibroblast      TRUE
## 
## 
## $`CL:0011115`
## $`CL:0011115`$name
## [1] "precursor cell"
## 
## $`CL:0011115`$descendents
## DataFrame with 2 rows and 2 columns
##                          name      set1
##                   <character> <logical>
## CL:0000047 neuronal stem cell      TRUE
## CL:0000576           monocyte      TRUE
## 
## 
## $`CL:0000066`
## $`CL:0000066`$name
## [1] "epithelial cell"
## 
## $`CL:0000066`$descendents
## DataFrame with 3 rows and 2 columns
##                        name      set1
##                 <character> <logical>
## CL:0000115 endothelial cell      TRUE
## CL:0000182       hepatocyte      TRUE
## CL:0000065   ependymal cell      TRUE
```

We can also use this function to synchronize multiple sets of terms to the same resolution.
Here, we consider the ImmGen dataset [@ImmGenRef], which provides highly resolved annotation of immune cell types.
The `findCommonAncestors()` function specifies the origins of the descendents for each LCA,
allowing us to focus on LCAs that have representatives in both sets of terms.


```r
ref2 <- ImmGenData(cell.ont="nonna")
lca2 <- findCommonAncestors(MouseRNA=ref$label.ont,
    ImmGen=ref2$label.ont, g=g, names=cl$name)
head(lca2)
```

```
## $`CL:0000126`
## $`CL:0000126`$name
## [1] "macroglial cell"
## 
## $`CL:0000126`$descendents
## DataFrame with 2 rows and 3 columns
##                       name  MouseRNA    ImmGen
##                <character> <logical> <logical>
## CL:0000127       astrocyte      TRUE     FALSE
## CL:0000128 oligodendrocyte      TRUE     FALSE
## 
## 
## $`CL:0000393`
## $`CL:0000393`$name
## [1] "electrically responsive cell"
## 
## $`CL:0000393`$descendents
## DataFrame with 2 rows and 3 columns
##                           name  MouseRNA    ImmGen
##                    <character> <logical> <logical>
## CL:0000540              neuron      TRUE     FALSE
## CL:0000746 cardiac muscle cell      TRUE     FALSE
## 
## 
## $`CL:0000623`
## $`CL:0000623`$name
## [1] "natural killer cell"
## 
## $`CL:0000623`$descendents
## DataFrame with 2 rows and 3 columns
##                              name  MouseRNA    ImmGen
##                       <character> <logical> <logical>
## CL:0000623    natural killer cell      TRUE      TRUE
## CL:0002438 NK1.1-positive natur..     FALSE      TRUE
## 
## 
## $`CL:0000813`
## $`CL:0000813`$name
## [1] "memory T cell"
## 
## $`CL:0000813`$descendents
## DataFrame with 2 rows and 3 columns
##                              name  MouseRNA    ImmGen
##                       <character> <logical> <logical>
## CL:0000897 CD4-positive, alpha-..     FALSE      TRUE
## CL:0000909 CD8-positive, alpha-..     FALSE      TRUE
## 
## 
## $`CL:0000815`
## $`CL:0000815`$name
## [1] "regulatory T cell"
## 
## $`CL:0000815`$descendents
## DataFrame with 2 rows and 3 columns
##                              name  MouseRNA    ImmGen
##                       <character> <logical> <logical>
## CL:0000792 CD4-positive, CD25-p..     FALSE      TRUE
## CL:0000815      regulatory T cell     FALSE      TRUE
## 
## 
## $`CL:0000819`
## $`CL:0000819`$name
## [1] "B-1 B cell"
## 
## $`CL:0000819`$descendents
## DataFrame with 2 rows and 3 columns
##                   name  MouseRNA    ImmGen
##            <character> <logical> <logical>
## CL:0000820 B-1a B cell     FALSE      TRUE
## CL:0000821 B-1b B cell     FALSE      TRUE
```

For example, we might notice that the mouse RNA-seq reference only has a single "T cell" term.
To synchronize resolution across references, 
we would need to roll up all of the ImmGen's finely resolved subsets into that LCA as shown below.
Of course, this results in some loss of precision and information;
whether this is an acceptable price for simpler interpretation is a decision that is left to the user.


```r
children <- lca2$`CL:0000084`$descendents
children
```

```
## DataFrame with 35 rows and 3 columns
##                              name  MouseRNA    ImmGen
##                       <character> <logical> <logical>
## CL:0000084                 T cell      TRUE      TRUE
## CL:0002427 resting double-posit..     FALSE      TRUE
## CL:0000809 double-positive, alp..     FALSE      TRUE
## CL:0002429 CD69-positive double..     FALSE      TRUE
## CL:0000624 CD4-positive, alpha-..     FALSE      TRUE
## ...                           ...       ...       ...
## CL:0002415 immature Vgamma1.1-p..     FALSE      TRUE
## CL:0002411 Vgamma1.1-positive, ..     FALSE      TRUE
## CL:0002416 mature Vgamma1.1-pos..     FALSE      TRUE
## CL:0002407 mature Vgamma2-posit..     FALSE      TRUE
## CL:0000815      regulatory T cell     FALSE      TRUE
```

```r
# Synchronization:
synced.mm <- ref$label.ont
synced.mm[synced.mm %in% rownames(children)] <- "CL:0000084"
synced.ig <- ref2$label.ont
synced.ig[synced.ig %in% rownames(children)] <- "CL:0000084"
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
 [1] igraph_1.2.5                SingleR_1.3.5              
 [3] SummarizedExperiment_1.19.5 DelayedArray_0.15.1        
 [5] matrixStats_0.56.0          Biobase_2.49.0             
 [7] GenomicRanges_1.41.5        GenomeInfoDb_1.25.1        
 [9] IRanges_2.23.9              S4Vectors_0.27.12          
[11] BiocGenerics_0.35.4         ontoProc_1.11.1            
[13] ontologyIndex_2.5           BiocStyle_2.17.0           
[15] rebook_0.99.0              

loaded via a namespace (and not attached):
 [1] httr_1.4.1                    BiocSingular_1.5.0           
 [3] AnnotationHub_2.21.0          bit64_0.9-7                  
 [5] DelayedMatrixStats_1.11.0     ontologyPlot_1.4             
 [7] paintmap_1.0                  shiny_1.4.0.2                
 [9] assertthat_0.2.1              interactiveDisplayBase_1.27.5
[11] BiocManager_1.30.10           BiocFileCache_1.13.0         
[13] blob_1.2.1                    GenomeInfoDbData_1.2.3       
[15] yaml_2.2.1                    BiocVersion_3.12.0           
[17] pillar_1.4.4                  RSQLite_2.2.0                
[19] lattice_0.20-41               glue_1.4.1                   
[21] digest_0.6.25                 promises_1.1.1               
[23] XVector_0.29.2                htmltools_0.4.0              
[25] httpuv_1.5.4                  Matrix_1.2-18                
[27] XML_3.99-0.3                  pkgconfig_2.0.3              
[29] bookdown_0.19                 zlibbioc_1.35.0              
[31] purrr_0.3.4                   xtable_1.8-4                 
[33] processx_3.4.2                later_1.1.0.1                
[35] BiocParallel_1.23.0           tibble_3.0.1                 
[37] generics_0.0.2                ellipsis_0.3.1               
[39] DT_0.13                       magrittr_1.5                 
[41] crayon_1.3.4                  CodeDepends_0.6.5            
[43] mime_0.9                      memoise_1.1.0                
[45] evaluate_0.14                 ps_1.3.3                     
[47] graph_1.67.1                  tools_4.0.0                  
[49] lifecycle_0.2.0               stringr_1.4.0                
[51] irlba_2.3.3                   AnnotationDbi_1.51.0         
[53] callr_3.4.3                   compiler_4.0.0               
[55] rsvd_1.0.3                    rlang_0.4.6                  
[57] grid_4.0.0                    RCurl_1.98-1.2               
[59] BiocNeighbors_1.7.0           rappdirs_0.3.1               
[61] htmlwidgets_1.5.1             bitops_1.0-6                 
[63] rmarkdown_2.2                 ExperimentHub_1.15.0         
[65] codetools_0.2-16              DBI_1.1.0                    
[67] curl_4.3                      R6_2.4.1                     
[69] knitr_1.28                    dplyr_1.0.0                  
[71] fastmap_1.0.1                 bit_1.1-15.2                 
[73] Rgraphviz_2.33.0              stringi_1.4.6                
[75] Rcpp_1.0.4.6                  vctrs_0.3.1                  
[77] dbplyr_1.4.4                  tidyselect_1.1.0             
[79] xfun_0.14                    
```
</div>
