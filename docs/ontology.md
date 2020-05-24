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

As previously discussed in Section \@ref(using-harmonized-labels),
*[SingleR](https://bioconductor.org/packages/3.12/SingleR)* maps the labels in its references to the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl).
The most obvious advantage of doing this is to provide a standardized vocabulary with which to describe cell types,
thus facilitating integrated analyses with multiple references.
However, another useful feature of the Cell Ontology is its hierarchical organization of terms,
allowing us to adjust cell type annotations to the desired resolution.
This represents a more dynamic alternative to the static `label.main` and `label.fine` options in each reference.

## Session information {-}


```r
prettySessionInfo()
```

```
## <button class="aaron-collapse">View session info</button>
## <div class="aaron-content">
## ```
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] BiocStyle_2.17.0 rebook_0.99.0   
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4.6        bookdown_0.19       codetools_0.2-16   
##  [4] XML_3.99-0.3        ps_1.3.3            digest_0.6.25      
##  [7] stats4_4.0.0        magrittr_1.5        evaluate_0.14      
## [10] graph_1.67.0        rlang_0.4.6         stringi_1.4.6      
## [13] callr_3.4.3         rmarkdown_2.1       tools_4.0.0        
## [16] stringr_1.4.0       processx_3.4.2      parallel_4.0.0     
## [19] xfun_0.13           yaml_2.2.1          compiler_4.0.0     
## [22] BiocGenerics_0.35.2 BiocManager_1.30.10 htmltools_0.4.0    
## [25] CodeDepends_0.6.5   knitr_1.28         
## ```
## </div>
```
