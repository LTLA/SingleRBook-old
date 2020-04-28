

# Using multiple references

In some cases, we may wish to use multiple references for annotation of a test dataset.
This yield a more comprehensive set of cell types that are not covered by any individual reference, especially when differences in resolution are also considered.
Use of multiple references is supported by simply passing multiple objects to the `ref=` and `label=` argument in `SingleR()`.
We demonstrate below by including another reference (from Blueprint-Encode) in our annotation of the @lamanno2016molecular dataset:


```r
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')

library(scater)
hESCs <- logNormCounts(hESCs)

library(SingleR)
bp.se <- BlueprintEncodeData()
hpca.se <- HumanPrimaryCellAtlasData()

pred.combined <- SingleR(test = hESCs, 
    ref = list(BP=bp.se, HPCA=hpca.se), 
    labels = list(bp.se$label.main, hpca.se$label.main))
```

The output is the same form as previously described, and we can easily gain access to the combined set of labels:


```r
table(pred.combined$labels)
```

```
## 
##            Astrocyte           Astrocytes         Chondrocytes 
##                   22                    2                    1 
## Embryonic_stem_cells         Erythrocytes                  HSC 
##                  128                    1                    1 
##            iPS_cells      Mesangial cells Neuroepithelial_cell 
##                  193                    1                 1016 
##              Neurons  Smooth_muscle_cells    Tissue_stem_cells 
##                  342                    7                    1
```

Our strategy is to perform annotation on each reference separately and then take the highest-scoring label across references.
This provides a light-weight approach to combining information from multiple references while avoiding batch effects and the need for up-front harmonization.
(Of course, the main practical difficulty of this approach is that the same cell type may have different labels across references, which will require some implicit harmonization during interpretation.)
Further comments on the justification behind the choice of this method can be found at `?"combine-predictions"`.


