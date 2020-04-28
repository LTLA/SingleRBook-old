# Using the built-in references



## Running the algorithm

*[SingleR](https://bioconductor.org/packages/3.11/SingleR)* provides several reference datasets (mostly derived from bulk RNA-seq or microarray data) through dedicated data retrieval functions.
For example, we obtain reference data from the Human Primary Cell Atlas using the `HumanPrimaryCellAtlasData()` function,
which returns a `SummarizedExperiment` object containing matrix of log-expression values with sample-level labels.


```r
library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
```

```
## class: SummarizedExperiment 
## dim: 19363 713 
## metadata(0):
## assays(1): logcounts
## rownames(19363): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
## rowData names(0):
## colnames(713): GSM112490 GSM112491 ... GSM92233 GSM92234
## colData names(3): label.main label.fine label.ont
```

Our test dataset will is taken from @lamanno2016molecular.  
For the sake of speed, we will only label the first 100 cells from this dataset.


```r
library(scRNAseq)
hESCs <- LaMannoBrainData('human-es')

# SingleR() expects log-counts, but the function will also happily take raw
# counts for the test dataset. The reference, however, must have log-values.
library(scater)
hESCs <- logNormCounts(hESCs)
```

We use our `hpca.se` reference to annotate each cell in `hESCs` via the `SingleR()` function, which uses the algorithm described above.
Note that the default marker detection method is to take the genes with the largest positive log-fold changes in the per-label medians for each gene.


```r
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, labels = hpca.se$label.main)
pred.hesc
```

```
## DataFrame with 100 rows and 5 columns
##                                         scores         first.labels
##                                       <matrix>          <character>
## 1772122_301_C02 0.347652:0.139036:0.109547:... Neuroepithelial_cell
## 1772122_180_E05 0.361187:0.155395:0.134934:... Neuroepithelial_cell
## 1772122_300_H02 0.446411:0.218052:0.190084:... Neuroepithelial_cell
## 1772122_180_B09 0.373512:0.172438:0.143537:... Neuroepithelial_cell
## 1772122_180_G04 0.357341:0.157275:0.126511:... Neuroepithelial_cell
## ...                                        ...                  ...
## 1772122_299_E07 0.371989:0.202363:0.169379:... Neuroepithelial_cell
## 1772122_180_D02 0.353314:0.146049:0.115864:... Neuroepithelial_cell
## 1772122_300_D09 0.348789:0.129193:0.136732:... Neuroepithelial_cell
## 1772122_298_F09 0.332361:0.173357:0.141439:... Neuroepithelial_cell
## 1772122_302_A11 0.324928:0.127518:0.101609:... Neuroepithelial_cell
##                       tuning.scores               labels        pruned.labels
##                         <DataFrame>          <character>          <character>
## 1772122_301_C02 0.1824402:0.0991116 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_180_E05 0.1375484:0.0647134              Neurons              Neurons
## 1772122_300_H02 0.2757982:0.1369690 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_180_B09 0.0851623:0.0819878 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_180_G04 0.1988415:0.1016622 Neuroepithelial_cell Neuroepithelial_cell
## ...                             ...                  ...                  ...
## 1772122_299_E07 0.1760025:0.0922504 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_180_D02 0.1967609:0.1124805 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_300_D09 0.0816424:0.0221368 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_298_F09 0.1872499:0.0671893 Neuroepithelial_cell Neuroepithelial_cell
## 1772122_302_A11 0.1560800:0.1051322            Astrocyte            Astrocyte
```

Each row of the output `DataFrame` contains prediction results for a single cell.
Labels are shown before fine-tuning (`first.labels`), after fine-tuning (`labels`) and after pruning (`pruned.labels`), along with the associated scores.
We summarize the distribution of labels across our subset of cells:


```r
table(pred.hesc$labels)
```

```
## 
##            Astrocyte Neuroepithelial_cell              Neurons 
##                   14                   81                    5
```

At this point, it is worth noting that *[SingleR](https://bioconductor.org/packages/3.11/SingleR)* is workflow/package agnostic.
The above example uses `SummarizedExperiment` objects, but the same functions will accept any (log-)normalized expression matrix.

## Available references

The [legacy SingleR package](https://github.com/dviraran/SingleR/tree/master/data) provides RDA files that contain normalized expression values and cell types labels based on bulk RNA-seq, microarray and single-cell RNA-seq data from:

* Blueprint [@blueprintRef] and Encode [@encodeRef],
* the Human Primary Cell Atlas [@hpcaRef],
* the murine [ImmGen](http://www.immgen.org/) [@ImmGenRef], and
* a collection of mouse data sets downloaded from GEO [@Benayoun2019].

The bulk RNA-seq and microarray data sets of the first three reference data sets were obtained from pre-sorted cell populations, i.e., the cell labels of these samples were mostly derived based on the respective sorting/purification strategy, not via *in silico* prediction methods.

Three additional reference datasets from bulk RNA-seq and microarray data for immune cells have also been prepared.
Each of these datasets were also obtained from pre-sorted cell populations:

* The [Database for Immune Cell Expression(/eQTLs/Epigenomics)](https://dice-database.org) [@diceRef],
* Novershtern Hematopoietic Cell Data - [GSE24759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24759) - formerly known as Differentiation Map [@dmapRef], and
* Monaco Immune Cell Data - [GSE107011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011)  [@monaco_immuneRef].

The characteristics of each dataset are summarized below:

| Retrieval function |  Organism  | Samples | Sample types |  No. of main labels  | No. of fine labels | Cell type focus |
|------------------|----------|----------|-------------|----------------------|------------|----------|
|`HumanPrimaryCellAtlasData()`| human | 713 | microarrays of sorted cell populations  | 37 |  157 | Non-specific |
|`BlueprintEncodeData()` |  human | 259 | RNA-seq | 24 | 43 | Non-specific |
|`DatabaseImmuneCellExpressionData()` | human | 1561 | RNA-seq | 5 | 15 | Immune |
|`NovershternHematopoieticData()` | human | 211 | microarrays of sorted cell populations | 17 | 38 | Hematopoietic & Immune |
|`MonacoImmuneData()` | human | 114 | RNA-seq | 11 | 29 | Immune |
|`ImmGenData()`|  mouse | 830  | microarrays of sorted cell populations | 20 | 253 | Hematopoietic & Immune |
|`MouseRNAseqData()`| mouse |358  |RNA-seq| 18  | 28 | Non-specific |

Details for each dataset can be viewed on the corresponding help page for its retrieval function (e.g., `?ImmGenData`).
The available sample types in each set can be viewed in the collapsible sections below.
The cell types in each dataset have also been manually mapped to the [Cell Ontology](https://www.ebi.ac.uk/ols/ontologies/cl), which provides a standardized vocabulary for comparison of labels across studies.

<details>
  <summary>`BlueprintEncodeData` Labels</summary>


|                                                               |label.main        |label.fine                    |label.ont  |
|:--------------------------------------------------------------|:-----------------|:-----------------------------|:----------|
|mature.neutrophil                                              |Neutrophils       |Neutrophils                   |CL:0000775 |
|CD14.positive..CD16.negative.classical.monocyte                |Monocytes         |Monocytes                     |CL:0000576 |
|megakaryocyte.erythroid.progenitor.cell                        |HSC               |MEP                           |CL:0000050 |
|CD4.positive..alpha.beta.T.cell                                |CD4+ T-cells      |CD4+ T-cells                  |CL:0000624 |
|regulatory.T.cell                                              |CD4+ T-cells      |Tregs                         |CL:0000815 |
|central.memory.CD4.positive..alpha.beta.T.cell                 |CD4+ T-cells      |CD4+ Tcm                      |CL:0000904 |
|effector.memory.CD4.positive..alpha.beta.T.cell                |CD4+ T-cells      |CD4+ Tem                      |CL:0000905 |
|central.memory.CD8.positive..alpha.beta.T.cell                 |CD8+ T-cells      |CD8+ Tcm                      |CL:0000907 |
|effector.memory.CD8.positive..alpha.beta.T.cell                |CD8+ T-cells      |CD8+ Tem                      |CL:0000913 |
|cytotoxic.CD56.dim.natural.killer.cell                         |NK cells          |NK cells                      |CL:0000623 |
|CD38.negative.naive.B.cell                                     |B-cells           |naive B-cells                 |CL:0000788 |
|memory.B.cell                                                  |B-cells           |Memory B-cells                |CL:0000787 |
|class.switched.memory.B.cell                                   |B-cells           |Class-switched memory B-cells |CL:0000972 |
|hematopoietic.stem.cell                                        |HSC               |HSC                           |CL:0000037 |
|hematopoietic.multipotent.progenitor.cell                      |HSC               |MPP                           |CL:0000837 |
|common.lymphoid.progenitor                                     |HSC               |CLP                           |CL:0000051 |
|granulocyte.monocyte.progenitor.cell                           |HSC               |GMP                           |CL:0000557 |
|macrophage                                                     |Macrophages       |Macrophages                   |CL:0000235 |
|CD8.positive..alpha.beta.T.cell                                |CD8+ T-cells      |CD8+ T-cells                  |CL:0000625 |
|erythroblast                                                   |Erythrocytes      |Erythrocytes                  |CL:0000232 |
|CD34.negative..CD41.positive..CD42.positive.megakaryocyte.cell |HSC               |Megakaryocytes                |CL:0000556 |
|common.myeloid.progenitor                                      |HSC               |CMP                           |CL:0000049 |
|inflammatory.macrophage                                        |Macrophages       |Macrophages M1                |CL:0000863 |
|alternatively.activated.macrophage                             |Macrophages       |Macrophages M2                |CL:0000890 |
|endothelial.cell.of.umbilical.vein..proliferating.             |Endothelial cells |Endothelial cells             |CL:0000115 |
|conventional.dendritic.cell                                    |DC                |DC                            |CL:0000451 |
|mature.eosinophil                                              |Eosinophils       |Eosinophils                   |CL:0000771 |
|plasma.cell                                                    |B-cells           |Plasma cells                  |CL:0000786 |
|articular.chondrocyte.of.knee.joint                            |Chondrocytes      |Chondrocytes                  |CL:0000138 |
|pericardium.fibroblast                                         |Fibroblasts       |Fibroblasts                   |CL:0000057 |
|smooth.muscle.cell.of.the.umbilical.artery                     |Smooth muscle     |Smooth muscle                 |CL:0000192 |
|epithelial.cell.of.proximal.tubule                             |Epithelial cells  |Epithelial cells              |CL:0000066 |
|melanocyte.of.skin                                             |Melanocytes       |Melanocytes                   |CL:0000148 |
|skeletal.muscle.tissue                                         |Skeletal muscle   |Skeletal muscle               |CL:0000188 |
|hair.follicular.keratinocyte                                   |Keratinocytes     |Keratinocytes                 |CL:0000312 |
|lung.microvascular.endothelial.cell                            |Endothelial cells |mv Endothelial cells          |CL:2000008 |
|regular.cardiac.myocyte                                        |Myocytes          |Myocytes                      |CL:0000187 |
|adipose.tissue.of.omentum                                      |Adipocytes        |Adipocytes                    |CL:0000136 |
|Purkinje.cell                                                  |Neurons           |Neurons                       |CL:0000540 |
|pericyte.cell                                                  |Pericytes         |Pericytes                     |CL:0000669 |
|subcutaneous.preadipocyte                                      |Adipocytes        |Preadipocytes                 |NA         |
|astrocyte                                                      |Astrocytes        |Astrocytes                    |CL:0000127 |
|mesangial.cell                                                 |Mesangial cells   |Mesangial cells               |CL:0000650 |
</details>

<details>
  <summary>`HumanPrimaryCellAtlasData` Labels</summary>


|                                   |label.main           |label.fine                                             |label.ont  |
|:----------------------------------|:--------------------|:------------------------------------------------------|:----------|
|GSM112490                          |DC                   |DC:monocyte-derived:immature                           |CL:0000840 |
|GSM112541                          |DC                   |DC:monocyte-derived:Galectin-1                         |CL:0000451 |
|GSM112665                          |DC                   |DC:monocyte-derived:LPS                                |CL:0000451 |
|GSM112668                          |DC                   |DC:monocyte-derived                                    |CL:0000451 |
|GSM116101                          |Smooth_muscle_cells  |Smooth_muscle_cells:bronchial:vit_D                    |CL:0002598 |
|GSM116104                          |Smooth_muscle_cells  |Smooth_muscle_cells:bronchial                          |CL:0002598 |
|GSM119354                          |Epithelial_cells     |Epithelial_cells:bronchial                             |CL:0002328 |
|GSM1209554_HH1763_UI33plus2_201004 |B_cell               |B_cell                                                 |CL:0000236 |
|GSM1209558_HH1713_u133plus2_011004 |Neutrophils          |Neutrophil                                             |CL:0000775 |
|GSM1209561_TW1681_u133plus2_061004 |T_cells              |T_cell:CD8+_Central_memory                             |CL:0000907 |
|GSM1209564_HH1765_UI33plus2_201004 |T_cells              |T_cell:CD8+                                            |CL:0000625 |
|GSM1209565_HH1769_UI33plus2_201004 |T_cells              |T_cell:CD4+                                            |CL:0000624 |
|GSM1209573_TW1678_u133plus2_061004 |T_cells              |T_cell:CD8+_effector_memory_RA                         |CL:0001062 |
|GSM1209577_TW1675_u133plus2_061004 |T_cells              |T_cell:CD8+_effector_memory                            |CL:0000913 |
|GSM1209581_TW1676_u133plus2_061004 |T_cells              |T_cell:CD8+_naive                                      |CL:0000900 |
|GSM1209585_HH1762_UI33plus2_201004 |Monocyte             |Monocyte                                               |CL:0000576 |
|GSM1209591_HH1719_u133plus2_011004 |Erythroblast         |Erythroblast                                           |CL:0000765 |
|GSM1209599_HH1715_u133plus2_011004 |BM & Prog.           |BM                                                     |NA         |
|GSM132921                          |DC                   |DC:monocyte-derived:rosiglitazone                      |CL:0000451 |
|GSM132922                          |DC                   |DC:monocyte-derived:AM580                              |CL:0000451 |
|GSM132926                          |DC                   |DC:monocyte-derived:rosiglitazone/AGN193109            |CL:0000451 |
|GSM140970                          |DC                   |DC:monocyte-derived:anti-DC-SIGN_2h                    |CL:0000451 |
|GSM141251                          |Endothelial_cells    |Endothelial_cells:HUVEC                                |CL:0002618 |
|GSM141252                          |Endothelial_cells    |Endothelial_cells:HUVEC:Borrelia_burgdorferi           |CL:0002618 |
|GSM141255                          |Endothelial_cells    |Endothelial_cells:HUVEC:IFNg                           |CL:0002618 |
|GSM143717                          |Endothelial_cells    |Endothelial_cells:lymphatic                            |CL:0002138 |
|GSM143728                          |Endothelial_cells    |Endothelial_cells:HUVEC:Serum_Amyloid_A                |CL:0002618 |
|GSM143907                          |Endothelial_cells    |Endothelial_cells:lymphatic:TNFa_48h                   |CL:0002138 |
|GSM153893                          |T_cells              |T_cell:effector                                        |CL:0000911 |
|GSM154081                          |T_cells              |T_cell:CCR10+CLA+1,25(OH)2_vit_D3/IL-12                |CL:0000084 |
|GSM154084                          |T_cells              |T_cell:CCR10-CLA+1,25(OH)2_vit_D3/IL-12                |CL:0000084 |
|GSM158468                          |Gametocytes          |Gametocytes:spermatocyte                               |CL:0000017 |
|GSM160532                          |DC                   |DC:monocyte-derived:A._fumigatus_germ_tubes_6h         |CL:0000451 |
|GSM172865                          |Neurons              |Neurons:ES_cell-derived_neural_precursor               |CL:0000031 |
|GSM173532                          |Keratinocytes        |Keratinocytes                                          |CL:0000312 |
|GSM173535                          |Keratinocytes        |Keratinocytes:IL19                                     |CL:0000312 |
|GSM173538                          |Keratinocytes        |Keratinocytes:IL20                                     |CL:0000312 |
|GSM173541                          |Keratinocytes        |Keratinocytes:IL22                                     |CL:0000312 |
|GSM173544                          |Keratinocytes        |Keratinocytes:IL24                                     |CL:0000312 |
|GSM173547                          |Keratinocytes        |Keratinocytes:IL26                                     |CL:0000312 |
|GSM173550                          |Keratinocytes        |Keratinocytes:KGF                                      |CL:0000312 |
|GSM173553                          |Keratinocytes        |Keratinocytes:IFNg                                     |CL:0000312 |
|GSM173555                          |Keratinocytes        |Keratinocytes:IL1b                                     |CL:0000312 |
|GSM178549                          |HSC_-G-CSF           |HSC_-G-CSF                                             |CL:0000037 |
|GSM181971                          |DC                   |DC:monocyte-derived:mature                             |CL:0000841 |
|GSM182001                          |Monocyte             |Monocyte:anti-FcgRIIB                                  |CL:0000576 |
|GSM183165                          |Macrophage           |Macrophage:monocyte-derived:IL-4/cntrl                 |CL:0000235 |
|GSM183217                          |Macrophage           |Macrophage:monocyte-derived:IL-4/Dex/cntrl             |CL:0000235 |
|GSM183392                          |Macrophage           |Macrophage:monocyte-derived:IL-4/Dex/TGFb              |CL:0000235 |
|GSM183483                          |Macrophage           |Macrophage:monocyte-derived:IL-4/TGFb                  |CL:0000235 |
|GSM189451                          |Monocyte             |Monocyte:leukotriene_D4                                |CL:0000576 |
|GSM198942                          |NK_cell              |NK_cell                                                |CL:0000623 |
|GSM198943                          |NK_cell              |NK_cell:IL2                                            |CL:0000623 |
|GSM225042                          |Embryonic_stem_cells |Embryonic_stem_cells                                   |CL:0002322 |
|GSM239260                          |Tissue_stem_cells    |Tissue_stem_cells:iliac_MSC                            |CL:0000134 |
|GSM239606                          |Chondrocytes         |Chondrocytes:MSC-derived                               |CL:0000138 |
|GSM239616                          |Osteoblasts          |Osteoblasts                                            |CL:0000062 |
|GSM250019                          |Tissue_stem_cells    |Tissue_stem_cells:BM_MSC                               |CL:0000134 |
|GSM260308                          |Osteoblasts          |Osteoblasts:BMP2                                       |CL:0000062 |
|GSM260663                          |Tissue_stem_cells    |Tissue_stem_cells:BM_MSC:BMP2                          |CL:0000134 |
|GSM260675                          |Tissue_stem_cells    |Tissue_stem_cells:BM_MSC:TGFb3                         |CL:0000134 |
|GSM260693                          |DC                   |DC:monocyte-derived:Poly(IC)                           |CL:0000451 |
|GSM260696                          |DC                   |DC:monocyte-derived:CD40L                              |CL:0000451 |
|GSM260699                          |DC                   |DC:monocyte-derived:Schuler_treatment                  |CL:0000451 |
|GSM264757                          |DC                   |DC:monocyte-derived:antiCD40/VAF347                    |CL:0000451 |
|GSM265494                          |Tissue_stem_cells    |Tissue_stem_cells:dental_pulp                          |CL:0002148 |
|GSM279572                          |T_cells              |T_cell:CD4+_central_memory                             |CL:0000904 |
|GSM279577                          |T_cells              |T_cell:CD4+_effector_memory                            |CL:0000905 |
|GSM279581                          |T_cells              |T_cell:CD4+_Naive                                      |CL:0000895 |
|GSM287216                          |Smooth_muscle_cells  |Smooth_muscle_cells:vascular                           |CL:0000359 |
|GSM287217                          |Smooth_muscle_cells  |Smooth_muscle_cells:vascular:IL-17                     |CL:0000359 |
|GSM289612                          |BM                   |BM                                                     |NA         |
|GSM290414                          |Platelets            |Platelets                                              |CL:0000233 |
|GSM299095                          |Epithelial_cells     |Epithelial_cells:bladder                               |CL:0000066 |
|GSM299556                          |Macrophage           |Macrophage:monocyte-derived                            |CL:0000235 |
|GSM299557                          |Macrophage           |Macrophage:monocyte-derived:M-CSF                      |CL:0000235 |
|GSM299558                          |Macrophage           |Macrophage:monocyte-derived:M-CSF/IFNg                 |CL:0000235 |
|GSM299559                          |Macrophage           |Macrophage:monocyte-derived:M-CSF/Pam3Cys              |CL:0000235 |
|GSM299560                          |Macrophage           |Macrophage:monocyte-derived:M-CSF/IFNg/Pam3Cys         |CL:0000235 |
|GSM300389                          |Macrophage           |Macrophage:monocyte-derived:IFNa                       |CL:0000235 |
|GSM304260                          |Gametocytes          |Gametocytes:oocyte                                     |CL:0000023 |
|GSM305433                          |Monocyte             |Monocyte:F._tularensis_novicida                        |CL:0000576 |
|GSM305786                          |Endothelial_cells    |Endothelial_cells:HUVEC:B._anthracis_LT                |CL:0002618 |
|GSM310429                          |B_cell               |B_cell:Germinal_center                                 |CL:0000844 |
|GSM310432                          |B_cell               |B_cell:Plasma_cell                                     |CL:0000786 |
|GSM310435                          |B_cell               |B_cell:Naive                                           |CL:0000788 |
|GSM310438                          |B_cell               |B_cell:Memory                                          |CL:0000787 |
|GSM320544                          |DC                   |DC:monocyte-derived:AEC-conditioned                    |CL:0000451 |
|GSM322374                          |Tissue_stem_cells    |Tissue_stem_cells:lipoma-derived_MSC                   |CL:0000134 |
|GSM322376                          |Tissue_stem_cells    |Tissue_stem_cells:adipose-derived_MSC_AM3              |CL:0000134 |
|GSM330314                          |Endothelial_cells    |Endothelial_cells:HUVEC:FPV-infected                   |CL:0002618 |
|GSM330315                          |Endothelial_cells    |Endothelial_cells:HUVEC:PR8-infected                   |CL:0002618 |
|GSM330316                          |Endothelial_cells    |Endothelial_cells:HUVEC:H5N1-infected                  |CL:0002618 |
|GSM343803                          |Macrophage           |Macrophage:monocyte-derived:S._aureus                  |CL:0000235 |
|GSM346941                          |Fibroblasts          |Fibroblasts:foreskin                                   |CL:1001608 |
|GSM347916                          |iPS_cells            |iPS_cells:skin_fibroblast-derived                      |NA         |
|GSM347919                          |iPS_cells            |iPS_cells:skin_fibroblast                              |NA         |
|GSM349848                          |T_cells              |T_cell:gamma-delta                                     |CL:0000798 |
|GSM350084                          |Monocyte             |Monocyte:CD14+                                         |CL:0001054 |
|GSM359332                          |Macrophage           |Macrophage:Alveolar                                    |CL:0000583 |
|GSM359758                          |Macrophage           |Macrophage:Alveolar:B._anthacis_spores                 |CL:0000583 |
|GSM361272                          |Neutrophils          |Neutrophil:inflam                                      |CL:0000775 |
|GSM366942                          |iPS_cells            |iPS_cells:PDB_fibroblasts                              |NA         |
|GSM367219                          |iPS_cells            |iPS_cells:PDB_1lox-17Puro-5                            |NA         |
|GSM367240                          |iPS_cells            |iPS_cells:PDB_1lox-17Puro-10                           |NA         |
|GSM367241                          |iPS_cells            |iPS_cells:PDB_1lox-21Puro-20                           |NA         |
|GSM367242                          |iPS_cells            |iPS_cells:PDB_1lox-21Puro-26                           |NA         |
|GSM367243                          |iPS_cells            |iPS_cells:PDB_2lox-5                                   |NA         |
|GSM367244                          |iPS_cells            |iPS_cells:PDB_2lox-22                                  |NA         |
|GSM367245                          |iPS_cells            |iPS_cells:PDB_2lox-21                                  |NA         |
|GSM367258                          |iPS_cells            |iPS_cells:PDB_2lox-17                                  |NA         |
|GSM372142                          |iPS_cells            |iPS_cells:CRL2097_foreskin                             |NA         |
|GSM372154                          |iPS_cells            |iPS_cells:CRL2097_foreskin-derived:d20_hepatic_diff    |NA         |
|GSM372157                          |iPS_cells            |iPS_cells:CRL2097_foreskin-derived:undiff.             |NA         |
|GSM381339                          |B_cell               |B_cell:CXCR4+_centroblast                              |CL:0000965 |
|GSM381340                          |B_cell               |B_cell:CXCR4-_centrocyte                               |CL:0000966 |
|GSM385338                          |Endothelial_cells    |Endothelial_cells:HUVEC:VEGF                           |CL:0002618 |
|GSM402707                          |iPS_cells            |iPS_cells:fibroblasts                                  |NA         |
|GSM402717                          |iPS_cells            |iPS_cells:fibroblast-derived:Direct_del._reprog        |NA         |
|GSM402806                          |iPS_cells            |iPS_cells:fibroblast-derived:Retroviral_transf         |NA         |
|GSM410672                          |Endothelial_cells    |Endothelial_cells:lymphatic:KSHV                       |CL:0002138 |
|GSM410678                          |Endothelial_cells    |Endothelial_cells:blood_vessel                         |CL:0000071 |
|GSM422109                          |Monocyte             |Monocyte:CD16-                                         |CL:0000576 |
|GSM422113                          |Monocyte             |Monocyte:CD16+                                         |CL:0000576 |
|GSM451153                          |Tissue_stem_cells    |Tissue_stem_cells:BM_MSC:osteogenic                    |CL:0000134 |
|GSM456349                          |Hepatocytes          |Hepatocytes                                            |CL:0000182 |
|GSM466515                          |Neutrophils          |Neutrophil:uropathogenic_E._coli_UTI89                 |CL:0000775 |
|GSM466516                          |Neutrophils          |Neutrophil:commensal_E._coli_MG1655                    |CL:0000775 |
|GSM469125                          |MSC                  |MSC                                                    |CL:0000134 |
|GSM469409                          |Neuroepithelial_cell |Neuroepithelial_cell:ESC-derived                       |CL:0002259 |
|GSM469411                          |Astrocyte            |Astrocyte:Embryonic_stem_cell-derived                  |CL:0000127 |
|GSM476783                          |Endothelial_cells    |Endothelial_cells:HUVEC:IL-1b                          |CL:0002618 |
|GSM483480                          |HSC_CD34+            |HSC_CD34+                                              |CL:0000037 |
|GSM488968                          |CMP                  |CMP                                                    |CL:0000049 |
|GSM488970                          |GMP                  |GMP                                                    |CL:0000557 |
|GSM488972                          |B_cell               |B_cell:immature                                        |CL:0000816 |
|GSM488974                          |MEP                  |MEP                                                    |CL:0000050 |
|GSM488976                          |Myelocyte            |Myelocyte                                              |CL:0002193 |
|GSM488978                          |Pre-B_cell_CD34-     |Pre-B_cell_CD34-                                       |CL:0000955 |
|GSM488980                          |Pro-B_cell_CD34+     |Pro-B_cell_CD34+                                       |CL:0002048 |
|GSM488982                          |Pro-Myelocyte        |Pro-Myelocyte                                          |CL:0000836 |
|GSM492834                          |Smooth_muscle_cells  |Smooth_muscle_cells:umbilical_vein                     |CL:0002588 |
|GSM500995                          |iPS_cells            |iPS_cells:foreskin_fibrobasts                          |NA         |
|GSM500996                          |iPS_cells            |iPS_cells:iPS:minicircle-derived                       |NA         |
|GSM501001                          |iPS_cells            |iPS_cells:adipose_stem_cells                           |NA         |
|GSM501004                          |iPS_cells            |iPS_cells:adipose_stem_cell-derived:lentiviral         |NA         |
|GSM501007                          |iPS_cells            |iPS_cells:adipose_stem_cell-derived:minicircle-derived |NA         |
|GSM501890                          |Fibroblasts          |Fibroblasts:breast                                     |CL:0002555 |
|GSM514669                          |Monocyte             |Monocyte:MCSF                                          |CL:0000576 |
|GSM514671                          |Monocyte             |Monocyte:CXCL4                                         |CL:0000576 |
|GSM53382                           |Neurons              |Neurons:adrenal_medulla_cell_line                      |CL:0000540 |
|GSM540714                          |Tissue_stem_cells    |Tissue_stem_cells:CD326-CD56+                          |CL:0000222 |
|GSM542578                          |NK_cell              |NK_cell:CD56hiCD62L+                                   |CL:0000623 |
|GSM547998                          |T_cells              |T_cell:Treg:Naive                                      |CL:0002677 |
|GSM549577                          |Neutrophils          |Neutrophil:LPS                                         |CL:0000775 |
|GSM549581                          |Neutrophils          |Neutrophil:GM-CSF_IFNg                                 |CL:0000775 |
|GSM556665                          |Monocyte             |Monocyte:S._typhimurium_flagellin                      |CL:0000576 |
|GSM92231                           |Neurons              |Neurons:Schwann_cell                                   |CL:0002573 |
</details>

<details>
  <summary>`DatabaseImmuneCellExpressionData` Labels</summary>


|         |label.main    |label.fine                       |label.ont  |
|:--------|:-------------|:--------------------------------|:----------|
|TPM_1    |B cells       |B cells, naive                   |CL:0000788 |
|TPM_1.1  |Monocytes     |Monocytes, CD14+                 |CL:0002057 |
|TPM_1.2  |Monocytes     |Monocytes, CD16+                 |CL:0002396 |
|TPM_1.3  |NK cells      |NK cells                         |CL:0000623 |
|TPM_1.4  |T cells, CD4+ |T cells, CD4+, memory TREG       |CL:0000792 |
|TPM_1.5  |T cells, CD4+ |T cells, CD4+, naive             |CL:0000895 |
|TPM_1.6  |T cells, CD4+ |T cells, CD4+, naive, stimulated |CL:0000896 |
|TPM_1.7  |T cells, CD4+ |T cells, CD4+, naive TREG        |CL:0001045 |
|TPM_1.8  |T cells, CD4+ |T cells, CD4+, TFH               |CL:0002038 |
|TPM_1.9  |T cells, CD4+ |T cells, CD4+, Th1               |CL:0000545 |
|TPM_1.10 |T cells, CD4+ |T cells, CD4+, Th1_17            |CL:0000492 |
|TPM_1.11 |T cells, CD4+ |T cells, CD4+, Th17              |CL:0000899 |
|TPM_1.12 |T cells, CD4+ |T cells, CD4+, Th2               |CL:0000546 |
|TPM_1.13 |T cells, CD8+ |T cells, CD8+, naive             |CL:0000900 |
|TPM_1.14 |T cells, CD8+ |T cells, CD8+, naive, stimulated |CL:0000906 |
</details>

<details>
  <summary>`NovershternHematopoieticData` Labels</summary>


|          |label.main      |label.fine                                 |label.ont  |
|:---------|:---------------|:------------------------------------------|:----------|
|GSM609632 |Basophils       |Basophils                                  |CL:0000767 |
|GSM609638 |B cells         |Naive B cells                              |CL:0000788 |
|GSM609643 |B cells         |Mature B cells class able to switch        |CL:0000970 |
|GSM609648 |B cells         |Mature B cells                             |CL:0000785 |
|GSM609653 |B cells         |Mature B cells class switched              |CL:0000972 |
|GSM609658 |CMPs            |Common myeloid progenitors                 |CL:0000049 |
|GSM609662 |Dendritic cells |Plasmacytoid Dendritic Cells               |CL:0000784 |
|GSM609667 |Dendritic cells |Myeloid Dendritic Cells                    |CL:0000782 |
|GSM609672 |Eosinophils     |Eosinophils                                |CL:0000771 |
|GSM609677 |Erythroid cells |Erythroid_CD34+ CD71+ GlyA-                |CL:0002003 |
|GSM609684 |Erythroid cells |Erythroid_CD34- CD71+ GlyA-                |CL:0002004 |
|GSM609691 |Erythroid cells |Erythroid_CD34- CD71+ GlyA+                |CL:0002021 |
|GSM609697 |Erythroid cells |Erythroid_CD34- CD71lo GlyA+               |CL:0002016 |
|GSM609704 |Erythroid cells |Erythroid_CD34- CD71- GlyA+                |CL:0002018 |
|GSM609710 |GMPs            |Granulocyte/monocyte progenitors           |CL:0000557 |
|GSM609714 |Granulocytes    |Colony Forming Unit-Granulocytes           |CL:0000094 |
|GSM609719 |Granulocytes    |Granulocytes (Neutrophilic Metamyelocytes) |CL:0000582 |
|GSM609723 |Granulocytes    |Granulocytes (Neutrophils)                 |CL:0000776 |
|GSM609727 |HSCs            |Hematopoietic stem cells_CD133+ CD34dim    |CL:0000037 |
|GSM609737 |HSCs            |Hematopoietic stem cells_CD38- CD34+       |CL:0001024 |
|GSM609741 |Megakaryocytes  |Colony Forming Unit-Megakaryocytic         |CL:0000556 |
|GSM609746 |Megakaryocytes  |Megakaryocytes                             |CL:0000556 |
|GSM609753 |MEPs            |Megakaryocyte/erythroid progenitors        |CL:0000050 |
|GSM609762 |Monocytes       |Colony Forming Unit-Monocytes              |CL:0000576 |
|GSM609766 |Monocytes       |Monocytes                                  |CL:0000576 |
|GSM609771 |NK cells        |Mature NK cells_CD56- CD16+ CD3-           |CL:0000623 |
|GSM609775 |NK cells        |Mature NK cells_CD56+ CD16+ CD3-           |CL:0000623 |
|GSM609780 |NK cells        |Mature NK cells_CD56- CD16- CD3-           |CL:0000623 |
|GSM609785 |NK T cells      |NK T cells                                 |CL:0000814 |
|GSM609789 |B cells         |Early B cells                              |CL:0002046 |
|GSM609793 |B cells         |Pro B cells                                |CL:0000826 |
|GSM609798 |CD8+ T cells    |CD8+ Effector Memory RA                    |CL:0001062 |
|GSM609802 |CD8+ T cells    |Naive CD8+ T cells                         |CL:0000900 |
|GSM609809 |CD8+ T cells    |CD8+ Effector Memory                       |CL:0000913 |
|GSM609815 |CD8+ T cells    |CD8+ Central Memory                        |CL:0000907 |
|GSM609822 |CD4+ T cells    |Naive CD4+ T cells                         |CL:0000895 |
|GSM609829 |CD4+ T cells    |CD4+ Effector Memory                       |CL:0000905 |
|GSM609836 |CD4+ T cells    |CD4+ Central Memory                        |CL:0000904 |
</details>

<details>
  <summary>`MonacoImmuneData` Labels</summary>


|                  |label.main      |label.fine                    |label.ont  |
|:-----------------|:---------------|:-----------------------------|:----------|
|DZQV_CD8_naive    |CD8+ T cells    |Naive CD8 T cells             |CL:0000900 |
|DZQV_CD8_CM       |CD8+ T cells    |Central memory CD8 T cells    |CL:0000907 |
|DZQV_CD8_EM       |CD8+ T cells    |Effector memory CD8 T cells   |CL:0000913 |
|DZQV_CD8_TE       |CD8+ T cells    |Terminal effector CD8 T cells |CL:0001062 |
|DZQV_MAIT         |T cells         |MAIT cells                    |CL:0000940 |
|DZQV_VD2+         |T cells         |Vd2 gd T cells                |CL:0000798 |
|DZQV_VD2-         |T cells         |Non-Vd2 gd T cells            |CL:0000798 |
|DZQV_TFH          |CD4+ T cells    |Follicular helper T cells     |CL:0002038 |
|DZQV_Treg         |CD4+ T cells    |T regulatory cells            |CL:0000815 |
|DZQV_Th1          |CD4+ T cells    |Th1 cells                     |CL:0000545 |
|DZQV_Th1/Th17     |CD4+ T cells    |Th1/Th17 cells                |CL:0000912 |
|DZQV_Th17         |CD4+ T cells    |Th17 cells                    |CL:0000899 |
|DZQV_Th2          |CD4+ T cells    |Th2 cells                     |CL:0000546 |
|DZQV_CD4_naive    |CD4+ T cells    |Naive CD4 T cells             |CL:0000895 |
|DZQV_Progenitor   |Progenitors     |Progenitor cells              |CL:0002043 |
|DZQV_B_naive      |B cells         |Naive B cells                 |CL:0000788 |
|DZQV_B_NSM        |B cells         |Non-switched memory B cells   |CL:0000970 |
|DZQV_B_Ex         |B cells         |Exhausted B cells             |CL:0000236 |
|DZQV_B_SM         |B cells         |Switched memory B cells       |CL:0000972 |
|DZQV_Plasmablasts |B cells         |Plasmablasts                  |CL:0000980 |
|DZQV_C_mono       |Monocytes       |Classical monocytes           |CL:0000860 |
|DZQV_I_mono       |Monocytes       |Intermediate monocytes        |CL:0002393 |
|DZQV_NC_mono      |Monocytes       |Non classical monocytes       |CL:0000875 |
|DZQV_NK           |NK cells        |Natural killer cells          |CL:0000623 |
|DZQV_pDC          |Dendritic cells |Plasmacytoid dendritic cells  |CL:0000784 |
|DZQV_mDC          |Dendritic cells |Myeloid dendritic cells       |CL:0000782 |
|DZQV_Neutrophils  |Neutrophils     |Low-density neutrophils       |CL:0000096 |
|DZQV_Basophils    |Basophils       |Low-density basophils         |CL:0000043 |
|925L_CD4_TE       |CD4+ T cells    |Terminal effector CD4 T cells |CL:0001044 |
</details>

<details>
  <summary>`ImmGenData` Labels</summary>


|                                                                                 |label.main        |label.fine                             |label.ont  |
|:--------------------------------------------------------------------------------|:-----------------|:--------------------------------------|:----------|
|GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL                  |Macrophages       |Macrophages (MF.11C-11B+)              |CL:0000235 |
|GSM1136122_EA07068_260300_MOGENE-1_0-ST-V1_MF.ALV.LU_1.CEL                       |Macrophages       |Macrophages (MF.ALV)                   |CL:0000583 |
|GSM1136125_EA07068_260307_MOGENE-1_0-ST-V1_MO.6+I-.BL_1.CEL                      |Monocytes         |Monocytes (MO.6+I-)                    |CL:0000576 |
|GSM1136126_EA07068_260303_MOGENE-1_0-ST-V1_MO.6+2+.MLN_1.CEL                     |Monocytes         |Monocytes (MO.6+2+)                    |CL:0000576 |
|GSM1282081_EA07068_147711_MOGENE-1_0-ST-V1_B.MEM.SP_1.CEL                        |B cells           |B cells (B.MEM)                        |CL:0000787 |
|GSM1282083_EA07068_122841_MOGENE-1_0-ST-V1_B1A.SP_3.CEL                          |B cells           |B cells (B1A)                          |CL:0000820 |
|GSM1282084_EA07068_267995_MOGENE-1_0-ST-V1_DC.11B+.AT_2.CEL                      |DC                |DC (DC.11B+)                           |CL:0002465 |
|GSM1282087_EA07068_267991_MOGENE-1_0-ST-V1_DC.11B-.AT_1.CEL                      |DC                |DC (DC.11B-)                           |CL:0000990 |
|GSM1282089_EA07068_210451_MOGENE-1_0-ST-V1_DN.SLN.CFA.D6_1.CEL                   |Stromal cells     |Stromal cells (DN.CFA)                 |CL:0000499 |
|GSM1282091_EA07068_210437_MOGENE-1_0-ST-V1_DN.SLN.V2_1.CEL                       |Stromal cells     |Stromal cells (DN)                     |CL:0000499 |
|GSM1282093_EA07068_267987_MOGENE-1_0-ST-V1_EO.AT_1.CEL                           |Eosinophils       |Eosinophils (EO)                       |CL:0000771 |
|GSM1282097_EA07068_204064_MOGENE-1_0-ST-V1_FRC.CAD11.WT.CEL                      |Fibroblasts       |Fibroblasts (FRC.CAD11.WT)             |CL:0000057 |
|GSM1282098_EA07068_210445_MOGENE-1_0-ST-V1_FRC.SLN.CFA.D6_1.CEL                  |Fibroblasts       |Fibroblasts (FRC.CFA)                  |CL:0000057 |
|GSM1282100_EA07068_210431_MOGENE-1_0-ST-V1_FRC.SLN.V2_1.CEL                      |Fibroblasts       |Fibroblasts (FRC)                      |CL:0000057 |
|GSM1282102_EA07068_256271_MOGENE-1_0-ST-V1_GN.BL_4.CEL                           |Neutrophils       |Neutrophils (GN)                       |CL:0000775 |
|GSM1282106_EA07068_210448_MOGENE-1_0-ST-V1_LEC.SLN.CFA.D6_2.CEL                  |Endothelial cells |Endothelial cells (LEC.CFA)            |CL:0000115 |
|GSM1282107_EA07068_210433_MOGENE-1_0-ST-V1_LEC.SLN.V2_1.CEL                      |Endothelial cells |Endothelial cells (LEC)                |CL:0000115 |
|GSM1282109_EA07068_267983_MOGENE-1_0-ST-V1_MF.AT._1.CEL                          |Macrophages       |Macrophages (MF)                       |CL:0000235 |
|GSM1282112_EA07068_201211_MOGENE-1_0-ST-V1_T.DP.69-.E17.TH_1.CEL                 |T cells           |T cells (T.DP.69-)                     |CL:0002427 |
|GSM1282115_EA07068_208625_MOGENE-1_0-ST-V1_T.DP.TH_1.CEL                         |T cells           |T cells (T.DP)                         |CL:0000809 |
|GSM1282118_EA07068_208628_MOGENE-1_0-ST-V1_T.DP69+.TH_1.CEL                      |T cells           |T cells (T.DP69+)                      |CL:0002429 |
|GSM1308350_EA07068_256264_MOGENE-1_0-ST-V1_MF.F480HI.GATA6KO.PC_1.CEL            |Macrophages       |Macrophages (MF.F480HI.GATA6KO)        |CL:0000235 |
|GSM1308353_EA07068_256261_MOGENE-1_0-ST-V1_MF.F480HI.CTRL.PC_1.CEL               |Macrophages       |Macrophages (MF.F480HI.CTRL)           |CL:0000235 |
|GSM1358373_EA07068_232185_MOGENE-1_0-ST-V1_CD4.1H_1.CEL                          |T cells           |T cells (T.CD4.1H)                     |CL:0000624 |
|GSM1358375_EA07068_232189_MOGENE-1_0-ST-V1_CD4.24H_1.CEL                         |T cells           |T cells (T.CD4.24H)                    |CL:0000624 |
|GSM1358377_EA07068_232191_MOGENE-1_0-ST-V1_CD4.48H_1.CEL                         |T cells           |T cells (T.CD4.48H)                    |CL:0000624 |
|GSM1358379_EA07068_232187_MOGENE-1_0-ST-V1_CD4.5H_1.CEL                          |T cells           |T cells (T.CD4.5H)                     |CL:0000624 |
|GSM1358381_EA07068_232193_MOGENE-1_0-ST-V1_CD4.96H_1.CEL                         |T cells           |T cells (T.CD4.96H)                    |CL:0000624 |
|GSM1358382_EA07068_232183_MOGENE-1_0-ST-V1_CD4.CTR_1.CEL                         |T cells           |T cells (T.CD4.CTR)                    |CL:0000624 |
|GSM1358384_EA07068_232186_MOGENE-1_0-ST-V1_CD8.1H_1.CEL                          |T cells           |T cells (T.CD8.1H)                     |CL:0000625 |
|GSM1358386_EA07068_232190_MOGENE-1_0-ST-V1_CD8.24H_1.CEL                         |T cells           |T cells (T.CD8.24H)                    |CL:0000625 |
|GSM1358388_EA07068_232192_MOGENE-1_0-ST-V1_CD8.48H_1.CEL                         |T cells           |T cells (T.CD8.48H)                    |CL:0000625 |
|GSM1358390_EA07068_232188_MOGENE-1_0-ST-V1_CD8.5H_1.CEL                          |T cells           |T cells (T.CD8.5H)                     |CL:0000625 |
|GSM1358392_EA07068_232194_MOGENE-1_0-ST-V1_CD8.96H_1.CEL                         |T cells           |T cells (T.CD8.96H)                    |CL:0000625 |
|GSM1358393_EA07068_232184_MOGENE-1_0-ST-V1_CD8.CTR_1.CEL                         |T cells           |T cells (T.CD8.CTR)                    |CL:0000625 |
|GSM1398469_EA07068_117717_MOGENE-1_0-ST-V1_MF.PPAR-.LU_2.CEL                     |Macrophages       |Macrophages (MFAR-)                    |CL:0000235 |
|GSM1398483_EA07068_260311_MOGENE-1_0-ST-V1_MO.LU_1.CEL                           |Monocytes         |Monocytes (MO)                         |CL:0000576 |
|GSM1585312_EA07068_339227_MOGENE-1_0-ST-V1_ILC1.CD127+.SP_1.CEL                  |ILC               |ILC (ILC1.CD127+)                      |CL:0001067 |
|GSM1585315_EA07068_339236_MOGENE-1_0-ST-V1_LIV.ILC1.DX5-_1.CEL                   |ILC               |ILC (LIV.ILC1.DX5-)                    |CL:0001067 |
|GSM1585318_EA07068_339248_MOGENE-1_0-ST-V1_LPL.NCR+ILC1_1.CEL                    |ILC               |ILC (LPL.NCR+ILC1)                     |CL:0001067 |
|GSM1585320_EA07068_339234_MOGENE-1_0-ST-V1_ILC2.SI_2.CEL                         |ILC               |ILC (ILC2)                             |CL:0001069 |
|GSM1585322_EA07068_339251_MOGENE-1_0-ST-V1_LPL.NCR+ILC3_1.CEL                    |ILC               |ILC (LPL.NCR+ILC3)                     |CL:0001071 |
|GSM1585325_EA07068_305553_MOGENE-1_0-ST-V1_ILC3.LTI.CD4+.SI_4.CEL                |ILC               |ILC (ILC3.LTI.CD4+)                    |CL:0001071 |
|GSM1585326_EA07068_305550_MOGENE-1_0-ST-V1_ILC3.LTI.CD4-.SI_4.CEL                |ILC               |ILC (ILC3.LTI.CD4-)                    |CL:0001071 |
|GSM1585329_EA07068_267952_MOGENE-1_0-ST-V1_ILC3.LTI.4+.SI_1.CEL                  |ILC               |ILC (ILC3.LTI.4+)                      |CL:0001071 |
|GSM1585330_EA07068_339254_MOGENE-1_0-ST-V1_NK.CD127-.SP_1.CEL                    |NK cells          |NK cells (NK.CD127-)                   |CL:0001065 |
|GSM1585333_EA07068_339239_MOGENE-1_0-ST-V1_LIV.NK.DX5+_1.CEL                     |ILC               |ILC (LIV.NK.DX5+)                      |CL:0001065 |
|GSM1585336_EA07068_339242_MOGENE-1_0-ST-V1_LPL.NCR+CNK_1.CEL                     |ILC               |ILC (LPL.NCR+CNK)                      |CL:0001065 |
|GSM2112407_EA07068_388554_MOGENE-1_0-ST-V1_BA.BL_1.CEL                           |Basophils         |Basophils (BA)                         |CL:0000767 |
|GSM2112413_EA07068_397997_MOGENE-1_0-ST-V1_Ep.5wk.MEC.Sca1+.Th_1.CEL             |Epithelial cells  |Epithelial cells (Ep.5wk.MEC.Sca1+)    |CL:0000066 |
|GSM2112415_EA07068_397999_MOGENE-1_0-ST-V1_Ep.5wk.MEChi.Th_2.CEL                 |Epithelial cells  |Epithelial cells (Ep.5wk.MEChi)        |CL:0000066 |
|GSM2112416_EA07068_397996_MOGENE-1_0-ST-V1_Ep.5wk.MEClo.Th_1.CEL                 |Epithelial cells  |Epithelial cells (Ep.5wk.MEClo)        |CL:0000066 |
|GSM2112418_EA07068_398003_MOGENE-1_0-ST-V1_Ep.8wk.CEC.Sca1+.Th_1.CEL             |Epithelial cells  |Epithelial cells (Ep.8wk.CEC.Sca1+)    |CL:0000066 |
|GSM2112420_EA07068_398002_MOGENE-1_0-ST-V1_Ep.8wk.CEChi.Th_1.CEL                 |Epithelial cells  |Epithelial cells (Ep.8wk.CEChi)        |CL:0000066 |
|GSM2112422_EA07068_398004_MOGENE-1_0-ST-V1_Ep.8wk.MEChi.Th_1.CEL                 |Epithelial cells  |Epithelial cells (Ep.8wk.MEChi)        |CL:0000066 |
|GSM2112424_EA07068_398005_MOGENE-1_0-ST-V1_Ep.8wk.MEClo.Th_1.CEL                 |Epithelial cells  |Epithelial cells (Ep.8wk.MEClo)        |CL:0000066 |
|GSM2112426_EA07068_388553_MOGENE-1_0-ST-V1_MC.ES_1.CEL                           |Mast cells        |Mast cells (MC.ES)                     |CL:0000097 |
|GSM2112428_EA07068_339312_MOGENE-1_0-ST-V1_MAST.PC_2.CEL                         |Mast cells        |Mast cells (MC)                        |CL:0000097 |
|GSM2112437_EA07068_354402_MOGENE-1_0-ST-V1_MC.TO_1.CEL                           |Mast cells        |Mast cells (MC.TO)                     |CL:0000097 |
|GSM2112440_EA07068_388549_MOGENE-1_0-ST-V1_MC.TR_1.CEL                           |Mast cells        |Mast cells (MC.TR)                     |CL:0000097 |
|GSM2112443_EA07068_449869_MOGENE-1_0-ST-V1_MC.DIGEST.PC_1.CEL                    |Mast cells        |Mast cells (MC.DIGEST)                 |CL:0000097 |
|GSM2112446_EA07068_201145_MOGENE-1_0-ST-V1_MECHI.GFP+.ADULT_6.CEL                |Epithelial cells  |Epithelial cells (MECHI.GFP+.ADULT)    |CL:0000066 |
|GSM2112449_EA07068_201151_MOGENE-1_0-ST-V1_MECHI.GFP+.ADULT.KO_1.CEL             |Epithelial cells  |Epithelial cells (MECHI.GFP+.ADULT.KO) |CL:0000066 |
|GSM2112452_EA07068_201148_MOGENE-1_0-ST-V1_MECHI.GFP-.ADULT_6.CEL                |Epithelial cells  |Epithelial cells (MECHI.GFP-.ADULT)    |CL:0000066 |
|GSM2112455_EA07068_307792_MOGENE-1_0-ST-V1_MF.480HI.LV.NAIVE_1.CEL               |Macrophages       |Macrophages (MF.480HI.NAIVE)           |CL:0000235 |
|GSM2112458_EA07068_307793_MOGENE-1_0-ST-V1_MF.480INT.LV.NAIVE_1.CEL              |Macrophages       |Macrophages (MF.480INT.NAIVE)          |CL:0000235 |
|GSM2112461_EA07068_235599_MOGENE-1_0-ST-V1_T.4EFF49D+11A+.SP.D8.LCMV.CEL         |T cells           |T cells (T.4EFF49D+11A+.D8.LCMV)       |CL:0001044 |
|GSM2112463_EA07068_235601_MOGENE-1_0-ST-V1_T.4MEM49D+11A+.SP.D30.LCMV.CEL        |T cells           |T cells (T.4MEM49D+11A+.D30.LCMV)      |CL:0000897 |
|GSM2112465_EA07068_235603_MOGENE-1_0-ST-V1_T.4NVE44-49D-11A-.SP.CEL              |T cells           |T cells (T.4NVE44-49D-11A-)            |CL:0000895 |
|GSM2112467_EA07068_349158_MOGENE-1_0-ST-V1_T.8EFF.TBET+.SP.OT1.D6LISOVA_1.CEL    |T cells           |T cells (T.8EFF.TBET+.OT1LISOVA)       |CL:0001050 |
|GSM2112470_EA07068_349161_MOGENE-1_0-ST-V1_T.8EFF.TBET-.SP.OT1.D6LISOVA_1.CEL    |T cells           |T cells (T.8EFF.TBET-.OT1LISOVA)       |CL:0001050 |
|GSM2112473_EA07068_311873_MOGENE-1_0-ST-V1_T.8EFFKLRG1+CD127-.SP.D8.LISOVA_2.CEL |T cells           |T cells (T.8EFFKLRG1+CD127-.D8.LISOVA) |CL:0001050 |
|GSM2112475_EA07068_311875_MOGENE-1_0-ST-V1_T.8MEMKLRG1-CD127+.SP.D8.LISOVA_1.CEL |T cells           |T cells (T.8MEMKLRG1-CD127+.D8.LISOVA) |CL:0000909 |
|GSM399362_EA07068_56648_MoGene_T.4+8int.Th_#1.cel                                |T cells           |T cells (T.4+8int)                     |CL:0002431 |
|GSM399365_EA07068_55678_MoGene_T.4FP3+25+.Sp_#2.cel                              |T cells           |T cells (T.4FP3+25+)                   |CL:0000792 |
|GSM399367_EA07068_56651_MoGene_T.4int8+.Th_#1.cel                                |T cells           |T cells (T.4int8+)                     |CL:0002430 |
|GSM399370_EA07068_52774_MoGene_T.4SP24-.Th_#1.cel                                |T cells           |T cells (T.4SP24-)                     |CL:0000624 |
|GSM399373_EA07068_52777_MoGene_T.4SP24int.Th_#1.cel                              |T cells           |T cells (T.4SP24int)                   |CL:0000624 |
|GSM399376_EA07068_52768_MoGene_T.4SP69+.Th_#1.cel                                |T cells           |T cells (T.4SP69+)                     |CL:0000896 |
|GSM399379_EA07068_52780_MoGene_T.8SP24-.Th_#1.cel                                |T cells           |T cells (T.8SP24-)                     |CL:0000625 |
|GSM399382_EA07068_52783_MoGene_T.8SP24int.Th_#1.cel                              |T cells           |T cells (T.8SP24int)                   |CL:0000625 |
|GSM399385_EA07068_52771_MoGene_T.8SP69+.Th_#1.cel                                |T cells           |T cells (T.8SP69+)                     |CL:0000906 |
|GSM399397_EA07068_56645_MoGene_T.DPbl.Th_#1.cel                                  |T cells           |T cells (T.DPbl)                       |CL:0002428 |
|GSM399400_EA07068_56642_MoGene_T.DPsm.Th_#1.cel                                  |T cells           |T cells (T.DPsm)                       |CL:0000809 |
|GSM399403_EA07068_52786_MoGene_T.ISP.Th_#1.cel                                   |T cells           |T cells (T.ISP)                        |CL:0000084 |
|GSM399438_EA07068_54191_MoGene_B.FrE.BM_#2.cel                                   |B cells           |B cells (B.FrE)                        |CL:0002054 |
|GSM399440_EA07068_54192_MoGene_B.FrF.BM_#2.cel                                   |B cells           |B cells (B.FrF)                        |CL:0002056 |
|GSM399448_EA07068_52806_MoGene_preB.FrD.BM_#1.cel                                |B cells           |B cells (preB.FrD)                     |CL:0002052 |
|GSM399450_EA07068_52803_MoGene_proB.FrBC.BM_#1.cel                               |B cells           |B cells (proB.FrBC)                    |CL:0002400 |
|GSM399452_EA07068_54189_MoGene_preB.FrC.BM_#2.cel                                |B cells           |B cells (preB.FrC)                     |CL:0002049 |
|GSM399454_EA07068_80000_MoGene_CD150-CD48-.BM#1.CEL                              |Stem cells        |Stem cells (SC.STSL)                   |CL:0000034 |
|GSM403986_EA07068_81316_MoGene_CD4+TESTNA.CEL                                    |T cells           |T cells (T.CD4+TESTNA)                 |CL:0000624 |
|GSM403987_EA07068_81315_MoGene_CD4+TESTDB.CEL                                    |T cells           |T cells (T.CD4+TESTDB)                 |CL:0000624 |
|GSM403988_EA07068_54833_MoGene_CD19CONTROL_#2.cel                                |B cells           |B cells (B.CD19CONTROL)                |CL:0000236 |
|GSM403994_EA07068_54832_MoGene_CD4CONTROL_#2.cel                                 |T cells           |T cells (T.CD4CONTROL)                 |CL:0000624 |
|GSM404000_EA07068_82676_MoGene_CD4TESTJS#1.CEL                                   |T cells           |T cells (T.CD4TESTJS)                  |CL:0000624 |
|GSM404003_EA07068_82674_MoGene_CD4TESTCJ#2.CEL                                   |T cells           |T cells (T.CD4TESTCJ)                  |CL:0000624 |
|GSM476654_EA07068_80001_MoGene_CD150-CD48-.BM#2.CEL                              |Stem cells        |Stem cells (SC.CD150-CD48-)            |CL:0000034 |
|GSM476655_EA07068_54199_MoGene_immTgd.vg2+.Th_#1.cel                             |Tgd               |Tgd (Tgd.imm.vg2+)                     |CL:0000799 |
|GSM476660_EA07068_56601_MoGene_immTgd.vg2.e17.Th_#2.cel                          |Tgd               |Tgd (Tgd.imm.vg2)                      |CL:0000799 |
|GSM476664_EA07068_56603_MoGene_matTgd.vg3.e17.Th_#1.cel                          |Tgd               |Tgd (Tgd.mat.vg3)                      |CL:0000800 |
|GSM476665_EA07068_56604_MoGene_matTgd.vg3.e17.Th.#2.cel                          |Tgd               |Tgd (Tgd.mat.vg3.)                     |CL:0000800 |
|GSM476672_EA07068_87590_MoGene_TGD.SP#1.CEL                                      |Tgd               |Tgd (Tgd)                              |CL:0000798 |
|GSM476678_EA07068_54193_MoGene_Tgd.vg2+.act.Sp_#1.cel                            |Tgd               |Tgd (Tgd.vg2+.act)                     |CL:0000798 |
|GSM476681_EA07068_54196_MoGene_Tgd.vg2-.act.Sp_#1.cel                            |Tgd               |Tgd (Tgd.vg2-.act)                     |CL:0000798 |
|GSM476684_EA07068_54550_MoGene_Tgd.vg2-.Sp_#1.cel                                |Tgd               |Tgd (Tgd.vg2-)                         |CL:0000798 |
|GSM538198_EA07068_56621_MoGene_B.Fo.PC_#1.CEL                                    |B cells           |B cells (B.Fo)                         |CL:0000843 |
|GSM538204_EA07068_80055_MoGene_B.FRE.FL#1.CEL                                    |B cells           |B cells (B.FRE)                        |CL:0000236 |
|GSM538207_EA07068_80057_MoGene_B.GC.SP#1.CEL                                     |B cells           |B cells (B.GC)                         |CL:0000844 |
|GSM538210_EA07068_56627_MoGene_B.MZ.Sp_#1.CEL                                    |B cells           |B cells (B.MZ)                         |CL:0000845 |
|GSM538213_EA07068_56630_MoGene_B.T1.Sp_#1.CEL                                    |B cells           |B cells (B.T1)                         |CL:0000958 |
|GSM538216_EA07068_56633_MoGene_B.T2.Sp_#1.CEL                                    |B cells           |B cells (B.T2)                         |CL:0000959 |
|GSM538219_EA07068_56636_MoGene_B.T3.Sp_#1.CEL                                    |B cells           |B cells (B.T3)                         |CL:0000960 |
|GSM538222_EA07068_56615_MoGene_B1a.PC_#1.CEL                                     |B cells           |B cells (B1a)                          |CL:0000820 |
|GSM538228_EA07068_56618_MoGene_B1b.PC_#1.CEL                                     |B cells           |B cells (B1b)                          |CL:0000821 |
|GSM538231_EA07068_87581_MoGene_DC2.LU#1.CEL                                      |DC                |DC (DC)                                |CL:0000451 |
|GSM538234_EA07068_96463_MoGene_DC.103+11B-.LV#1.CEL                              |DC                |DC (DC.103+11B-)                       |CL:0002506 |
|GSM538263_EA07068_96434_MoGene_DC.8-4-11B+.MLN#4.CEL                             |DC                |DC (DC.8-4-11B+)                       |CL:0002454 |
|GSM538280_EA07068_111375_MoGene_DC.LC.SK#4.CEL                                   |DC                |DC (DC.LC)                             |CL:0000451 |
|GSM538285_EA07068_96472_MoGene_NK.49CI+.SP#1@N2.CEL                              |NK cells          |NK cells (NK.49CI+)                    |CL:0000623 |
|GSM538288_EA07068_96475_MoGene_NK.49CI-.SP#1@N2.CEL                              |NK cells          |NK cells (NK.49CI-)                    |CL:0000623 |
|GSM538291_EA07068_96478_MoGene_NK.B2M-.SP#1.CEL                                  |NK cells          |NK cells (NK.B2M-)                     |CL:0000623 |
|GSM538294_EA07068_93784_MoGene_NK.DAP10-.SP#1.CEL                                |NK cells          |NK cells (NK.DAP10-)                   |CL:0000623 |
|GSM538297_EA07068_99792_MoGene_NK.DAP12-.SP#1.CEL                                |NK cells          |NK cells (NK.DAP12-)                   |CL:0000623 |
|GSM538300_EA07068_99749_MoGene_NK.H+.MCMV1.SP#1.CEL                              |NK cells          |NK cells (NK.H+.MCMV1)                 |CL:0000623 |
|GSM538303_EA07068_99755_MoGene_NK.H+.MCMV7.SP#1.CEL                              |NK cells          |NK cells (NK.H+.MCMV7)                 |CL:0000623 |
|GSM538309_EA07068_87578_MoGene_NK.H+MCMV1#1.CEL                                  |NK cells          |NK cells (NK.H+MCMV1)                  |CL:0000623 |
|GSM538312_EA07068_90292_MoGene_NK.MCMV7#1.CEL                                    |NK cells          |NK cells (NK.MCMV7)                    |CL:0000623 |
|GSM538315_EA07068_86161_MoGene_NK.SP#7.CEL                                       |NK cells          |NK cells (NK)                          |CL:0000623 |
|GSM538318_EA07068_91097_MoGene_NKT.4+.LV#1.CEL                                   |NKT               |NKT (NKT.4+)                           |CL:0000923 |
|GSM538325_EA07068_91101_MoGene_NKT.4-.LV#1.CEL                                   |NKT               |NKT (NKT.4-)                           |CL:0000924 |
|GSM538332_EA07068_91103_MoGene_NKT.44+NK1.1+.TH#1.CEL                            |NKT               |NKT (NKT.44+NK1.1+)                    |CL:0002438 |
|GSM538335_EA07068_91105_MoGene_NKT.44+NK1.1-.TH#1.CEL                            |NKT               |NKT (NKT.44+NK1.1-)                    |CL:0002041 |
|GSM538338_EA07068_96453_MoGene_NKT.44-NK1.1-.TH#1.CEL                            |NKT               |NKT (NKT.44-NK1.1-)                    |CL:0002040 |
|GSM538340_EA07068_80056_MoGene_PREB.FRD.FL#1.CEL                                 |B cells           |B cells (preB.FRD)                     |CL:0000817 |
|GSM538343_EA07068_52801_MoGene_proB.CLP.BM_#1.CEL                                |B cells           |B cells (proB.CLP)                     |CL:0000051 |
|GSM538346_EA07068_88784_MoGene_CLP#5.CEL                                         |Stem cells        |Stem cells (proB.CLP)                  |CL:0000051 |
|GSM538351_EA07068_52802_MoGene_proB.FrA.BM_#1.CEL                                |B cells           |B cells (proB.FrA)                     |CL:0002045 |
|GSM538352_EA07068_81297_MoGene_PROB.FRA.BM#4.CEL                                 |B cells           |B cells (proB.FRA)                     |CL:0002045 |
|GSM538353_EA07068_88783_MoGene_FRA#5.CEL                                         |B cells, pro      |B cells, pro (proB.FrA)                |CL:0002045 |
|GSM538362_EA07068_85523_MoGene_T.4MEM.LN#1.CEL                                   |T cells           |T cells (T.4MEM)                       |CL:0000897 |
|GSM538365_EA07068_58854_MoGene_T.4Mem.Sp_#1.CEL                                  |T cells           |T cells (T.4Mem)                       |CL:0000897 |
|GSM538368_EA07068_96415_MoGene_T.4MEM44H62L.LN#1.CEL                             |T cells           |T cells (T.4MEM44H62L)                 |CL:0000897 |
|GSM538374_EA07068_52756_MoGene_T.4Nve.LN_#1.CEL                                  |T cells           |T cells (T.4Nve)                       |CL:0000895 |
|GSM538380_EA07068_83933_MoGene_T.4NVE.PP#1.CEL                                   |T cells           |T cells (T.4NVE)                       |CL:0000895 |
|GSM538385_EA07068_80031_MoGene_AG#8.CEL                                          |T cells           |T cells (T.8EFF.OT1.D15.VSVOVA)        |CL:0001050 |
|GSM538387_EA07068_85512_MoGene_T.8EFF.SP.OT1.D5.VSVOVA#1.CEL                     |T cells           |T cells (T.8EFF.OT1.D5.VSVOVA)         |CL:0001050 |
|GSM538389_EA07068_80026_MoGene_AG#1.CEL                                          |T cells           |T cells (T.8EFF.OT1.VSVOVA)            |CL:0001050 |
|GSM538392_EA07068_80029_MoGene_AG#5.CEL                                          |T cells           |T cells (T.8EFF.OT1.D8.VSVOVA)         |CL:0001050 |
|GSM538395_EA07068_85520_MoGene_T.8MEM.LN#1.CEL                                   |T cells           |T cells (T.8MEM)                       |CL:0000909 |
|GSM538398_EA07068_58857_MoGene_T.8Mem.Sp_#1.CEL                                  |T cells           |T cells (T.8Mem)                       |CL:0000909 |
|GSM538401_EA07068_85518_MoGene_T.8MEM.SP.OT1.D106.VSVOVA#2.CEL                   |T cells           |T cells (T.8MEM.OT1.D106.VSVOVA)       |CL:0000909 |
|GSM538403_EA07068_80032_MoGene_AG#9.CEL                                          |T cells           |T cells (T.8EFF.OT1.D45VSV)            |CL:0001050 |
|GSM538406_EA07068_52759_MoGene_T.8Nve.LN_#1.CEL                                  |T cells           |T cells (T.8Nve)                       |CL:0000900 |
|GSM538412_EA07068_83936_MoGene_T.8NVE.PP#1.CEL                                   |T cells           |T cells (T.8NVE)                       |CL:0000900 |
|GSM538418_EA07068_81298_MoGene_PROB.FRBC.BM#4.CEL                                |B cells           |B cells (proB.FRBC)                    |CL:0000826 |
|GSM605753_EA07068_58851_MoGene_T.4.LN.BDC_#2.CEL                                 |T cells           |T cells (T.4)                          |CL:0000624 |
|GSM605756_EA07068_58845_MoGene_T.4.Pa.BDC_#2.CEL                                 |T cells           |T cells (T.4.Pa)                       |CL:0000624 |
|GSM605758_EA07068_58848_MoGene_T.4.PLN.BDC_#1.CEL                                |T cells           |T cells (T.4.PLN)                      |CL:0000624 |
|GSM605766_EA07068_55683_MoGene_T.4FP3-.Sp_#1.CEL                                 |T cells           |T cells (T.4FP3-)                      |CL:0000624 |
|GSM605787_EA07068_96412_MoGene_TGD.VG2+.SP#4.CEL                                 |Tgd               |Tgd (Tgd.VG2+)                         |CL:0000798 |
|GSM605790_EA07068_54559_MoGene_Tgd.vg2+.Sp.TCRbko_#1.CEL                         |Tgd               |Tgd (Tgd.vg2+.TCRbko)                  |CL:0000798 |
|GSM605796_EA07068_54202_MoGene_Tgd.vg2-.Sp.TCRbko_#1.CEL                         |Tgd               |Tgd (Tgd.vg2-.TCRbko)                  |CL:0000798 |
|GSM605802_EA07068_56609_MoGene_Tgd.vg5+.act.IEL_#1.CEL                           |Tgd               |Tgd (Tgd.vg5+.act)                     |CL:0000798 |
|GSM605804_EA07068_81294_MoGene_TGD.VG5+.ACT.IEL.#4.CEL                           |Tgd               |Tgd (Tgd.VG5+.ACT)                     |CL:0000798 |
|GSM605805_EA07068_81291_MoGene_TGD.VG5+.IEL.#4.CEL                               |Tgd               |Tgd (Tgd.VG5+)                         |CL:0000798 |
|GSM605808_EA07068_56606_MoGene_Tgd.vg5-.act.IEL_#1.CEL                           |Tgd               |Tgd (Tgd.vg5-.act)                     |CL:0000798 |
|GSM605811_EA07068_81288_MoGene_TGD.VG5-.IEL.#4.CEL                               |Tgd               |Tgd (Tgd.VG5-)                         |CL:0000798 |
|GSM605814_EA07068_108027_MoGene_NK.49H+.SP#1.CEL                                 |NK cells          |NK cells (NK.49H+)                     |CL:0000623 |
|GSM605817_EA07068_108030_MoGene_NK.49H-.SP#1.CEL                                 |NK cells          |NK cells (NK.49H-)                     |CL:0000623 |
|GSM605828_EA07068_108118_MoGene_DC.8+.TH#1.CEL                                   |DC                |DC (DC.8+)                             |CL:0001000 |
|GSM605831_EA07068_108115_MoGene_DC.8-.TH#1.CEL                                   |DC                |DC (DC.8-)                             |CL:0002460 |
|GSM605836_EA07068_96401_MoGene_DC.8-4-11B-.MLN#6@N2.CEL                          |DC                |DC (DC.8-4-11B-)                       |CL:0000998 |
|GSM605840_EA07068_105309_MoGene_DC.PDC.8+.SP#1.CEL                               |DC                |DC (DC.PDC.8+)                         |CL:0002456 |
|GSM605843_EA07068_105312_MoGene_DC.PDC.8-.SP#1.CEL                               |DC                |DC (DC.PDC.8-)                         |CL:0002455 |
|GSM605850_EA07068_105224_MoGene_MF.II-480HI.PC#1.CEL                             |Macrophages       |Macrophages (MF.II-480HI)              |CL:0000235 |
|GSM605853_EA07068_105221_MoGene_MF.RP.SP#1.CEL                                   |Macrophages       |Macrophages (MF.RP)                    |CL:0000235 |
|GSM605856_EA07068_105233_MoGene_MF.THIO5.II+480INT.PC#1.CEL                      |Macrophages       |Macrophages (MFIO5.II+480INT)          |CL:0000235 |
|GSM605859_EA07068_105242_MoGene_MF.THIO5.II+480LO.PC#1.CEL                       |Macrophages       |Macrophages (MFIO5.II+480LO)           |CL:0000235 |
|GSM605862_EA07068_105239_MoGene_MF.THIO5.II-480HI.PC#1.CEL                       |Macrophages       |Macrophages (MFIO5.II-480HI)           |CL:0000235 |
|GSM605865_EA07068_105236_MoGene_MF.THIO5.II-480INT.PC#1.CEL                      |Macrophages       |Macrophages (MFIO5.II-480INT)          |CL:0000235 |
|GSM605868_EA07068_96442_MoGene_MO.6C+II+.BL#1.CEL                                |Monocytes         |Monocytes (MO.6C+II+)                  |CL:0002470 |
|GSM605872_EA07068_96439_MoGene_MO.6C+II-.BL#1.CEL                                |Monocytes         |Monocytes (MO.6C+II-)                  |CL:0002469 |
|GSM605878_EA07068_96448_MoGene_MO.6C-II+.BL#1.CEL                                |Monocytes         |Monocytes (MO.6C-II+)                  |CL:0002473 |
|GSM605884_EA07068_96447_MoGene_MO.6C-II-.BL#3.CEL                                |Monocytes         |Monocytes (MO.6C-II-)                  |CL:0002471 |
|GSM605886_EA07068_96450_MoGene_MO.6C-IIINT.BL#1.CEL                              |Monocytes         |Monocytes (MO.6C-IIINT)                |CL:0002472 |
|GSM605891_EA07068_82682_MoGene_T.8EFF.SP.OT1.D10LIS.CEL                          |T cells           |T cells (T.8EFF.OT1.D10LIS)            |CL:0001050 |
|GSM605892_EA07068_85511_MoGene_T.8EFF.SP.OT1.D10.LISOVA#2.CEL                    |T cells           |T cells (T.8EFF.OT1.D10.LISOVA)        |CL:0001050 |
|GSM605894_EA07068_82683_MoGene_T.8EFF.SP.OT1.D15LIS.CEL                          |T cells           |T cells (T.8EFF.OT1.D15LIS)            |CL:0001050 |
|GSM605895_EA07068_85510_MoGene_T.8EFF.SP.OT1.D15.LISOVA#2.CEL                    |T cells           |T cells (T.8EFF.OT1.D15.LISOVA)        |CL:0001050 |
|GSM605898_EA07068_82680_MoGene_T.8EFF.SP.OT1.D6LISO.CEL                          |T cells           |T cells (T.8EFF.OT1LISO)               |CL:0001050 |
|GSM605899_EA07068_85549_MoGene_T.8EFF.SP.OT1.D6.LISOVA#2.CEL                     |T cells           |T cells (T.8EFF.OT1.LISOVA)            |CL:0001050 |
|GSM605901_EA07068_82681_MoGene_T.8EFF.SP.OT1.D8LISO.CEL                          |T cells           |T cells (T.8EFF.OT1.D8LISO)            |CL:0001050 |
|GSM605902_EA07068_85509_MoGene_T.8EFF.SP.OT1.D8.LISOVA#2.CEL                     |T cells           |T cells (T.8EFF.OT1.D8.LISOVA)         |CL:0001050 |
|GSM605904_EA07068_85517_MoGene_T.8MEM.SP.OT1.D100.LISOVA#1.CEL                   |T cells           |T cells (T.8MEM.OT1.D100.LISOVA)       |CL:0000909 |
|GSM605907_EA07068_85516_MoGene_T.8MEM.SP.OT1.D45.LISOVA#1.CEL                    |T cells           |T cells (T.8MEM.OT1.D45.LISOVA)        |CL:0000909 |
|GSM605909_EA07068_105264_MoGene_T.8NVE.SP.OT1#3.CEL                              |T cells           |T cells (T.8NVE.OT1)                   |CL:0000900 |
|GSM777019_EA07068_124592_MOGENE-1_0-ST-V1_B.FO.LN_1.CEL                          |B cells           |B cells (B.FO)                         |CL:0000843 |
|GSM777032_EA07068_108045_MoGene_BEC.MLN_3.CEL                                    |Endothelial cells |Endothelial cells (BEC)                |CL:0000115 |
|GSM777041_EA07068_81324_MoGene_EP.MECHI.TH_2.CEL                                 |Epithelial cells  |Epithelial cells (EP.MECHI)            |CL:0000066 |
|GSM777043_EA07068_81329_MoGene_FI.MTS15+.TH_1.CEL                                |Fibroblasts       |Fibroblasts (FI.MTS15+)                |CL:0000057 |
|GSM777046_EA07068_110672_MoGene_FI.SK_1.CEL                                      |Fibroblasts       |Fibroblasts (FI)                       |CL:0000057 |
|GSM777067_EA07068_121816_MOGENE-1_0-ST-V1_ST.31-38-44-.SLN_1.CEL                 |Stromal cells     |Stromal cells (ST.31-38-44-)           |CL:0000499 |
|GSM791102_EA07068_142883_MOGENE-1_0-ST-V1_SC.LT34F.BM_1.CEL                      |Stem cells        |Stem cells (SC.LT34F)                  |CL:0000034 |
|GSM791105_EA07068_140220_MOGENE-1_0-ST-V1_SC.MDP.BM_1.CEL                        |Stem cells        |Stem cells (SC.MDP)                    |CL:0002009 |
|GSM791108_EA07068_130473_MOGENE-1_0-ST-V1_SC.MEP.BM_1.CEL                        |Stem cells        |Stem cells (SC.MEP)                    |CL:0000050 |
|GSM791110_EA07068_130475_MOGENE-1_0-ST-V1_SC.MPP34F.BM_1.CEL                     |Stem cells        |Stem cells (SC.MPP34F)                 |CL:0000837 |
|GSM791112_EA07068_130477_MOGENE-1_0-ST-V1_SC.ST34F.BM_1.CEL                      |Stem cells        |Stem cells (SC.ST34F)                  |CL:0000034 |
|GSM791114_EA07068_140217_MOGENE-1_0-ST-V1_SC.CDP.BM_1.CEL                        |Stem cells        |Stem cells (SC.CDP)                    |CL:0000034 |
|GSM791117_EA07068_130471_MOGENE-1_0-ST-V1_SC.CMP.BM.DR_1.CEL                     |Stem cells        |Stem cells (SC.CMP.DR)                 |CL:0000049 |
|GSM791119_EA07068_111380_MoGene_GMP.BM_1.CEL                                     |Stem cells        |Stem cells (GMP)                       |CL:0000557 |
|GSM791124_EA07068_54184_MoGene_MLP.BM__1.cel                                     |Stem cells        |Stem cells (MLP)                       |CL:0000037 |
|GSM791126_EA07068_80048_MoGene_LTHSC.FL_1.CEL                                    |Stem cells        |Stem cells (LTHSC)                     |CL:0000034 |
|GSM791134_EA07068_110598_MoGene_T.DN2-3.TH_2.CEL                                 |T cells           |T cells (T.DN2-3)                      |CL:0002489 |
|GSM791136_EA07068_110595_MoGene_T.DN2.TH_4.CEL                                   |T cells           |T cells (T.DN2)                        |CL:0000806 |
|GSM791139_EA07068_117726_MOGENE-1_0-ST-V1_T.DN2A.TH_1.CEL                        |T cells           |T cells (T.DN2A)                       |CL:0000806 |
|GSM791141_EA07068_117728_MOGENE-1_0-ST-V1_T.DN2B.TH_1.CEL                        |T cells           |T cells (T.DN2B)                       |CL:0000806 |
|GSM791143_EA07068_110601_MoGene_T.DN3-4.TH_1.CEL                                 |T cells           |T cells (T.DN3-4)                      |CL:0002489 |
|GSM791146_EA07068_110599_MoGene_T.DN3A.TH_1.CEL                                  |T cells           |T cells (T.DN3A)                       |CL:0000807 |
|GSM791149_EA07068_110600_MoGene_T.DN3B.TH_1.CEL                                  |T cells           |T cells (T.DN3B)                       |CL:0000807 |
|GSM791152_EA07068_110653_MoGene_T.DN1-2.TH_3.CEL                                 |T cells           |T cells (T.DN1-2)                      |CL:0002489 |
|GSM791154_EA07068_110602_MoGene_T.DN4.TH_4.CEL                                   |T cells           |T cells (T.DN4)                        |CL:0000808 |
|GSM854258_EA07068_116124_MOGENE-1_0-ST-V1_DC.103-11B+.SALM3.SI_1.CEL             |Macrophages       |Macrophages (MF.103-11B+.SALM3)        |CL:0000235 |
|GSM854262_EA07068_105273_MoGene_DC.103-11B+.SI_1.CEL                             |Macrophages       |Macrophages (MF.103-11B+)              |CL:0000235 |
|GSM854269_EA07068_142689_MOGENE-1_0-ST-V1_DC.103-11B+24+.LU_1_N2.CEL             |DC                |DC (DC.103-11B+24+)                    |CL:0002505 |
|GSM854271_EA07068_142687_MOGENE-1_0-ST-V1_DC.103-11B+24-.LU_1_N2.CEL             |Macrophages       |Macrophages (MF.103-11B+24-)           |CL:0000235 |
|GSM854273_EA07068_140199_MOGENE-1_0-ST-V1_DC.103-11B+F4-80LO.KD_1.CEL            |DC                |DC (DC.103-11B+F4-80LO.KD)             |CL:0002505 |
|GSM854276_EA07068_116128_MOGENE-1_0-ST-V1_DC.11CLOSER.SALM3.SI_1.CEL             |Macrophages       |Macrophages (MF.11CLOSER.SALM3)        |CL:0000235 |
|GSM854280_EA07068_105277_MoGene_DC.11CLOSER.SI_1.CEL                             |Macrophages       |Macrophages (MF.11CLOSER)              |CL:0000235 |
|GSM854283_EA07068_108771_MoGene_DC.103CLOSER.SI_4.CEL                            |Macrophages       |Macrophages (MF.103CLOSER)             |CL:0000235 |
|GSM854294_EA07068_105226_MoGene_DC.II+480LO.PC_1.CEL                             |Macrophages       |Macrophages (MF.II+480LO)              |CL:0000235 |
|GSM854303_EA07068_121819_MOGENE-1_0-ST-V1_GN.ARTH.BM_1.CEL                       |Neutrophils       |Neutrophils (GN.ARTH)                  |CL:0000775 |
|GSM854309_EA07068_124598_MOGENE-1_0-ST-V1_GN.THIO.PC_1.CEL                       |Neutrophils       |Neutrophils (GN.Thio)                  |CL:0000775 |
|GSM854312_EA07068_121825_MOGENE-1_0-ST-V1_GN.URAC.PC_1.CEL                       |Neutrophils       |Neutrophils (GN.URAC)                  |CL:0000775 |
|GSM854315_EA07068_140214_MOGENE-1_0-ST-V1_MF.169+11CHI.SLN_1.CEL                 |Macrophages       |Macrophages (MF.169+11CHI)             |CL:0000235 |
|GSM854322_EA07068_140211_MOGENE-1_0-ST-V1_MF.MEDL.SLN_1.CEL                      |Macrophages       |Macrophages (MF.MEDL)                  |CL:0000235 |
|GSM854324_EA07068_140209_MOGENE-1_0-ST-V1_MF.SBCAPS.SLN_2.CEL                    |Macrophages       |Macrophages (MF.SBCAPS)                |CL:0000235 |
|GSM854326_EA07068_111383_MoGene_MICROGLIA.CNS_1.CEL                              |Microglia         |Microglia (Microglia)                  |CL:0000129 |
|GSM854335_EA07068_110652_MoGene_T.ETP.TH_6.CEL                                   |T cells           |T cells (T.ETP)                        |CL:0002425 |
|GSM920616_EA07068_108089_MoGene_IMMTGD.VG1+.TH.B6_1.CEL                          |Tgd               |Tgd (Tgd.imm.VG1+)                     |CL:0002414 |
|GSM920619_EA07068_108092_MoGene_IMMTGD.VG1+VD6+.TH.B6_1.CEL                      |Tgd               |Tgd (Tgd.imm.VG1+VD6+)                 |CL:0002415 |
|GSM920622_EA07068_108084_MoGene_MATTGD.VG1+.TH.B6_1.CEL                          |Tgd               |Tgd (Tgd.mat.VG1+)                     |CL:0002411 |
|GSM920624_EA07068_108086_MoGene_MATTGD.VG1+VD6+.TH.B6_1.CEL                      |Tgd               |Tgd (Tgd.mat.VG1+VD6+)                 |CL:0002416 |
|GSM920627_EA07068_114326_MOGENE-1_0-ST-V1_MATTGD.VG2+.TH_1.CEL                   |Tgd               |Tgd (Tgd.mat.VG2+)                     |CL:0002407 |
|GSM920629_EA07068_140223_MOGENE-1_0-ST-V1_TGD.VG3+24AHI.E17.TH_1.CEL             |Tgd               |Tgd (Tgd.VG3+24AHI)                    |CL:0000798 |
|GSM920632_EA07068_142881_MOGENE-1_0-ST-V1_TGD.VG5+24AHI.TH_1.CEL                 |Tgd               |Tgd (Tgd.VG5+24AHI)                    |CL:0000798 |
|GSM920634_EA07068_130429_MOGENE-1_0-ST-V1_T.8EFF.SP.OT1.12HR.LISOVA_1.CEL        |T cells           |T cells (T.8EFF.OT1.12HR.LISOVA)       |CL:0001050 |
|GSM920637_EA07068_130430_MOGENE-1_0-ST-V1_T.8EFF.SP.OT1.24HR.LISOVA_1.CEL        |T cells           |T cells (T.8EFF.OT1.24HR.LISOVA)       |CL:0001050 |
|GSM920640_EA07068_130432_MOGENE-1_0-ST-V1_T.8EFF.SP.OT1.48HR.LISOVA_1.CEL        |T cells           |T cells (T.8EFF.OT1.48HR.LISOVA)       |CL:0001050 |
|GSM920642_EA07068_105198_MoGene_B614WABDTREG_1.CEL                               |T cells           |T cells (T.Tregs)                      |CL:0000815 |
|GSM920648_EA07068_201208_MOGENE-1_0-ST-V1_TGD.VG2+24AHI.E17.TH_1.CEL             |Tgd               |Tgd (Tgd.VG2+24AHI)                    |CL:0000798 |
|GSM920651_EA07068_201205_MOGENE-1_0-ST-V1_TGD.VG4+24AHI.E17.TH_1.CEL             |Tgd               |Tgd (Tgd.VG4+24AHI)                    |CL:0000798 |
|GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL             |Tgd               |Tgd (Tgd.VG4+24ALO)                    |CL:0000798 |
</details>

<details>
  <summary>`MouseRNAseqData` Labels</summary>


|                                |label.main        |label.fine            |label.ont  |
|:-------------------------------|:-----------------|:---------------------|:----------|
|ERR525589Aligned                |Adipocytes        |Adipocytes            |CL_0000136 |
|PGE_young_EAligned              |Neurons           |aNSCs                 |CL:0000047 |
|SRR1033783Aligned               |Astrocytes        |Astrocytes            |CL:0000127 |
|SRR2938973Aligned               |Astrocytes        |Astrocytes activated  |CL:0000127 |
|SRR1033795Aligned               |Endothelial cells |Endothelial cells     |CL:0000115 |
|SRR1536428Aligned               |Erythrocytes      |Erythrocytes          |CL:0000232 |
|SRR1390714Aligned               |Fibroblasts       |Fibroblasts           |CL:0000057 |
|SRR1015752Aligned               |Fibroblasts       |Fibroblasts activated |CL:0000057 |
|SRR832851Aligned                |Fibroblasts       |Fibroblasts senescent |CL:0000057 |
|SRR1536401Aligned               |Granulocytes      |Granulocytes          |CL:0000094 |
|SRR1536397Aligned               |Macrophages       |Macrophages           |CL:0000235 |
|SRR1033793Aligned               |Microglia         |Microglia             |CL:0000129 |
|SRR2082382Aligned               |Microglia         |Microglia activated   |CL:0000129 |
|SRR1536407Aligned               |Monocytes         |Monocytes             |CL:0000576 |
|SRR1033785Aligned               |Neurons           |Neurons               |CL:0000540 |
|SRR2938959Aligned               |Neurons           |Neurons activated     |CL:0000540 |
|SRR1536422Aligned               |NK cells          |NK cells              |CL:0000623 |
|E_young_CAligned                |Neurons           |NPCs                  |CL:0002319 |
|SRR1033791Aligned               |Oligodendrocytes  |Oligodendrocytes      |CL:0000128 |
|PG_young_DAligned               |Neurons           |qNSCs                 |CL:0000047 |
|SRR1536413Aligned               |T cells           |T cells               |CL:0000084 |
|SRR2040609Aligned               |Dendritic cells   |Dendritic cells       |CL:0000451 |
|Cardiomyocyte_pseudo_Bulk       |Cardiomyocytes    |Cardiomyocytes        |CL:0000746 |
|Hepatocyte_pooled_Bulk2         |Hepatocytes       |Hepatocytes           |CL:0000182 |
|SRR1536411Aligned               |B cells           |B cells               |CL:0000236 |
|Ependymal_Striatum_pseudoBulk_1 |Epithelial cells  |Ependymal             |CL:0000065 |
|OPCs_pseudoBulk_1               |Oligodendrocytes  |OPCs                  |CL:0002453 |
|SRR1044039Aligned               |Macrophages       |Macrophages activated |CL:0000890 |
</details>

