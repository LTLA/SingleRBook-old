--- 
title: "Assigning cell types with SingleR"
date: "2020-06-10"
site: bookdown::bookdown_site
documentclass: book
bibliography: ref.bib
biblio-style: apalike
link-citations: yes
description: "The SingleR book. Because sometimes, a vignette just isn't enough."
---



# Preface {-}

<a href="https://bioconductor.org"><img src="https://github.com/Bioconductor/BiocStickers/raw/master/Bioconductor/Bioconductor-serial.gif" width="200" alt="Bioconductor Sticker" align="right" style="margin: 0 1em 0 1em" /></a> 

Imagine a world without a reference genome.
Whenever you receive new sequencing data - and we'll just talk about transcriptomic data for now -
you'll have to run it through an assembler to characterize the genes that are being expressed.
If you're lucky enough to get a decent assembly that isn't confused by contamination or homologous sequences,
it falls on you to pick through the sequence motifs to figure out what what the gene actually is. 
And if you're at a conference and you hear the name of a gene you're working on,
how can you be sure that everyone's talking about the same thing?

Now, this particular scenario only exists history books but it's easy to see the parallels with single-cell data analyses;
simply replace reads with cells, assemblies with clusters, and genes with cell types.
In the single-cell analysis field, a typical practitioner will hope that their clusters are reasonable proxies for the biological states of interest (a strong assumption indeed!) and that their manual annotation of the clusters is accurate (which is itself dependent on domain expertise and coffee consumption).
The current process of clustering, looking at a handful of markers and making a guesstimate of the cell type is best described as "artisanal" - which is not inherently bad, but one would hope that a mature technology would not require so much manual intervention for its routine use.

If hard work is a disease, then automation is the cure, and this case is no exception.
Automated cell type annotation methods match cells in a new dataset against curated reference profiles of known cell types, assigning each cell to the type that its expression profile is most similar to. 
This allows users to skip the mundane annotation of their data and jump directly to the interesting questions - does my cell type change in abundance or expression across treatments? Is there interesting substructure within an existing population?
In this respect, automated annotation methods are the single-cell field's equivalent to genome aligners,
and one hopes that the former will become standard procedure for single-cell data analysis.

This book covers the use of *[SingleR](https://bioconductor.org/packages/3.12/SingleR)*, one implementation of an automated annotation method.
If you want a survey of different annotation methods - this book is not for you.
If you want to create hand-crafted cluster definitions - this book is not for you.
(Read the [other one](https://osca.bioconductor.org) instead.)
If you want to use the pre-Bioconductor version of the package - this book is not for you.
But if you're tired of manually annotating your single-cell data and you want to do something better with your life, then read on.
