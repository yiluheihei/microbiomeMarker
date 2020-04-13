
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R build
status](https://github.com/yiluheihei/microbiomeMarker/workflows/R-CMD-check/badge.svg)](https://github.com/yiluheihei/microbiomeMarker/actions)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/yiluheihei/microbiomeMarker/blob/master/LICENSE.md)
[![Codecov test
coverage](https://codecov.io/gh/yiluheihei/microbiomeMarker/branch/master/graph/badge.svg)](https://codecov.io/gh/yiluheihei/microbiomeMarker?branch=master)
[![DOI](https://zenodo.org/badge/215731961.svg)](https://zenodo.org/badge/latestdoi/215731961)
<!-- badges: end -->

## Motivation

**The aim of this package is to build a unified toolbox in R for
mcirobiome biomarker discovery by integrating various existing
methods.**

Many statistical methods have been proposed to discovery the microbiome
biomaker by compare the taxon abundance between different classes. Some
methods developed specifically for microbial community, such as linear
discriminant analysis (LDA) effect size (LEfSe) (Segata et al. 2011),
metagenomeSeq (Paulson et al. 2013); and some methods developed
specifically for RNA-Seq data, such as DESeq2 (Love, Huber, and Anders
2014) and edgeR (Robinson, McCarthy, and Smyth 2009), have been proposed
for microbiome biomarker discovery. We usually use several methods for
microbiome biomarker discovery and compare the results, which requires
multiple tools developed in different programming, even in different OS.

**microbiomeMarker** take the `phyloseq-class` object in package
[phyloseq](https://github.com/joey711/phyloseq) as input, since
**phyloseq** is the most popular R package in microbiome analysis and
with phyloseq you can easily import taxon abundance and phylogenetic
tree of taxon output from common microbiome bioinformatics platforms,
such as [DADA2](https://benjjneb.github.io/dada2/) and
[qiime2](https://qiime2.org/).

## Installation

You can install the package directly from github

``` r
if (!require(remotes)) install.packages("remotes")
remotes::install_github("yiluheihei/microbiomeMarker")
```

## LEfSe

Curently, LEfSe is the most used tool for microbiome biomarker
discovery, and the first method to integrate to **microbiomeMarker** is
LEfSe.

### lefse analysis

``` r
library(microbiomeMarker)
#> Registered S3 method overwritten by 'treeio':
#>   method     from
#>   root.phylo ape
library(ggplot2)

# sample data from lefse python script. The dataset contains 30 abundance 
# profiles (obtained processing the 16S reads with RDP) belonging to 10 rag2 
# (control) and 20 truc (case) mice
data("spontaneous_colitis")
lefse_out <- lefse(
  spontaneous_colitis, 
  normalization = 1e6, 
  class = "class", 
  multicls_strat = TRUE
)
# lefse return a microbioMarker class inherits from phyloseq
lefse_out
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 29 microbiome markers with 5 variables ]
#> otu_table()   OTU Table:         [ 132 taxa and  30 samples ]
#> tax_table()   Taxonomy Table:    [ 132 taxa by 1 taxonomic ranks ]
```

The microbiome biomarker information was stored in a new data structure
`marker_table-class` inherited from `data.frame`, and you can access it
by using `marker_table()`.

``` r
marker_table(lefse_out)
#>                                                                                                  feature
#> 1                                                                                 Bacteria|Bacteroidetes
#> 2                                                                     Bacteria|Bacteroidetes|Bacteroidia
#> 3                                                       Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales
#> 4            Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium
#> 5                                    Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae
#> 6                            Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae
#> 7                                               Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales
#> 8                                                                                Bacteria|Actinobacteria
#> 9                                                                 Bacteria|Actinobacteria|Actinobacteria
#> 10                       Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella
#> 11                                                                Bacteria|Firmicutes|Bacilli|Bacillales
#> 12                                              Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae
#> 13                               Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus
#> 14                   Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides
#> 15                                Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia
#> 16             Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia
#> 17                            Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter
#> 18                                                          Bacteria|Firmicutes|Clostridia|Clostridiales
#> 19                                                                        Bacteria|Firmicutes|Clostridia
#> 20                                                                                   Bacteria|Firmicutes
#> 21                                          Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae
#> 22                                          Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae
#> 23                            Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter
#> 24                                                           Bacteria|Proteobacteria|Gammaproteobacteria
#> 25 Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella
#> 26                                          Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae
#> 27                            Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus
#> 28             Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter
#> 29           Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia
#>    enrich_group log_max_mean      lda      p_value
#> 1          rag2     5.451241 5.178600 1.553428e-02
#> 2          rag2     5.433686 5.178501 1.375221e-02
#> 3          rag2     5.433686 5.178501 1.375221e-02
#> 4          rag2     5.082944 5.044767 1.217981e-04
#> 5          rag2     4.987349 4.886991 1.320110e-03
#> 6          rag2     4.789752 4.750839 1.217981e-04
#> 7          rag2     4.789752 4.750839 1.217981e-04
#> 8          rag2     4.800377 4.743824 6.001591e-04
#> 9          rag2     4.800377 4.743824 6.001591e-04
#> 10         rag2     4.761159 4.645092 1.320110e-03
#> 11         rag2     4.045278 3.839820 4.863507e-03
#> 12         rag2     4.007237 3.788714 7.256517e-03
#> 13         rag2     3.990162 3.770290 8.270483e-03
#> 14         rag2     3.652599 3.454102 7.056806e-03
#> 15         rag2     3.259877 3.225737 6.151322e-06
#> 16         rag2     2.982350 2.884262 4.204414e-04
#> 17         rag2     2.616605 2.571572 7.997905e-05
#> 18         truc     5.826884 5.455265 5.098324e-04
#> 19         truc     5.827542 5.453962 5.098324e-04
#> 20         truc     5.893569 5.391350 3.659080e-04
#> 21         truc     5.326873 4.943669 2.400614e-03
#> 22         truc     4.945682 4.512788 3.688840e-03
#> 23         truc     4.033392 3.639159 4.773320e-02
#> 24         truc     3.273543 3.310445 4.368471e-02
#> 25         truc     3.216059 3.259396 1.240061e-02
#> 26         truc     3.354679 3.156800 6.378989e-03
#> 27         truc     3.350311 3.149519 6.378989e-03
#> 28         truc     3.103242 2.914455 3.430807e-02
#> 29         truc     2.002354 2.021579 2.378651e-02
```

### Visualization of the result of lefse analysis

Bar plot for output of lefse:

``` r
lefse_barplot(lefse_out, label_level = 1) +
  scale_fill_manual(values = c("rag2" = "blue", "truc" = "red"))
```

![](man/figures/README-lefse-barplot-1.png)<!-- -->

Cladogram plot for output of lefse

``` r
lefse_cladogram(lefse_out, color = c("blue", "red"))
#> Warning: `data_frame()` is deprecated as of tibble 1.1.0.
#> Please use `tibble()` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
```

![](man/figures/README-lefse-cladogram-1.png)<!-- -->

## Welcome

**microbiomeMarker is still a newborn, and only contains lefse methods
right now. Your suggestion and contribution will be highly
appreciated.**

## Citation

Kindly cite as follows: Yang Cao (2020). microbiomeMarker: microbiome
biomarker analysis. R package version 0.0.1.9000.
<https://github.com/yiluheihei/microbiomeMarker>. DOI:
[10.5281/zenodo.3749415](https://doi.org/10.5281/zenodo.3749415).

## Acknowledgement

  - [lefse python
    script](https://bitbucket.org/biobakery/biobakery/wiki/lefse), The
    main lefse code are translated from **lefse python script**,
  - [microbiomeViz](https://github.com/lch14forever/microbiomeViz),
    cladogram visualization of lefse is modified from **microbiomeViz**.
  - [phyloseq](https://github.com/joey711/phyloseq), the main data
    structures used in **microbiomeMarker** are from or inherit from
    `phyloseq-class` in package **phyloseq**.

## Reference

<div id="refs" class="references">

<div id="ref-Love_2014">

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12). Springer Science; Business Media LLC.
<https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Paulson_2013">

Paulson, Joseph N, O Colin Stine, H’ector Corrada Bravo, and Mihai Pop.
2013. “Differential Abundance Analysis for Microbial Marker-Gene
Surveys.” *Nature Methods* 10 (12). Springer Science; Business Media
LLC: 1200–1202. <https://doi.org/10.1038/nmeth.2658>.

</div>

<div id="ref-Robinson_2009">

Robinson, M. D., D. J. McCarthy, and G. K. Smyth. 2009. “edgeR: A
Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1). Oxford University Press
(OUP): 139–40. <https://doi.org/10.1093/bioinformatics/btp616>.

</div>

<div id="ref-Segata_2011">

Segata, Nicola, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa
Miropolsky, Wendy S Garrett, and Curtis Huttenhower. 2011. “Metagenomic
Biomarker Discovery and Explanation.” *Genome Biology* 12 (6). Springer
Science; Business Media LLC: R60.
<https://doi.org/10.1186/gb-2011-12-6-r60>.

</div>

</div>
