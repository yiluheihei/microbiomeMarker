
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R build
status](https://github.com/yiluheihei/microbiomeMarker/workflows/R-CMD-check/badge.svg)](https://github.com/yiluheihei/microbiomeMarker/actions)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/yiluheihei/microbiomeMarker/blob/master/LICENSE.md)
<!-- badges: end -->

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
**The aim of this package is to build a unified toolbox in R for
mcirobiome biomarker discovery by integrating various existing
methods.**

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

``` r
library(microbiomeMarker)

# sample data
data("oxygen")
lefse_out <- lefse(
  oxygen, 
  normalization = 1e6, 
  summarize = "lefse", 
  class = "oxygen_availability", 
  subclass = "body_site", 
  multicls_strat = TRUE
)
#> Warning: Setting row names on a tibble is deprecated.
lefse_out
#>                                                                                         feature
#> 1                                                                       Bacteria|Actinobacteria
#> 2                                                        Bacteria|Actinobacteria|Actinobacteria
#> 3                                        Bacteria|Actinobacteria|Actinobacteria|Actinomycetales
#> 4                   Bacteria|Actinobacteria|Actinobacteria|Actinomycetales|Propionibacteriaceae
#> 5 Bacteria|Actinobacteria|Actinobacteria|Actinomycetales|Propionibacteriaceae|Propionibacterium
#> 6                               Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae
#> 7                   Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides
#> 8                                                   Bacteria|Firmicutes|Bacilli|Lactobacillales
#>   enrich_group log_max_mean      lda      p_value
#> 1      High_O2     6.426614 6.256256 4.148229e-09
#> 2      High_O2     5.826868 5.658028 4.148229e-09
#> 3      High_O2     5.826094 5.656212 3.841252e-09
#> 4      High_O2     5.697155 5.546813 9.483488e-10
#> 5      High_O2     5.697009 5.546671 1.034228e-09
#> 6       Low_O2     5.613862 5.431316 3.034294e-07
#> 7       Low_O2     5.613862 5.431316 3.034294e-07
#> 8       Mid_O2     5.655506 5.502118 5.889748e-08
```

Bar plot for output of lefse:

``` r
lefse_barplot(lefse_out, label_level = 1)
```

<img src="man/figures/README-lefse-barplot-1.png" width="100%" />

**Cladogram plot for output of lefse is in developing**

## Welcome

**microbiomeMarker is still a newborn, and only contains lefse methods
right now. Your suggestion and contribution will be highly
appreciated.**

## Reference

<div id="refs" class="references">

<div id="ref-Love_2014">

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12). <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Paulson_2013">

Paulson, Joseph N, O Colin Stine, H’ector Corrada Bravo, and Mihai Pop.
2013. “Differential Abundance Analysis for Microbial Marker-Gene
Surveys.” *Nature Methods* 10 (12): 1200–1202.
<https://doi.org/10.1038/nmeth.2658>.

</div>

<div id="ref-Robinson_2009">

Robinson, M. D., D. J. McCarthy, and G. K. Smyth. 2009. “edgeR: A
Bioconductor Package for Differential Expression Analysis of Digital
Gene Expression Data.” *Bioinformatics* 26 (1): 139–40.
<https://doi.org/10.1093/bioinformatics/btp616>.

</div>

<div id="ref-Segata_2011">

Segata, Nicola, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa
Miropolsky, Wendy S Garrett, and Curtis Huttenhower. 2011. “Metagenomic
Biomarker Discovery and Explanation.” *Genome Biology* 12 (6): R60.
<https://doi.org/10.1186/gb-2011-12-6-r60>.

</div>

</div>
