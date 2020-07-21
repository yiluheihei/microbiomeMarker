
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
#>                                                                                                                                                                                                     feature
#> Bacteria|Bacteroidetes                                                                                                                                                               Bacteria|Bacteroidetes
#> Bacteria|Bacteroidetes|Bacteroidia                                                                                                                                       Bacteria|Bacteroidetes|Bacteroidia
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales                                                                                                           Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium                     Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae                                                                     Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae                                                     Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales                                                                                           Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales
#> Bacteria|Actinobacteria                                                                                                                                                             Bacteria|Actinobacteria
#> Bacteria|Actinobacteria|Actinobacteria                                                                                                                               Bacteria|Actinobacteria|Actinobacteria
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella                                             Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella
#> Bacteria|Firmicutes|Bacilli|Bacillales                                                                                                                               Bacteria|Firmicutes|Bacilli|Bacillales
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae                                                                                           Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus                                                             Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides                                     Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia                                                               Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia                         Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter                                                       Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter
#> Bacteria|Firmicutes|Clostridia|Clostridiales                                                                                                                   Bacteria|Firmicutes|Clostridia|Clostridiales
#> Bacteria|Firmicutes|Clostridia                                                                                                                                               Bacteria|Firmicutes|Clostridia
#> Bacteria|Firmicutes                                                                                                                                                                     Bacteria|Firmicutes
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae                                                                                   Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae                                                                                   Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter                                                       Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter
#> Bacteria|Proteobacteria|Gammaproteobacteria                                                                                                                     Bacteria|Proteobacteria|Gammaproteobacteria
#> Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae                                                                                   Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus                                                       Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus
#> Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter                         Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter
#> Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia                     Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia
#>                                                                                                       enrich_group
#> Bacteria|Bacteroidetes                                                                                        rag2
#> Bacteria|Bacteroidetes|Bacteroidia                                                                            rag2
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales                                                              rag2
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium                   rag2
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae                                           rag2
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae                                   rag2
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales                                                      rag2
#> Bacteria|Actinobacteria                                                                                       rag2
#> Bacteria|Actinobacteria|Actinobacteria                                                                        rag2
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella                               rag2
#> Bacteria|Firmicutes|Bacilli|Bacillales                                                                        rag2
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae                                                      rag2
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus                                       rag2
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides                           rag2
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia                                        rag2
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia                     rag2
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter                                    rag2
#> Bacteria|Firmicutes|Clostridia|Clostridiales                                                                  truc
#> Bacteria|Firmicutes|Clostridia                                                                                truc
#> Bacteria|Firmicutes                                                                                           truc
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae                                                  truc
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae                                                  truc
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter                                    truc
#> Bacteria|Proteobacteria|Gammaproteobacteria                                                                   truc
#> Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella         truc
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae                                                  truc
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus                                    truc
#> Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter                     truc
#> Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia                   truc
#>                                                                                                       log_max_mean
#> Bacteria|Bacteroidetes                                                                                    5.451241
#> Bacteria|Bacteroidetes|Bacteroidia                                                                        5.433686
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales                                                          5.433686
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium               5.082944
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae                                       4.987349
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae                               4.789752
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales                                                  4.789752
#> Bacteria|Actinobacteria                                                                                   4.800377
#> Bacteria|Actinobacteria|Actinobacteria                                                                    4.800377
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella                           4.761159
#> Bacteria|Firmicutes|Bacilli|Bacillales                                                                    4.045278
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae                                                  4.007237
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus                                   3.990162
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides                       3.652599
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia                                    3.259877
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia                 2.982350
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter                                2.616605
#> Bacteria|Firmicutes|Clostridia|Clostridiales                                                              5.826884
#> Bacteria|Firmicutes|Clostridia                                                                            5.827542
#> Bacteria|Firmicutes                                                                                       5.893569
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae                                              5.326873
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae                                              4.945682
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter                                4.033392
#> Bacteria|Proteobacteria|Gammaproteobacteria                                                               3.273543
#> Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella     3.216059
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae                                              3.354679
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus                                3.350311
#> Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter                 3.103242
#> Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia               2.002354
#>                                                                                                            lda
#> Bacteria|Bacteroidetes                                                                                5.178600
#> Bacteria|Bacteroidetes|Bacteroidia                                                                    5.178501
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales                                                      5.178501
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium           5.044767
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae                                   4.886991
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae                           4.750839
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales                                              4.750839
#> Bacteria|Actinobacteria                                                                               4.743824
#> Bacteria|Actinobacteria|Actinobacteria                                                                4.743824
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella                       4.645092
#> Bacteria|Firmicutes|Bacilli|Bacillales                                                                3.839820
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae                                              3.788714
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus                               3.770290
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides                   3.454102
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia                                3.225737
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia             2.884262
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter                            2.571572
#> Bacteria|Firmicutes|Clostridia|Clostridiales                                                          5.455265
#> Bacteria|Firmicutes|Clostridia                                                                        5.453962
#> Bacteria|Firmicutes                                                                                   5.391350
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae                                          4.943669
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae                                          4.512788
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter                            3.639159
#> Bacteria|Proteobacteria|Gammaproteobacteria                                                           3.310445
#> Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella 3.259396
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae                                          3.156800
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus                            3.149519
#> Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter             2.914455
#> Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia           2.021579
#>                                                                                                            p_value
#> Bacteria|Bacteroidetes                                                                                1.553428e-02
#> Bacteria|Bacteroidetes|Bacteroidia                                                                    1.375221e-02
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales                                                      1.375221e-02
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Bifidobacterium           1.217981e-04
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae                                   1.320110e-03
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae                           1.217981e-04
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales                                              1.217981e-04
#> Bacteria|Actinobacteria                                                                               6.001591e-04
#> Bacteria|Actinobacteria|Actinobacteria                                                                6.001591e-04
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Barnesiella                       1.320110e-03
#> Bacteria|Firmicutes|Bacilli|Bacillales                                                                4.863507e-03
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae                                              7.256517e-03
#> Bacteria|Firmicutes|Bacilli|Bacillales|Staphylococcaceae|Staphylococcus                               8.270483e-03
#> Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Porphyromonadaceae|Parabacteroides                   7.056806e-03
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae|Roseburia                                6.151322e-06
#> Bacteria|Actinobacteria|Actinobacteria|Bifidobacteriales|Bifidobacteriaceae|Metascardovia             4.204414e-04
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Papillibacter                            7.997905e-05
#> Bacteria|Firmicutes|Clostridia|Clostridiales                                                          5.098324e-04
#> Bacteria|Firmicutes|Clostridia                                                                        5.098324e-04
#> Bacteria|Firmicutes                                                                                   3.659080e-04
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Lachnospiraceae                                          2.400614e-03
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae                                          3.688840e-03
#> Bacteria|Firmicutes|Clostridia|Clostridiales|Ruminococcaceae|Oscillibacter                            4.773320e-02
#> Bacteria|Proteobacteria|Gammaproteobacteria                                                           4.368471e-02
#> Bacteria|Proteobacteria|Gammaproteobacteria|Enterobacteriales|Enterobacteriaceae|Escherichia/Shigella 1.240061e-02
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae                                          6.378989e-03
#> Bacteria|Firmicutes|Bacilli|Lactobacillales|Streptococcaceae|Streptococcus                            6.378989e-03
#> Bacteria|Actinobacteria|Actinobacteria|Coriobacteriales|Coriobacteriaceae|Asaccharobacter             3.430807e-02
#> Bacteria|Proteobacteria|Deltaproteobacteria|Desulfovibrionales|Desulfovibrionaceae|Lawsonia           2.378651e-02
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
#> Warning: `tbl_df()` is deprecated as of dplyr 1.0.0.
#> Please use `tibble::as_tibble()` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
```

![](man/figures/README-lefse-cladogram-1.png)<!-- -->

## statistical analysis (stamp)

STAMP (Parks et al. 2014) is a widely-used graphical software package
that provides “best pratices” in choose appropriate statistical methods
for microbial taxonomic and functional analysis. Users can tests for
both two groups or multiple groups, and effect sizes and confidence
intervals are supported that allows critical assessment of the
biological relevancy of test results. Here, **microbiomeMarker** also
integrates the statistical methods used in STAMP for microbial
comparison analysis between two-groups and multiple-groups.

### Statitical analysis between two groups

Function `test_two_groups()` is developed for statistical test between
two groups, and three test methods are provided: welch test, t test and
white test.

``` r
data("enterotypes_arumugam")
# take welch test for example
two_group_welch <- test_two_groups(
  enterotypes_arumugam, 
  group = "Gender", 
  rank_name = "Genus", 
  method = "welch.test"
)

# three significantly differential genera (marker)
two_group_welch
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 3 microbiome markers with 9 variables ]
#> otu_table()   OTU Table:         [ 248 taxa and  39 samples ]
#> tax_table()   Taxonomy Table:    [ 248 taxa by 2 taxonomic ranks ]
# details of result of the three markers
marker_table(two_group_welch)
#>              feature     pvalue F_mean_rel_freq M_mean_rel_freq     diff_mean
#> 1         Parvimonas 0.03281399    0.0003822353     0.001722092 -0.0013398567
#> 2 Peptostreptococcus 0.01714937    0.0039598236     0.010654869 -0.0066950454
#> 3     Heliobacterium 0.02940341    0.0002076471     0.001061864 -0.0008542172
#>       ci_lower      ci_upper ratio_proportion pvalue_corrected
#> 1 -0.002560839 -1.188741e-04        0.2219599       0.03281399
#> 2 -0.012106408 -1.283683e-03        0.3716445       0.01714937
#> 3 -0.001616006 -9.242875e-05        0.1955495       0.02940341
```

### Statistical analysis multiple groups

Function `test_multiple_groups()` is constructed for statistical test
for multiple groups, two test method are provided: anova and kruskal
test.

``` r
# three groups
ps <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)

multiple_group_anova <-  test_multiple_groups(
  ps, 
  group = "Enterotype", 
  rank_name = "Genus",
  method = "anova"
)

# 22 markers (significantly differential genera)
multiple_group_anova
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 22 microbiome markers with 7 variables ]
#> otu_table()   OTU Table:         [ 248 taxa and  32 samples ]
#> tax_table()   Taxonomy Table:    [ 248 taxa by 2 taxonomic ranks ]
marker_table(multiple_group_anova)
#>               feature       pvalue pvalue_corrected effect_size
#> 1         Bacteroides 8.396825e-10     8.396825e-10   0.7633661
#> 2          Bartonella 7.531192e-03     7.531192e-03   0.2861996
#> 3            Brucella 3.063042e-02     3.063042e-02   0.2136846
#> 4       Granulibacter 1.354378e-02     1.354378e-02   0.2567165
#> 5         Macrococcus 1.522491e-02     1.522491e-02   0.2506944
#> 6      Rhodospirillum 1.198159e-03     1.198159e-03   0.3711917
#> 7       Lactobacillus 7.720685e-03     7.720685e-03   0.2849752
#> 8       Streptococcus 4.575949e-02     4.575949e-02   0.1916129
#> 9    Stenotrophomonas 1.430707e-02     1.430707e-02   0.2539007
#> 10        Selenomonas 2.785320e-02     2.785320e-02   0.2188220
#> 11    Subdoligranulum 2.938189e-02     2.938189e-02   0.2159381
#> 12       Ruminococcus 3.283096e-02     3.283096e-02   0.2099133
#> 13 Peptostreptococcus 1.447573e-02     1.447573e-02   0.2532975
#> 14          Catonella 1.866784e-02     1.866784e-02   0.2400848
#> 15          Bulleidia 2.361928e-02     2.361928e-02   0.2276547
#> 16    Catenibacterium 3.962497e-02     3.962497e-02   0.1995979
#> 17         Holdemania 2.684419e-02     2.684419e-02   0.2208074
#> 18    Parabacteroides 1.314233e-02     1.314233e-02   0.2582573
#> 19         Prevotella 9.036411e-09     9.036411e-09   0.7212332
#> 20          Alistipes 3.922758e-02     3.922758e-02   0.2001541
#> 21          Scardovia 2.742042e-02     2.742042e-02   0.2196652
#> 22       Unclassified 1.731342e-04     1.731342e-04   0.4497271
#>    Enterotype 1:mean_rel_freq_percent Enterotype 2:mean_rel_freq_percent
#> 1                         34.95870756                       6.819226e+00
#> 2                          0.00000000                       1.202167e-03
#> 3                          0.00000000                       7.033333e-04
#> 4                          0.00000000                       8.008333e-04
#> 5                          0.00000000                       6.470000e-04
#> 6                          0.00000000                       1.073167e-03
#> 7                          0.01077725                       3.101887e-01
#> 8                          0.21598412                       9.916455e-01
#> 9                          0.00248725                       0.000000e+00
#> 10                         0.00000000                       1.567833e-03
#> 11                         0.28972975                       2.877475e-01
#> 12                         0.33682250                       5.873363e-01
#> 13                         0.01626638                       4.286667e-03
#> 14                         0.01895975                       5.800833e-03
#> 15                         0.00015000                       2.309000e-03
#> 16                         0.00780900                       2.068505e-01
#> 17                         0.12643087                       3.562817e-02
#> 18                         1.94900562                       8.111590e-01
#> 19                         0.32851937                       2.203184e+01
#> 20                         1.33913364                       1.057579e+00
#> 21                         0.00000000                       3.720167e-03
#> 22                        33.07299754                       5.003951e+01
#>    Enterotype 3:mean_rel_freq_percent
#> 1                        8.913236e+00
#> 2                        0.000000e+00
#> 3                        0.000000e+00
#> 4                        0.000000e+00
#> 5                        2.072222e-05
#> 6                        7.238889e-05
#> 7                        1.979034e-02
#> 8                        1.294123e-01
#> 9                        1.435556e-04
#> 10                       7.883333e-05
#> 11                       4.773734e-01
#> 12                       1.909550e+00
#> 13                       5.057118e-03
#> 14                       1.088767e-02
#> 15                       3.867778e-04
#> 16                       5.548968e-02
#> 17                       1.127455e-01
#> 18                       8.803286e-01
#> 19                       1.147909e+00
#> 20                       3.136127e+00
#> 21                       1.687222e-04
#> 22                       5.361500e+01
```

The result of multiple group statistic specified whether the means of
all groups is equal or not. To identify which pairs of groups may differ
from each other, post-hoc test must be performed.

``` r
pht <- posthoc_test(ps, group = "Enterotype" , rank_name = "Genus")
pht
#> postHocTest-class object
#> Pairwise test result of 248  features,  DataFrameList object, each DataFrame has five variables:
#>         comparions    : pair groups to test which separated by '-'
#>         diff_mean_prop: difference in mean proportions
#>         pvalue        : post hoc test p values
#>         ci_lower_prop : lower confidence interval
#>         ci_upper_prop : upper confidence interval
#> Posthoc multiple comparisons of means  using  tukey  method

# 22 significantly differential genera
markers <- marker_table(multiple_group_anova)$feature
markers
#>                  sp1                 sp12                 sp15 
#>        "Bacteroides"         "Bartonella"           "Brucella" 
#>                 sp22                 sp25                 sp27 
#>      "Granulibacter"        "Macrococcus"     "Rhodospirillum" 
#>                 sp93                sp113                sp124 
#>      "Lactobacillus"      "Streptococcus"   "Stenotrophomonas" 
#>                sp149                sp153                sp157 
#>        "Selenomonas"    "Subdoligranulum"       "Ruminococcus" 
#>                sp174                sp180                sp196 
#> "Peptostreptococcus"          "Catonella"          "Bulleidia" 
#>                sp197                sp200                sp203 
#>    "Catenibacterium"         "Holdemania"    "Parabacteroides" 
#>                sp205                sp206                sp241 
#>         "Prevotella"          "Alistipes"          "Scardovia" 
#>                sp248 
#>       "Unclassified"
# take a marker Bacteroides for example, we will show Bacteroides differ from 
# between Enterotype 2-Enterotype 1 and Enterotype 3-Enterotype 2.
pht@result$Bacteroides
#> DataFrame with 3 rows and 5 columns
#>                  comparions diff_mean_prop      pvalue ci_lower_prop
#>                 <character>      <numeric>   <numeric>     <numeric>
#> 1 Enterotype 2-Enterotype 1      -28.13948 4.77015e-08     -37.13469
#> 2 Enterotype 3-Enterotype 1      -26.04547 1.63635e-09     -33.12286
#> 3 Enterotype 3-Enterotype 2        2.09401 7.88993e-01      -5.75765
#>   ci_upper_prop
#>       <numeric>
#> 1     -19.14428
#> 2     -18.96808
#> 3       9.94567
```

Visualization of post test result of a given feature.

``` r
# visualize the post hoc test result of Bacteroides
plot_postHocTest(pht, feature = "Bacteroides")
```

![](man/figures/README-plot-posthoctest-1.png)<!-- -->

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
*Genome Biology* 15 (12). <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Parks_2014">

Parks, Donovan H., Gene W. Tyson, Philip Hugenholtz, and Robert G.
Beiko. 2014. “STAMP: Statistical Analysis of Taxonomic and Functional
Profiles.” *Bioinformatics* 30 (21): 3123–4.
<https://doi.org/10.1093/bioinformatics/btu494>.

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
