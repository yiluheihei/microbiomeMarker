
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

**If you think this project is helpful to you, you can give this project
a star** :star:

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
2014) and edgeR \[@{Robinson\_2009\], have been proposed for microbiome
biomarker discovery. We usually use several methods for microbiome
biomarker discovery and compare the results, which requires multiple
tools developed in different programming, even in different OS.

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

## Data import

Since [phyloseq](https://github.com/joey711/phyloseq) objects are a
great data-standard for microbiome data in R, the core functions in
**microbiomeMarker** take `phylosq` object as input. Conveniently,
**microbiomeMarker** provides features to import external data files
form two common tools of microbiome analysis,
[qiime2](http://qiime.org/) and
[dada2](https://benjjneb.github.io/dada2).

### Import from dada2

The output of the [dada2](https://benjjneb.github.io/dada2) pipeline is
a feature table of amplicon sequence variants (an ASV table): A matrix
with rows corresponding to samples and columns to ASVs, in which the
value of each entry is the number of times that ASV was observed in that
sample. This table is analogous to the traditional OTU table.
Conveniently, taxa names are saved as

``` r
library(microbiomeMarker)
#> Registered S3 method overwritten by 'treeio':
#>   method     from
#>   root.phylo ape

seq_tab <- readRDS(system.file("extdata", "dada2_seqtab.rds",
  package= "microbiomeMarker"))
tax_tab <- readRDS(system.file("extdata", "dada2_taxtab.rds",
 package= "microbiomeMarker"))
sam_tab <- read.table(system.file("extdata", "dada2_samdata.txt",
 package= "microbiomeMarker"), sep = "\t", header = TRUE, row.names = 1)
ps <- import_dada2(seq_tab = seq_tab, tax_tab = tax_tab, sam_tab = sam_tab)
ps
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 232 taxa and 20 samples ]
#> sample_data() Sample Data:       [ 20 samples by 4 sample variables ]
#> tax_table()   Taxonomy Table:    [ 232 taxa by 6 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 232 reference sequences ]
```

### Import from qiime2

[qiime2](http://qiime.org/) is the most widely used software for
metagenomic analysis. User can import the feature table, taxonomic
table, phylogenetic tree, representative sequence and sample metadata
from qiime2 using `import_qiime2()`.

``` r
otuqza_file <- system.file("extdata", "table.qza",package = "microbiomeMarker")
taxaqza_file <- system.file("extdata", "taxonomy.qza",package = "microbiomeMarker")
sample_file <- system.file(
  "extdata", "sample-metadata.tsv",
  package = "microbiomeMarker"
)
treeqza_file <- system.file("extdata", "tree.qza",package = "microbiomeMarker")
ps <- import_qiime2(
  otu_qza = otuqza_file, taxa_qza = taxaqza_file,
  sam_tab = sample_file, tree_qza = treeqza_file
)
ps
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 770 taxa and 34 samples ]
#> sample_data() Sample Data:       [ 34 samples by 9 sample variables ]
#> tax_table()   Taxonomy Table:    [ 770 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 770 tips and 768 internal nodes ]
```

### Other import functions reexport from phyloseq

**microbiomeMarker** reexports three import functions from **phyloseq**,
including `import_biom()`, `import_qiime()` and `import_mothur()`, to
help users to import data from [biom file](http://biom-format.org/), and
output from [qiime](http://www.qiime.org/) and
[mothur](http://www.mothur.org/). More details on these three import
functions can be see from
[here](https://joey711.github.io/phyloseq/import-data.html#the_import_family_of_functions).

Users can also import the external files into `phyloseq` object
manually. For more details on how to create `phyloseq` object from
manually imported data, please see [this
tutorial](http://joey711.github.io/phyloseq/import-data.html#manual).

## LEfSe

Curently, LEfSe is the most used tool for microbiome biomarker
discovery, and the first method to integrate to **microbiomeMarker** is
LEfSe.

### lefse analysis

``` r
library(ggplot2)

# sample data from lefse python script. The dataset contains 30 abundance 
# profiles (obtained processing the 16S reads with RDP) belonging to 10 rag2 
# (control) and 20 truc (case) mice
data("spontaneous_colitis")
# add prefix of ranks
spontaneous_colitis <- microbiomeMarker:::add_prefix_summarized(
  spontaneous_colitis,
  c("k", "p", "c", "o", "f", "g")
)
mm <- lefse(
  spontaneous_colitis, 
  normalization = 1e6, 
  class = "class", 
  multicls_strat = TRUE
)
# lefse return a microbioMarker class inherits from phyloseq
mm
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 29 microbiome markers with 5 variables ]
#> otu_table()   OTU Table:         [ 132 taxa and  30 samples ]
#> tax_table()   Taxonomy Table:    [ 132 taxa by 1 taxonomic ranks ]
```

The microbiome biomarker information was stored in a new data structure
`marker_table-class` inherited from `data.frame`, and you can access it
by using `marker_table()`.

``` r
head(marker_table(mm))
#>                                                                                                         feature
#> 1                                                                                  k__Bacteria|p__Bacteroidetes
#> 2                                                                   k__Bacteria|p__Bacteroidetes|c__Bacteroidia
#> 3                                                  k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales
#> 4 k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium
#> 5                            k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae
#> 6                    k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae
#>   enrich_group log_max_mean      lda      p_value
#> 1         rag2     5.451241 5.178600 0.0155342816
#> 2         rag2     5.433686 5.178501 0.0137522075
#> 3         rag2     5.433686 5.178501 0.0137522075
#> 4         rag2     5.082944 5.044767 0.0001217981
#> 5         rag2     4.987349 4.886991 0.0013201097
#> 6         rag2     4.789752 4.750839 0.0001217981
```

### Visualization of the result of lefse analysis

Bar plot for output of lefse:

``` r
lefse_barplot(mm, label_level = 1) +
  scale_fill_manual(values = c("rag2" = "blue", "truc" = "red"))
```

![](man/figures/README-lefse-barplot-1.png)<!-- -->

## statistical analysis (stamp)

STAMP (Parks et al. 2014) is a widely-used graphical software package
that provides “best pratices” in choose appropriate statisticalmethods
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
  method = "welch.test"
)

# three significantly differential genera (marker)
two_group_welch
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 3 microbiome markers with 10 variables ]
#> otu_table()   OTU Table:         [ 244 taxa and  39 samples ]
#> tax_table()   Taxonomy Table:    [ 244 taxa by 1 taxonomic ranks ]
# details of result of the three markers
head(marker_table(two_group_welch))
#>                                                                 feature
#> p__Firmicutes|g__Heliobacterium         p__Firmicutes|g__Heliobacterium
#> p__Firmicutes|g__Parvimonas                 p__Firmicutes|g__Parvimonas
#> p__Firmicutes|g__Peptostreptococcus p__Firmicutes|g__Peptostreptococcus
#>                                     enrich_group     pvalue F_mean_rel_freq
#> p__Firmicutes|g__Heliobacterium                M 0.02940341    0.0001038235
#> p__Firmicutes|g__Parvimonas                    M 0.03281399    0.0001911176
#> p__Firmicutes|g__Peptostreptococcus            M 0.01714937    0.0019799118
#>                                     M_mean_rel_freq     diff_mean      ci_lower
#> p__Firmicutes|g__Heliobacterium        0.0005309321 -0.0004271086 -0.0008080028
#> p__Firmicutes|g__Parvimonas            0.0008610460 -0.0006699283 -0.0012804196
#> p__Firmicutes|g__Peptostreptococcus    0.0053274345 -0.0033475227 -0.0060532040
#>                                          ci_upper ratio_proportion
#> p__Firmicutes|g__Heliobacterium     -4.621437e-05        0.1955495
#> p__Firmicutes|g__Parvimonas         -5.943705e-05        0.2219599
#> p__Firmicutes|g__Peptostreptococcus -6.418413e-04        0.3716445
#>                                     pvalue_corrected
#> p__Firmicutes|g__Heliobacterium           0.02940341
#> p__Firmicutes|g__Parvimonas               0.03281399
#> p__Firmicutes|g__Peptostreptococcus       0.01714937
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
  method = "anova"
)

# 22 markers (significantly differential genera)
multiple_group_anova
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 24 microbiome markers with 8 variables ]
#> otu_table()   OTU Table:         [ 238 taxa and  32 samples ]
#> tax_table()   Taxonomy Table:    [ 238 taxa by 1 taxonomic ranks ]
head(marker_table(multiple_group_anova))
#>                                                                 feature
#> p__Bacteroidetes                                       p__Bacteroidetes
#> p__Unclassified                                         p__Unclassified
#> p__Actinobacteria|g__Scardovia           p__Actinobacteria|g__Scardovia
#> p__Bacteroidetes|g__Alistipes             p__Bacteroidetes|g__Alistipes
#> p__Bacteroidetes|g__Bacteroides         p__Bacteroidetes|g__Bacteroides
#> p__Bacteroidetes|g__Parabacteroides p__Bacteroidetes|g__Parabacteroides
#>                                     enrich_group       pvalue pvalue_corrected
#> p__Bacteroidetes                    Enterotype 1 3.196070e-06     3.196070e-06
#> p__Unclassified                     Enterotype 3 1.731342e-04     1.731342e-04
#> p__Actinobacteria|g__Scardovia      Enterotype 2 2.742042e-02     2.742042e-02
#> p__Bacteroidetes|g__Alistipes       Enterotype 3 3.922758e-02     3.922758e-02
#> p__Bacteroidetes|g__Bacteroides     Enterotype 1 8.396825e-10     8.396825e-10
#> p__Bacteroidetes|g__Parabacteroides Enterotype 1 1.314233e-02     1.314233e-02
#>                                     effect_size
#> p__Bacteroidetes                      0.5821619
#> p__Unclassified                       0.4497271
#> p__Actinobacteria|g__Scardovia        0.2196652
#> p__Bacteroidetes|g__Alistipes         0.2001541
#> p__Bacteroidetes|g__Bacteroides       0.7633661
#> p__Bacteroidetes|g__Parabacteroides   0.2582573
#>                                     Enterotype 1:mean_rel_freq_percent
#> p__Bacteroidetes                                            19.3073387
#> p__Unclassified                                             16.5364988
#> p__Actinobacteria|g__Scardovia                               0.0000000
#> p__Bacteroidetes|g__Alistipes                                0.6695668
#> p__Bacteroidetes|g__Bacteroides                             17.4793538
#> p__Bacteroidetes|g__Parabacteroides                          0.9745028
#>                                     Enterotype 2:mean_rel_freq_percent
#> p__Bacteroidetes                                          15.372374444
#> p__Unclassified                                           25.019756042
#> p__Actinobacteria|g__Scardovia                             0.001860083
#> p__Bacteroidetes|g__Alistipes                              0.528789506
#> p__Bacteroidetes|g__Bacteroides                            3.409612826
#> p__Bacteroidetes|g__Parabacteroides                        0.405579500
#>                                     Enterotype 3:mean_rel_freq_percent
#> p__Bacteroidetes                                          7.046051e+00
#> p__Unclassified                                           2.680750e+01
#> p__Actinobacteria|g__Scardovia                            8.436111e-05
#> p__Bacteroidetes|g__Alistipes                             1.568063e+00
#> p__Bacteroidetes|g__Bacteroides                           4.456618e+00
#> p__Bacteroidetes|g__Parabacteroides                       4.401643e-01
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
#>  [1] "p__Bacteroidetes"                     
#>  [2] "p__Unclassified"                      
#>  [3] "p__Actinobacteria|g__Scardovia"       
#>  [4] "p__Bacteroidetes|g__Alistipes"        
#>  [5] "p__Bacteroidetes|g__Bacteroides"      
#>  [6] "p__Bacteroidetes|g__Parabacteroides"  
#>  [7] "p__Bacteroidetes|g__Prevotella"       
#>  [8] "p__Firmicutes|g__Bulleidia"           
#>  [9] "p__Firmicutes|g__Catenibacterium"     
#> [10] "p__Firmicutes|g__Catonella"           
#> [11] "p__Firmicutes|g__Holdemania"          
#> [12] "p__Firmicutes|g__Lactobacillus"       
#> [13] "p__Firmicutes|g__Macrococcus"         
#> [14] "p__Firmicutes|g__Peptostreptococcus"  
#> [15] "p__Firmicutes|g__Ruminococcus"        
#> [16] "p__Firmicutes|g__Selenomonas"         
#> [17] "p__Firmicutes|g__Streptococcus"       
#> [18] "p__Firmicutes|g__Subdoligranulum"     
#> [19] "p__Proteobacteria|g__Bartonella"      
#> [20] "p__Proteobacteria|g__Brucella"        
#> [21] "p__Proteobacteria|g__Granulibacter"   
#> [22] "p__Proteobacteria|g__Rhodospirillum"  
#> [23] "p__Proteobacteria|g__Stenotrophomonas"
#> [24] "p__Unclassified|g__Unclassified"
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

## Visulatiton

### Cladogram plot

``` r
plot_cladogram(mm, color = c("blue", "red"))
```

![](man/figures/README-cladogram-1.png)<!-- -->

It’s recommended to use a named vector to set the colors of enriched
group:

``` r
plot_cladogram(mm, color = c(truc = "blue", rag2 = "red"))
```

![](man/figures/README-cladogram-color-1.png)<!-- -->

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

-   [lefse python
    script](https://bitbucket.org/biobakery/biobakery/wiki/lefse), The
    main lefse code are translated from **lefse python script**,
-   [microbiomeViz](https://github.com/lch14forever/microbiomeViz),
    cladogram visualization of lefse is modified from **microbiomeViz**.
-   [phyloseq](https://github.com/joey711/phyloseq), the main data
    structures used in **microbiomeMarker** are from or inherit from
    `phyloseq-class` in package **phyloseq**.
-   [MicrobiotaProcess](https://github.com/YuLab-SMU/MicrobiotaProcess),
    function `import_dada2()` and `import_qiime2()` are modified from
    the `MicrobiotaProcess::import_dada2()`.
-   [qiime2R](https://github.com/jbisanz/qiime2R), `import_qiime2()` are
    refer to the functions in qiime2R.

## Question

If you have any question, please file an issue on the issue tracker
following the instructions in the issue template:

Please briefly describe your problem, what output actually happened, and
what output you expect.

Please provide a minimal reproducible example. For more details on how
to make a great minimal reproducible example, see
<https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example>
and <https://www.tidyverse.org/help/#reprex>.

    Brief description of the problem

    # insert minimal reprducible example here

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Love_2014" class="csl-entry">

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12). <https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Parks_2014" class="csl-entry">

Parks, Donovan H., Gene W. Tyson, Philip Hugenholtz, and Robert G.
Beiko. 2014. “STAMP: Statistical Analysis of Taxonomic and Functional
Profiles.” *Bioinformatics* 30 (21): 3123–24.
<https://doi.org/10.1093/bioinformatics/btu494>.

</div>

<div id="ref-Paulson_2013" class="csl-entry">

Paulson, Joseph N, O Colin Stine, H’ector Corrada Bravo, and Mihai Pop.
2013. “Differential Abundance Analysis for Microbial Marker-Gene
Surveys.” *Nature Methods* 10 (12): 1200–1202.
<https://doi.org/10.1038/nmeth.2658>.

</div>

<div id="ref-Segata_2011" class="csl-entry">

Segata, Nicola, Jacques Izard, Levi Waldron, Dirk Gevers, Larisa
Miropolsky, Wendy S Garrett, and Curtis Huttenhower. 2011. “Metagenomic
Biomarker Discovery and Explanation.” *Genome Biology* 12 (6): R60.
<https://doi.org/10.1186/gb-2011-12-6-r60>.

</div>

</div>
