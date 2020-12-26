
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

### Import from tab-delimited input file of biobakery lefse

For [biobakey lefse](https://huttenhower.sph.harvard.edu/lefse/) (a
[Galaxy module](http://huttenhower.sph.harvard.edu/galaxy), a Conda
formula, a Docker image, and included in bioBakery (VM and cloud).), the
input file must be a tab-delimited text, consists of a list of numerical
features, the class vector and optionally the subclass and subject
vectors. The features can be read counts directly or abundance
floating-point values more generally, and the first field is the name of
the feature. Class, subclass and subject vectors have a name (the first
field) and a list of non-numerical strings. [biobakery
lefse](https://huttenhower.sph.harvard.edu/lefse/). User can import the
input file suitable for [biobakery
lefse](https://huttenhower.sph.harvard.edu/lefse/) to `phyloseq` object
using `import_biobakery_lefse_in()`

``` r
file <- system.file(
  "extdata",
  "hmp_small_aerobiosis.txt",
  package = "microbiomeMarker"
)
# six level of taxonomic ranks,
# meta data: row 1 represents class (oxygen_availability),
# row 2 represents subclass (body_site), row 3 represents subject (subject_id)
hmp_oxygen <- import_biobakery_lefse_in(
  file,
  ranks_prefix = c("k", "p", "c", "o", "f", "g"),
  meta_rows = 1:3,
)
hmp_oxygen
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 928 taxa and 55 samples ]
#> sample_data() Sample Data:       [ 55 samples by 3 sample variables ]
#> tax_table()   Taxonomy Table:    [ 928 taxa by 1 taxonomic ranks ]
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
#>                                                                                                               feature
#> marker1                                                                                  k__Bacteria|p__Bacteroidetes
#> marker2                                                                   k__Bacteria|p__Bacteroidetes|c__Bacteroidia
#> marker3                                                  k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales
#> marker4 k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium
#> marker5                            k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae
#> marker6                    k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae
#>         enrich_group log_max_mean      lda      p_value
#> marker1         rag2     5.451241 5.178600 0.0155342816
#> marker2         rag2     5.433686 5.178501 0.0137522075
#> marker3         rag2     5.433686 5.178501 0.0137522075
#> marker4         rag2     5.082944 5.044767 0.0001217981
#> marker5         rag2     4.987349 4.886991 0.0013201097
#> marker6         rag2     4.789752 4.750839 0.0001217981
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
#>                                     feature enrich_group     pvalue
#> marker1     p__Firmicutes|g__Heliobacterium            M 0.02940341
#> marker2         p__Firmicutes|g__Parvimonas            M 0.03281399
#> marker3 p__Firmicutes|g__Peptostreptococcus            M 0.01714937
#>               F_mean       M_mean     diff_mean      ci_lower      ci_upper
#> marker1 1.038235e-06 5.309321e-06 -4.271086e-06 -8.080028e-06 -4.621437e-07
#> marker2 1.911176e-06 8.610460e-06 -6.699283e-06 -1.280420e-05 -5.943705e-07
#> marker3 1.979912e-05 5.327434e-05 -3.347523e-05 -6.053204e-05 -6.418413e-06
#>             ratio pvalue_corrected
#> marker1 0.1955495       0.02940341
#> marker2 0.2219599       0.03281399
#> marker3 0.3716445       0.01714937
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

# 24 markers
multiple_group_anova
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 24 microbiome markers with 8 variables ]
#> otu_table()   OTU Table:         [ 238 taxa and  32 samples ]
#> tax_table()   Taxonomy Table:    [ 238 taxa by 1 taxonomic ranks ]
head(marker_table(multiple_group_anova))
#>                                     feature enrich_group       pvalue
#> marker1                    p__Bacteroidetes Enterotype 1 3.196070e-06
#> marker2                     p__Unclassified Enterotype 3 1.731342e-04
#> marker3      p__Actinobacteria|g__Scardovia Enterotype 2 2.742042e-02
#> marker4       p__Bacteroidetes|g__Alistipes Enterotype 3 3.922758e-02
#> marker5     p__Bacteroidetes|g__Bacteroides Enterotype 1 8.396825e-10
#> marker6 p__Bacteroidetes|g__Parabacteroides Enterotype 1 1.314233e-02
#>         pvalue_corrected effect_size Enterotype 1:mean_abundance
#> marker1     3.196070e-06   0.5821619                 0.193073387
#> marker2     1.731342e-04   0.4497271                 0.165364988
#> marker3     2.742042e-02   0.2196652                 0.000000000
#> marker4     3.922758e-02   0.2001541                 0.006695668
#> marker5     8.396825e-10   0.7633661                 0.174793538
#> marker6     1.314233e-02   0.2582573                 0.009745028
#>         Enterotype 2:mean_abundance Enterotype 3:mean_abundance
#> marker1                1.537237e-01                7.046051e-02
#> marker2                2.501976e-01                2.680750e-01
#> marker3                1.860083e-05                8.436111e-07
#> marker4                5.287895e-03                1.568063e-02
#> marker5                3.409613e-02                4.456618e-02
#> marker6                4.055795e-03                4.401643e-03
```

The result of multiple group statistic specified whether the means of
all groups is equal or not. To identify which pairs of groups may differ
from each other, post-hoc test must be performed.

``` r
pht <- posthoc_test(ps, group = "Enterotype")
pht
#> postHocTest-class object
#> Pairwise test result of 238  features,  DataFrameList object, each DataFrame has five variables:
#>         comparions    : pair groups to test which separated by '-'
#>         diff_mean: difference in mean proportions
#>         pvalue        : post hoc test p values
#>         ci_lower : lower confidence interval
#>         ci_upper : upper confidence interval
#> Posthoc multiple comparisons of means  using  tukey  method

# 24 significantly differential genera
markers <- marker_table(multiple_group_anova)$feature
markers
#>                        p__Bacteroidetes                         p__Unclassified 
#>                      "p__Bacteroidetes"                       "p__Unclassified" 
#>          p__Actinobacteria|g__Scardovia           p__Bacteroidetes|g__Alistipes 
#>        "p__Actinobacteria|g__Scardovia"         "p__Bacteroidetes|g__Alistipes" 
#>         p__Bacteroidetes|g__Bacteroides     p__Bacteroidetes|g__Parabacteroides 
#>       "p__Bacteroidetes|g__Bacteroides"   "p__Bacteroidetes|g__Parabacteroides" 
#>          p__Bacteroidetes|g__Prevotella              p__Firmicutes|g__Bulleidia 
#>        "p__Bacteroidetes|g__Prevotella"            "p__Firmicutes|g__Bulleidia" 
#>        p__Firmicutes|g__Catenibacterium              p__Firmicutes|g__Catonella 
#>      "p__Firmicutes|g__Catenibacterium"            "p__Firmicutes|g__Catonella" 
#>             p__Firmicutes|g__Holdemania          p__Firmicutes|g__Lactobacillus 
#>           "p__Firmicutes|g__Holdemania"        "p__Firmicutes|g__Lactobacillus" 
#>            p__Firmicutes|g__Macrococcus     p__Firmicutes|g__Peptostreptococcus 
#>          "p__Firmicutes|g__Macrococcus"   "p__Firmicutes|g__Peptostreptococcus" 
#>           p__Firmicutes|g__Ruminococcus            p__Firmicutes|g__Selenomonas 
#>         "p__Firmicutes|g__Ruminococcus"          "p__Firmicutes|g__Selenomonas" 
#>          p__Firmicutes|g__Streptococcus        p__Firmicutes|g__Subdoligranulum 
#>        "p__Firmicutes|g__Streptococcus"      "p__Firmicutes|g__Subdoligranulum" 
#>         p__Proteobacteria|g__Bartonella           p__Proteobacteria|g__Brucella 
#>       "p__Proteobacteria|g__Bartonella"         "p__Proteobacteria|g__Brucella" 
#>      p__Proteobacteria|g__Granulibacter     p__Proteobacteria|g__Rhodospirillum 
#>    "p__Proteobacteria|g__Granulibacter"   "p__Proteobacteria|g__Rhodospirillum" 
#>   p__Proteobacteria|g__Stenotrophomonas         p__Unclassified|g__Unclassified 
#> "p__Proteobacteria|g__Stenotrophomonas"       "p__Unclassified|g__Unclassified"
# take a marker "p__Bacteroidetes|g__Bacteroides"  
# for example, we will show "p__Bacteroidetes|g__Bacteroides"  differ from 
# between Enterotype 2-Enterotype 1 and Enterotype 3-Enterotype 2.
pht@result$"p__Bacteroidetes|g__Bacteroides"
#> DataFrame with 3 rows and 5 columns
#>               comparions  diff_mean      pvalue   ci_lower   ci_upper
#>              <character>  <numeric>   <numeric>  <numeric>  <numeric>
#> 1 Enterotype 2-Enterot.. -0.2813948 4.77015e-08 -0.3713469 -0.1914428
#> 2 Enterotype 3-Enterot.. -0.2604547 1.63635e-09 -0.3312286 -0.1896808
#> 3 Enterotype 3-Enterot..  0.0209401 7.88993e-01 -0.0575765  0.0994567
```

Visualization of post test result of a given feature.

``` r
# visualize the post hoc test result of Bacteroides
plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Bacteroides")
```

![](man/figures/README-plot-posthoctest-1.png)<!-- -->

## metagenomeSeq

``` r
mm_mgs <- run_metagenomeseq(
  pediatric_ibd, 
  norm = "CSS",
  "Class", 
  "Control", 
  "CD", 
  p_value_cutoff = 0.1, 
  p_adjust = "fdr"
)
#> Default value being used.
mm_mgs
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 11 microbiome markers with 15 variables ]
#> otu_table()   OTU Table:         [ 786 taxa and  43 samples ]
#> tax_table()   Taxonomy Table:    [ 786 taxa by 1 taxonomic ranks ]
```

## DESeq2

``` r
mm_des <- run_deseq2(
  pediatric_ibd, 
  "Class", 
  "Control", 
  "CD", 
  p_value_cutoff = 0.05, 
  p_adjust = "fdr"
)
#> converting counts to integer mode
mm_des
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 36 microbiome markers with 8 variables ]
#> otu_table()   OTU Table:         [ 786 taxa and  43 samples ]
#> tax_table()   Taxonomy Table:    [ 786 taxa by 1 taxonomic ranks ]
```

## edgeR

``` r
mm_edger <- run_edger(
  pediatric_ibd, 
  "Class", "Control", "CD", 
  p_value_cutoff = 0.1,
  p_adjust = "fdr"
)
mm_edger
#> microbiomeMarker-class inherited from phyloseq-class
#> marker_table  Marker Table:      [ 34 microbiome markers with 8 variables ]
#> otu_table()   OTU Table:         [ 786 taxa and  43 samples ]
#> tax_table()   Taxonomy Table:    [ 786 taxa by 1 taxonomic ranks ]
```

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
