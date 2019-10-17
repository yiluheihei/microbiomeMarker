## Human Moving Picture from MicrobiomeAnalyst server

library(MicrobiomeAnalystR)
library(magrittr)

download.file("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/treebiom.zip",
  "data-raw/caporaso.zip"
)
unzip("data-raw/caporaso.zip", exdir = "data-raw/")
file.rename("data-raw/treebiom/", "data-raw/caporaso/")

mbset <- Init.mbSetObj() %>%
  SetModuleType("mdp") %>%
  Read16SAbundData("data-raw/caporaso/otu_table_mc2_w_tax_no_pynast_failures.biom","biom","GreengenesID","F") %>%
  ReadSampleTable("data-raw/caporaso/map.txt") %>%
  ReadTreeFile("data-raw/caporaso/rep_set.tre") %>%
  SanityCheckData("biom") %>%
  PlotLibSizeView("norm_libsizes_0","png") %>%
  CreatePhyloseqObj("biom","GreengenesID","F") %>%
  ApplyAbundanceFilter("prevalence", 4, 0.2) %>%
  ApplyVarianceFilter("iqr", 0.1) %>%
  PerformNormalization("none", "colsum", "none")

caporaso_phyloseq <- mbset$dataSet$norm.phyobj

usethis::use_data(caporaso_phyloseq)
