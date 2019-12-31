## Human Moving Picture from MicrobiomeAnalyst server

library(MicrobiomeAnalystR)
library(phyloseq)

download.file("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/data/treebiom.zip",
  "data-raw/caporaso.zip"
)
unzip("data-raw/caporaso.zip", exdir = "data-raw/")
file.rename("data-raw/treebiom/", "data-raw/caporaso/")

ps <- import_biom(
  "data-raw/caporaso/otu_table_mc2_w_tax_no_pynast_failures.biom",
  treefilename = "data-raw/caporaso/rep_set.tre",
)

colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species")

sampledata <- read.delim("data-raw/caporaso/map.txt", row.names = 1) %>%
  sample_data()
caporaso_phyloseq <- merge_phyloseq(ps, sampledata)

usethis::use_data(caporaso_phyloseq, overwrite = TRUE)
