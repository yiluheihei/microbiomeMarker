# availabel taxonomic ranks, Summarize represents summarized tax
availabel_ranks <- c("Summarize", "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species")
availabel_ranks <- factor(
  availabel_ranks,
  levels = availabel_ranks
)

usethis::use_data(availabel_ranks, internal = TRUE)
