# available taxonomic ranks, Summarize represents summarized tax
available_ranks <- c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
)
available_ranks <- factor(
    available_ranks,
    levels = available_ranks
)

usethis::use_data(available_ranks, internal = TRUE, overwrite = TRUE)
