library(ggplot2)
# phyloseq to the input of lefse galaxy -----------------------------------
lefse_in <- summarize_taxa(cid.phy)

meta <- sample_data(cid.phy)
all(names(lefse_in) == c("taxa", meta$Sample_ID))
names(lefse_in) <- c("taxa", meta$Consistency)

readr::write_tsv(lefse_in, "dev_test/lefse/lefse_in.tsv")
# format_input.py lefse_in.tsv lefse_in_formated.tsv -o 1000000
# run_lefse.py lefse_in_formated.tsv lefse_out.tsv

# filter enrich taxa
lefse_out <- readr::read_tsv("dev_test/lefse/lefse_out.tsv") %>%
  dplyr::filter(!is.na(X3)) %>%
  dplyr::mutate(X3 = sub(" ", "_", X3))
readr::write_tsv(
  lefse_out,
  "dev_test/lefse/lefse_out_filterd.tsv",
  col_names = FALSE
)
# plot_res.py --format pdf lefse_out_filterd.tsv lefse_res.pdf
names(lefse_out) <- c("taxa", "log_mean", "group", "lda", "p")
lefse_out <- dplyr::filter(
  lefse_out,
  lda > 4,
  group %in% c("formed_stool", "liquid")
) %>%
  dplyr::arrange(group, lda)
ggplot(lefse_out) +
  geom_col(aes(taxa, lda, fill = group)) +
  scale_x_discrete(limits = lefse_out$taxa) +
  coord_flip()

