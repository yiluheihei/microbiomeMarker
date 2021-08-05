# lefse - lda
mm_lefse <- run_lefse(
  kostic_crc,
  wilcoxon_cutoff = 0.01,
  group = "DIAGNOSIS",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4,
)

set.seed(2020)

# two groups --------------------------------------------------------------

# welch test - diff_mean
mm_welch <- run_test_two_groups(enterotypes_arumugam, "Gender")

# t test - diff_mean
mm_t <- run_test_two_groups(
  enterotypes_arumugam,
  group = "Gender",
  method = "t.test"
)

# white test - diff_mean
# mm_white <- test_two_groups(
#   enterotypes_arumugam,
#   group = "Gender",
#   method = "white.test",
#   nperm = 50
# )


# multiple groups ---------------------------------------------------------

enterotype <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)

# eta squared
mm_anova <- run_test_multiple_groups(
  enterotype,
  group = "Enterotype",
  method = "anova",
  effect_size_cutoff = 0.7
)

mm_kruskal <-  run_test_multiple_groups(
  enterotype,
  group = "Enterotype",
  method = "kruskal"
)

mm_mgs <- run_metagenomeseq(
  pediatric_ibd,
  "Class",
  contrast = c("CD","Control"),
  pvalue_cutoff = 0.1,
  p_adjust = "fdr"
)

mm_des <- run_deseq2(
  pediatric_ibd,
  "Class",
  contrast = c("Control", "CD"),
  pvalue_cutoff = 0.05,
  p_adjust = "fdr"
)

mm_edger <- run_edger(
  pediatric_ibd,
  "Class",
  c("CD", "Control"),
  pvalue_cutoff = 0.1,
  p_adjust = "fdr"
)

# 5 digits of decimal places to make sure the decimal results are equal
round_DF <- function(DF) {
  round2 <- function(x) {
    ifelse(
      x <= 1e-5,
      as.numeric(formatC(x, format = "g", digits = 5)),
      as.numeric(formatC(x, format = "f", digits = 5))
    )
  }

  df_rounded <- purrr::map_if(as.data.frame(DF), is.numeric, round2) %>%
    dplyr::bind_cols() %>%
    as.data.frame()
  row.names(df_rounded) <- row.names(DF)

  df_rounded
}


# a small phyloseq object for test
sample_ps <- phyloseq::subset_taxa(
  pediatric_ibd,
  Genus %in% c("g__Akkermansia", "g__Lactobacillus")
)
