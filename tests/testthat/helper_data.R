lefse_out <- lefse(
  oxygen,
  norm = 1e6,
  bootstrap_n = 5,
  # summarize = FALSE,
  class = "oxygen_availability",
  subclass = "body_site"
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
