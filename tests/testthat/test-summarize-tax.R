context("tax summarize")

test_that("check whether phyloseq tax summarized or not", {
  expect_true(check_tax_summarize(oxygen))
})

test_that("check the row names of summarize_taxa is the taxa name", {
  sep = "|"
  taxa <- summarize_taxa(pediatric_ibd, sep = sep) %>%
    row.names()

  expect_true(any(grepl(sep, taxa, fixed = TRUE)))
})
