context("tax summarize")

test_that("check whether phyloseq tax summarized or not", {
  expect_true(check_tax_summarize(oxygen))
  expect_false(check_tax_summarize(pediatric_ibd))
})

test_that("check the summarize_taxa", {
  skip_on_cran()
  sep = "|"
  taxa <- summarize_taxa(pediatric_ibd, sep = "|")

  expect_true(any(grepl(sep, row.names(taxa), fixed = TRUE)))
  expect_equal(nrow(taxa), 828L)

  expect_error(
    summarize_taxa(enterotypes_arumugam, level = "abc"),
    "`level` must in the ranks"
  )
})
