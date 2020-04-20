context("tax summarize")

test_that("check whether phyloseq tax summarized or not", {
  expect_true(check_tax_summarize(oxygen))
  expect_false(check_tax_summarize(pediatric_ibd))
})

test_that("check the summarize_taxa", {
  skip_on_cran()
  sep = "|"
  taxa <- summarize_taxa(pediatric_ibd, sep = sep)
  taxa <- round(taxa, 7)

  expect_true(any(grepl(sep, row.names(taxa), fixed = TRUE)))
  expect_known_output(
    head(taxa, 10),
    test_path("out/test-summarize-taxa-out.txt"),
    print = TRUE
  )
})
