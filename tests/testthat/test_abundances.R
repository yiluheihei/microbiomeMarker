context("`abundances` methods")

test_that("abundances",{
  abd <- abundances(enterotypes_arumugam, "log10")
  abd_t <- transform_abundances(enterotypes_arumugam, "log10")

  expect_warning(abundances(enterotypes_arumugam, "log10"), "contains zeroes")
  expect_identical(abd, as(otu_table(abd_t), "matrix"))
})
