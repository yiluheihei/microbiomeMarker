context("assignment methods")

test_that("otu_table<- works correctly", {
  expect_silent(otu_table(lefse_out) <- otu_table(lefse_out))
})
