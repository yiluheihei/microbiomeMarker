context("assignment methods")

test_that("otu_table<- works correctly", {
  expect_silent(otu_table(mm_lefse) <- otu_table(mm_lefse))
})
