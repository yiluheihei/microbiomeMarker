context("assignment methods")

test_that("otu_table<- works correctly", {
  expect_silent(otu_table(mm_lefse) <- otu_table(mm_lefse))
  expect_silent(otu_table(mm_lefse) <- mm_lefse)
  ps_assign <- phyloseq(otu_table(mm_lefse), tax_table(mm_lefse))
  expect_silent(otu_table(mm_lefse) <- ps_assign)
})
