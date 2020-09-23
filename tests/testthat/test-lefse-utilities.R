context("test for lefse utilities")

test_that("check whether taxa is summarized or not", {
  expect_true(check_tax_summarize(oxygen))
  expect_true(check_tax_summarize(otu_table(oxygen)))
  expect_false(check_tax_summarize(pediatric_ibd))
  expect_false(check_tax_summarize(otu_table(pediatric_ibd)))
})
