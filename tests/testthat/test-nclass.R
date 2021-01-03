context("nclass method")

test_that("nclass method", {
  expect_equal(nclass(mm_lefse), 2)
  expect_equal(nclass(marker_table(mm_lefse)), 2)
})
