context("nmarker method")

test_that("nmarker method", {
  expect_equal(nmarker(mm_lefse), 20)
  expect_equal(nmarker(marker_table(mm_lefse)), 20)
})
