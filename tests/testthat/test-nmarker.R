context("nmarker method")

test_that("nmarker method", {
  expect_equal(nmarker(lefse_out), 51)
  expect_equal(nmarker(marker_table(lefse_out)), 51)
})
