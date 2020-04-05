context("nclass method")

test_that("nclass method", {
  expect_equal(nclass(lefse_out), 3)
  expect_equal(nclass(marker_table(lefse_out)), 3)
})
