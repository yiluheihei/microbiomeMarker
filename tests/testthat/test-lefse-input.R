context("lefse input preprocessing")

test_that("add missing levels: keep abundance is lower than 1", {
  oxygen_feature <- otu_table(oxygen)
  feature <- add_missing_levels(oxygen_feature)

  expect_true(max(feature) <= 1)
  expect_true(min(feature) >= 0)
})

