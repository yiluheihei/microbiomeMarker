context("lefse bar visualization")

test_that("feature label in bar plot", {
  feature <- "Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae"
  short_feature <- "Bacteria|Bacteroidetes"

  expect_equal(get_feature_label(feature, 1), "Bacteroidaceae")
  expect_equal(get_feature_label(feature, 2), "Bacteroidales|Bacteroidaceae")
  expect_equal(
    get_feature_label(feature, 7),
    "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
  )
  expect_equal(
    get_feature_label(feature, 0),
    "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
  )
  expect_equal(get_feature_label(short_feature, 0), "Bacteria|Bacteroidetes")
})

