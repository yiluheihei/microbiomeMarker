context("test cladogram")

test_that("Generate unique id for short annotation label", {
  uid <- get_unique_id(500)
  expect_equal(uid[26], "z")
  expect_equal(uid[30], "ad")
  expect_equal(uid[260], "iz")
  expect_equal(uid[500], "sf")
})

test_that("drop the levels of groups (no marker) if the enrich_group is a factor", {
  marker <- readRDS("data/generate_cladogram_annotation.rds")
  group_n <- length(unique(marker$enrich_group))
  expect_error(
    generate_cladogram_annotation(marker, color = colors()[1:group_n]),
    NA
  )
})
