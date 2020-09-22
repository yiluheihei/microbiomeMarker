context("test cladogram")

test_that("Generate unique id for short annotation label", {
  uid <- get_unique_id(500)
  expect_equal(uid[26], "z")
  expect_equal(uid[30], "ad")
  expect_equal(uid[260], "iz")
  expect_equal(uid[500], "sf")
})
