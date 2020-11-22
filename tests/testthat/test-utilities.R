context("test utilities")

test_that("check upper first letter", {
  expect_equal(
    upper_firstletter(c("abc", "ABC", "Abc")),
    c("Abc", "Abc", "Abc")
  )
})

test_that(" check whether all names of taxonomic ranks include in available_ranks", {
  expect_true(check_rank_names(oxygen))
  expect_true(check_rank_names(caporaso))
})
