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

test_that("taxa prefix",
  expect_identical(get_prefix(c("Phylum", "Genus")), c("p__", "g__"))
)

test_that("var in sample_data", {
  expect_error(
    check_var_in_meta("abc", sample_data(ecam)),
    "variable of `sample_meta`"
  )
  expect_silent(check_var_in_meta("delivery", sample_data(ecam)))
})
