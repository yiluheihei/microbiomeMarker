context("test utilities")

test_that("check upper first letter", {
  expect_equal(
    upper_firstletter(c("abc", "ABC", "Abc")),
    c("Abc", "Abc", "Abc")
  )
})
