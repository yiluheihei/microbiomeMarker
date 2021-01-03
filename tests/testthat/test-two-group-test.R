context("test two group tests")

test_that("ration", {
  abd1 <- rep(0, 6)
  abd2 <- rep(0, 6)
  expect_equal(calc_ratio(abd1, abd2), 0)

  abd1 <- c(0.00014940, 0.05774625, 0.05419578, 0.05414964, 0.02028051)
  abd2 <- rep(0, 5)
  expect_equal(round(calc_ratio(abd1, abd2), 5), 1.00278)

  abd2 <- rev(abd1)
  expect_equal(round(calc_ratio(abd1, abd2), 5), 1)
})

test_that("test two group result", {
  skip_on_cran()

  expect_known_output(
    mm_welch,
    test_path("out/test-two-group-test-welch.txt"),
    print = TRUE
  )

  expect_known_output(
    mm_t,
    test_path("out/test-two-group-test-t.txt"),
    print = TRUE
  )

  expect_known_output(
    mm_white,
    test_path("out/test-two-group-test-white.txt"),
    print = TRUE
  )

})
