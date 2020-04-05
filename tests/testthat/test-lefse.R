context("test lefse output")

test_that("lefse output of oxygen", {
  skip_on_cran()

  expect_known_output(
    lefse_out,
    test_path("out/test-lefse-out.txt"),
    print = TRUE
  )
})
