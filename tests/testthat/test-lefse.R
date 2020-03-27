context("test lefse output")

test_that("lefse output of oxygen", {
  skip_on_cran()
  lefse_out <- lefse(
    oxygen,
    normalization = 1e6,
    bootstrap_n = 5,
    summarize = "lefse",
    class = "oxygen_availability",
    subclass = "body_site"
  )

  expect_known_output(
    lefse_out,
    test_path("out/test-lefse-out.txt"),
    print = TRUE
  )
})

