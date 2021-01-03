context("test lefse output")

test_that("lefse output of oxygen", {
  skip_on_cran()

  expect_known_output(
    mm_lefse,
    test_path("out/test-lefse-out.txt"),
    print = TRUE
  )

  expect_known_output(
    marker_table(mm_lefse),
    test_path("out/test-lefse-out_marker.txt"),
    print = TRUE
  )

})
