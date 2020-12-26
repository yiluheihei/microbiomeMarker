context("deseq2 algorithm")

test_that("deseq2 algorithm", {
  res <- run_deseq2(
    oxygen,
    "oxygen_availability",
    "High_O2", "Low_O2",
    p_adjust = "none"
  )

  expect_output_file(
    round_DF(marker_table(res)),
    test_path("out/test-deseq2.txt"),
    print = TRUE
  )
})
