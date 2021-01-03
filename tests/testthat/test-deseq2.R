context("deseq2 algorithm")

test_that("deseq2 algorithm", {
  skip_on_bioc()

  expect_output_file(
    round_DF(marker_table(mm_des)),
    test_path("out/test-deseq2.txt"),
    print = TRUE
  )
})
