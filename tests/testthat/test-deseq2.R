context("deseq2 algorithm")

test_that("deseq2 algorithm", {
  skip_on_bioc()
  mm_des <- run_deseq2(
    pediatric_ibd,
    "Class",
    "Control",
    "CD",
    p_value_cutoff = 0.05,
    p_adjust = "fdr"
  )

  expect_output_file(
    round_DF(marker_table(mm_des)),
    test_path("out/test-deseq2.txt"),
    print = TRUE
  )
})
