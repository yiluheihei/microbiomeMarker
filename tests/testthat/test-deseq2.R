context("deseq2 algorithm")

test_that("deseq2 algorithm", {
  skip_on_bioc()

  expect_output_file(
    round_DF(marker_table(mm_des)),
    test_path("out/test-deseq2.txt"),
    print = TRUE
  )

  expect_error(
    run_deseq2(
      pediatric_ibd,
      "Class",
      pvalue_cutoff = 0.05,
      p_adjust = "fdr"
    ),
    "`contrast` is requried"
  )

  # # test for multiple groups comparison
  # ps_test <- phyloseq::subset_samples(
  #   cid_ying,
  #   Consistency %in% c("formed stool", "liquid", "semi-formed")
  # )
  # mm_test <- run_deseq2(
  #   ps_test,
  #   "Consistency",
  #   pvalue_cutoff = 0.05,
  #   p_adjust = "fdr"
  # )
  #
  # # two groups comparison for multiple groups
  # mm_test2 <- run_deseq2(
  #   ps_test,
  #   "Consistency",
  #   contrast = c("liquid", "semi-formed"),
  #   pvalue_cutoff = 0.05,
  #   p_adjust = "fdr"
  # )
})
