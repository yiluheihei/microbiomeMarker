test_that("limma voom", {
  mm_lv <- run_limma_voom(
    enterotype,
    "Enterotype",
    contrast = c("Enterotype 3", "Enterotype 2"),
    pvalue_cutoff = 0.05,
    p_adjust = "fdr"
  )

  expect_output_file(
    round_DF(marker_table(mm_lv)),
    test_path("out/test-lv.txt"),
    print = TRUE
  )

})
