test_that("result of edger", {
    data(pediatric_ibd)
    mm_edger <- run_edger(
        pediatric_ibd,
        "Class",
        pvalue_cutoff = 0.1,
        p_adjust = "fdr"
    )
    expect_snapshot(print(marker_table(mm_edger), digits = 5))
})
