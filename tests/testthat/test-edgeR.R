context("edgeR")

test_that("result of edger", {
    data(pediatric_ibd)
    mm_edger <- run_edger(
        pediatric_ibd,
        "Class",
        pvalue_cutoff = 0.1,
        p_adjust = "fdr"
    )
    expect_output_file(
        print(marker_table(mm_edger), digits = 5),
        test_path("out/test-edger.txt"),
        print = TRUE
    )
})
