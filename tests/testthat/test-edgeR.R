context("edgeR")

test_that("check the norm factors in edgeR", {
    # skip_on_bioc()
    # library(edgeR)
    # normalize, summarize, and then convert to metagenomeSeq
    # ps <- pediatric_ibd
    # ps <- preprocess_ps(ps)
    # ps_normed <- norm_tmm(ps)
    # ps_summarized <- summarize_taxa(ps_normed)
    # dge <- phyloseq2edgeR(ps_summarized)
    # abd1 <- dge$counts
    #
    # # run_edger
    # abd2 <- as(otu_table(mm_dge), "matrix")
    # expect_equal(abd1, abd2)
})

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

    # ps_test <- phyloseq::subset_samples(
    #   cid_ying,
    #   Consistency %in% c("formed stool", "liquid", "semi-formed")
    # )
    # mm_test <- run_edger(
    #   ps_test,
    #   "Consistency",
    #   method  = "QLFT",
    #   pvalue_cutoff = 0.05,
    #   p_adjust = "fdr"
    # )
    #
    # # two groups comparison for multiple groups
    # mm_test2 <- run_edger(
    #   ps_test,
    #   "Consistency",
    #   contrast = c("liquid", "semi-formed"),
    #   pvalue_cutoff = 0.05,
    #   p_adjust = "fdr"
    # )
})
