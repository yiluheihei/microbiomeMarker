test_that("limma voom", {
    data(enterotypes_arumugam)
    enterotype <- phyloseq::subset_samples(
        enterotypes_arumugam,
        Enterotype %in% c("Enterotype 1", "Enterotype 2", "Enterotype 3")
    )
    mm_lv <- run_limma_voom(
        enterotype,
        "Enterotype",
        contrast = c("Enterotype 3", "Enterotype 2"),
        pvalue_cutoff = 0.05,
        p_adjust = "fdr"
    )
    
    expect_snapshot(print(marker_table(mm_lv), digits = 5))
})
