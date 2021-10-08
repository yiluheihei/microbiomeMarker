test_that("supervised machine learning method workds properly", {
    data(enterotypes_arumugam)
    ps_small <- phyloseq::subset_taxa(
        enterotypes_arumugam,
        Phylum %in% c("Firmicutes", "Bacteroidetes")
    )
    set.seed(2021)
    mm_lr <- run_sl(
        ps_small,
        group = "Gender",
        taxa_rank = "Genus",
        nfolds = 2,
        nrepeats = 1,
        top_n = 15,
        norm = "TSS",
        method = "LR",
    )
    expect_identical(nmarker(mm_lr), 15L)
    expect_identical(mm_lr@norm_method, "TSS")
    expect_identical(mm_lr@diff_method, "logistic regression")
})
