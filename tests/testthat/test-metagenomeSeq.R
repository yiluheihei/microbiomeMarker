context("metagenomeSeq")

test_that("check the norm factors in metagenomeSeq", {
    # skip_on_bioc()
    # library(metagenomeSeq)
    # # normalize, summarize, and then convert to metagenomeSeq
    # ps <- pediatric_ibd
    # ps_normed <- norm_css(ps)
    # ps_summarized <- summarize_taxa(ps_normed)
    # mgs <- phyloseq2metagenomeSeq(ps_summarized)
    #
    # # convert to metagenomeSeq, convert to metagenomeSeq and setting the norm
    # # factors
    # mgs2 <- phyloseq2metagenomeSeq(ps)
    # nf <- metagenomeSeq::calcNormFactors(mgs2)
    # ps_summarized2 <- summarize_taxa(ps)
    # mgs_summarized <- phyloseq2metagenomeSeq(ps_summarized2)
    # pData(mgs_summarized@expSummary$expSummary)$normFactors <- nf
    #
    # abd1 <- MRcounts(mgs)
    # abd2 <- MRcounts(mgs_summarized, norm = TRUE)
    # expect_equal(abd1, abd2)
})

test_that("result of metagenomeSeq", {
    # data(pediatric_ibd)
    # mm_mgs <- run_metagenomeseq(
    #   pediatric_ibd,
    #   "Class",
    #   # contrast = c("CD","Control"),
    #   pvalue_cutoff = 0.1,
    #   p_adjust = "fdr"
    # )
    #
    # expect_output_file(
    #   print(marker_table(mm_mgs), digits = 5),
    #   test_path("out/test-metagenomeSeq.txt"),
    #   print = TRUE
    # )

    ps <- phyloseq::phyloseq(
        otu_table = otu_table(
            matrix(
                1:12,
                nrow = 2,
                dimnames = list(
                    c("feature1", "feature2"),
                    paste0("sample", 1:6)
                )
            ),
            taxa_are_rows = TRUE
        ),
        tax_table = tax_table(
            matrix(
                c("taxa1", "taxa2"),
                nrow = 2,
                dimnames = list(c("feature1", "feature2"), c("Species"))
            )
        ),
        sam_data = sample_data(
            data.frame(
                group = rep(c("group1", "group2", "group3"), 2),
                row.names = paste0("sample", 1:6)
            )
        )
    )
    expect_error(
        run_metagenomeseq(ps, "group"),
        "ZILN method do not allows"
    )
    expect_error(
        run_metagenomeseq(ps, "group", contrast = c("group1", "group2")),
        "ZILN method do not allows"
    )
})

test_that("get enrich group of a featrue of multiple groups comparison", {
    group_pairs <- list(
        c("a", "b"),
        c("a", "c"),
        c("b", "c")
    )
    logFC_pairs <- c(1, 1, 1)

    expect_identical(get_mgs_enrich_group(group_pairs, logFC_pairs), "a")
})
