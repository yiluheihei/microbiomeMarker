context("metagenomeSeq")

test_that("result of metagenomeSeq", {
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
