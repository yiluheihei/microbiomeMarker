context("tax summarize")

ps <- phyloseq::phyloseq(
    otu_table = otu_table(
        matrix(
            sample(100, 40),
            nrow = 4,
            dimnames = list(
                paste0("otu", 1:4),
                paste0("sample", 1:10)
            )
        ),
        taxa_are_rows = TRUE
    ),
    tax_table = tax_table(
        matrix(
            c(rep("g1", 4), rep(c("s1", "s2"), 2)),
            nrow = 4, byrow = FALSE,
            dimnames = list(paste0("otu", 1:4), c("Genus", "Species"))
        )
    ),
    sam_data = sample_data(
        data.frame(
            group = rep(c("group1", "group2"), 5),
            row.names = paste0("sample", 1:10)
        )
    )
)
ps_summarized <- summarize_taxa(ps)

test_that("check whether phyloseq tax summarized or not", {
    expect_false(check_tax_summarize(ps))
    expect_true(check_tax_summarize(ps_summarized))
})

test_that("check the summarize_taxa", {
    sep <- "|"
    taxa <- otu_table(ps_summarized)

    expect_true(any(grepl(sep, row.names(taxa), fixed = TRUE)))
    expect_equal(nrow(taxa), 3L)

    expect_error(
        summarize_taxa(ps, level = "abc"),
        "`level` must in the ranks"
    )
})

test_that("extract prefix of names of taxonomic ranks", {
    expect_identical(extract_prefix(ps), c("g", "s"))
    expect_identical(
        extract_prefix(rank_names(ps)),
        c("g", "s")
    )
})
