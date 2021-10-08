context("`abundances` methods")

ps <- phyloseq::phyloseq(
    otu_table = otu_table(
        matrix(
            sample(100, 40),
            nrow = 2,
            dimnames = list(
                c("feature1", "feature2"),
                paste0("sample", 1:20)
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
            group = rep(c("group1", "group2"), 10),
            row.names = paste0("sample", 1:20)
        )
    )
)

test_that("abundances", {
    abd <- abundances(ps, "log10p")
    abd_t <- transform_abundances(ps, "log10p")

    expect_identical(abd, as(otu_table(abd_t), "matrix"))
})

test_that("abundances normalization", {
    abd_norm <- normalize(ps, "TSS")
    expect_true(all(colSums(otu_table(abd_norm)) == 1))
    expect_is(abd_norm, "phyloseq")
    expect_is(normalize(otu_table(ps)), "otu_table")
})
