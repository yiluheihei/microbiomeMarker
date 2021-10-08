context("test abundance transformation")

data(enterotypes_arumugam)
ps_t <- transform_abundances(enterotypes_arumugam, "log10p")
otutable_t <- transform_abundances(otu_table(enterotypes_arumugam), "log10p")
mm <- microbiomeMarker(
    marker_table = marker_table(data.frame(
        feature = c("speciesA", "speciesB"),
        enrich_group = c("groupA", "groupB"),
        ef_logFC = c(-2, 2),
        pvalue = c(0.01, 0.01),
        padj = c(0.01, 0.01),
        row.names = c("marker1", "marker2")
    )),
    otu_table = otu_table(matrix(
        c(4, 1, 1, 4),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("speciesA", "speciesB"), c("sample1", "sample2"))
    ),
    taxa_are_rows = TRUE
    ),
    tax_table = tax_table(matrix(
        c("speciesA", "speciesB"),
        nrow = 2,
        dimnames = list(c("speciesA", "speciesB"), "Species")
    ))
)
mm_transed <- transform_abundances(mm, "log10p")

test_that("return a object: the same class with the argument `object`", {
    expect_true(inherits(ps_t, "phyloseq"))
    expect_true(inherits(otutable_t, "otu_table"))
    expect_true(inherits(mm_transed, "microbiomeMarker"))

    expect_equal(otu_table(ps_t), otutable_t)
})

test_that("transformation", {
    # transformed using log10(1 + x) if data contains zero
    expect_warning(
        t_log10 <- transform_abundances(enterotypes_arumugam, "log10"),
        "Using log10(1 + x) instead",
        fixed = TRUE
    )
    t_log10p <- transform_abundances(enterotypes_arumugam, "log10p")
    expect_identical(t_log10, t_log10p)
})
