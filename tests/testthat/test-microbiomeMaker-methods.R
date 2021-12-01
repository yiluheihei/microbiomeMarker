test_that("nmarker method", {
    mm <- microbiomeMarker(
        marker_table = marker_table(data.frame(
            feature = c("speciesA", "speciesB"),
            enrich_group = c("groupA", "groupB"),
            ef_logFC = c(-2, 2),
            pvalue = c(0.01, 0.01),
            padj = c(0.01, 0.01),
            row.names = c("marker1", "marker2")
        )),
        norm_method = "TSS",
        diff_method = "DESeq2",
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
        )),
        sam_data = sample_data(data.frame(
            group = c("groupA", "groupB"),
            row.names = c("sample1", "sample2")
        ))
    )
    expect_identical(nmarker(mm), 2L)
    expect_identical(nmarker(marker_table(mm)), 2L)
    
    marker_table(mm) <- NULL
    expect_identical(nmarker(mm), 0L)
    
})
