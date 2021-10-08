context("assignment methods")

marker <- marker_table(
    data.frame(
        feature = paste0("sp", 1:5),
        enrich_group = c("cr", "er", "cr", "cr", "er"),
        stringsAsFactors = FALSE
    )
)
otu <- otu_table(
    data.frame(
        s1 = runif(10),
        s2 = runif(10)
    ),
    taxa_are_rows = TRUE
)
tax <- tax_table(data.frame(feature = paste0("sp", 1:10)) %>% as.matrix())

test_mm <- microbiomeMarker(
    marker = marker,
    otu_table = otu,
    tax_table = tax
)

test_that("otu_table<- works correctly", {
    expect_silent(otu_table(test_mm) <- otu_table(test_mm))
    expect_silent(otu_table(test_mm) <- test_mm)
    ps_assign <- phyloseq(otu_table(test_mm), tax_table(test_mm))
    expect_silent(otu_table(test_mm) <- ps_assign)
})

test_that("marker_table<- works correcly", {
    new_marker <- marker[1, ]
    expect_silent(marker_table(test_mm) <- new_marker)
})
