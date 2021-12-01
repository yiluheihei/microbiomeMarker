test_that("import_picrust2 works", {
    sam_tab <- system.file(
        "extdata", "picrust2_metadata.tsv",
        package = "microbiomeMarker")
    feature_tab <- system.file(
        "extdata", "path_abun_unstrat_descrip.tsv.gz",
        package = "microbiomeMarker")
    ps <- import_picrust2(feature_tab, sam_tab, trait = "PATHWAY")
    
    expect_identical(rank_names(ps), c("Picrust_trait", "Picrust_description"))
})
