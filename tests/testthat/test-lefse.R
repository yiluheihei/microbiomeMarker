# lefse - lda
data(kostic_crc)
kostic_crc_small <- phyloseq::subset_taxa(
    kostic_crc,
    Phylum == "Firmicutes"
)
mm_lefse <- withr::with_seed(
    2020,
    run_lefse(kostic_crc_small,
        wilcoxon_cutoff = 0.01,
        group = "DIAGNOSIS",
        kw_cutoff = 0.01,
        multigrp_strat = TRUE,
        lda_cutoff = 4
    )
)

test_that("lefse output of oxygen", {
    expect_snapshot(mm_lefse)
    expect_snapshot(marker_table(mm_lefse))
})

test_that("create phyloseq object from microbiomeMarker object", {
    ps <- create_ps_from_mm(mm_lefse)
    expect_true(all(marker_table(mm_lefse)$feature %in% taxa_names(ps)))
    
    ps2 <- create_ps_from_mm(mm_lefse, only_marker = FALSE)
    expect_identical(taxa_names(ps2), taxa_names(mm_lefse))
})