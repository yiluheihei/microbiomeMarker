context("import qiime2 output into phyloseq object")

test_that("whether the row.names of feature table is dna sequence or not", {
    expect_false(is_dna_seq("3597a2689efaf5525ce460494a8ac383"))
    expect_true(is_dna_seq(paste0(
        "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTG",
        "GTGATTTAAGTCAGCGGTGAAAGTTTGTGGCTCAACCATAAAATTGCCGTTGAA",
        "ACTGGGTTACTTGAGTGTGTTTGAGGTAGGCGGAATGCGTGG"
    )))
})
