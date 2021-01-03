context("edgeR")

test_that("check the norm factors in edgeR", {
  skip_on_bioc()
  # library(edgeR)
  # normalize, summarize, and then convert to metagenomeSeq
  # ps <- pediatric_ibd
  # ps <- preprocess_ps(ps)
  # ps_normed <- norm_tmm(ps)
  # ps_summarized <- summarize_taxa(ps_normed)
  # dge <- phyloseq2edgeR(ps_summarized)
  # abd1 <- dge$counts
  #
  # # run_edger
  # abd2 <- as(otu_table(mm_dge), "matrix")
  # expect_equal(abd1, abd2)
})

test_that("result of edger",{
  expect_output_file(
    round_DF(marker_table(mm_edger)),
    test_path("out/test-edger.txt"),
    print = TRUE
  )
})
