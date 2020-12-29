context("metagenomeSeq")

test_that("check the norm factors in metagenomeSeq", {
  skip_on_bioc()
  # library(metagenomeSeq)
  # # normalize, summarize, and then convert to metagenomeSeq
  # ps <- pediatric_ibd
  # ps_normed <- norm_css(ps)
  # ps_summarized <- summarize_taxa(ps_normed)
  # mgs <- phyloseq2metagenomeSeq(ps_summarized)
  #
  # # convert to metagenomeSeq, convert to metagenomeSeq and setting the norm
  # # factors
  # mgs2 <- phyloseq2metagenomeSeq(ps)
  # nf <- metagenomeSeq::calcNormFactors(mgs2)
  # ps_summarized2 <- summarize_taxa(ps)
  # mgs_summarized <- phyloseq2metagenomeSeq(ps_summarized2)
  # pData(mgs_summarized@expSummary$expSummary)$normFactors <- nf
  #
  # abd1 <- MRcounts(mgs)
  # abd2 <- MRcounts(mgs_summarized, norm = TRUE)
  # expect_equal(abd1, abd2)
})

test_that("result of metagenomeSeq", {
  mm_mgs <- run_metagenomeseq(
    pediatric_ibd,
    norm = "CSS",
    "Class",
    "Control",
    "CD",
    p_value_cutoff = 0.1,
    p_adjust = "fdr"
  )

  expect_output_file(
    round_DF(marker_table(mm_mgs)),
    test_path("out/test-metagenomeSeq.txt"),
    print = TRUE
  )
})
