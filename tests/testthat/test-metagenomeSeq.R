context("metagenomeSeq")

test_that("check the norm factors in metagenomeSeq", {
  skip_on_bioc()
  library(metagenomeSeq)
  # normalize, summarize, and then convert to metagenomeSeq
  ps <- pediatric_ibd
  ps_normed <- norm_css(ps)
  ps_summarized <- summarize_taxa(ps_normed)
  mgs <- phyloseq2metagenomeSeq(ps_summarized)

  # convert to metagenomeSeq, convert to metagenomeSeq and setting the norm
  # factors
  mgs2 <- phyloseq2metagenomeSeq(ps)
  nf <- metagenomeSeq::calcNormFactors(mgs2)
  ps_summarized2 <- summarize_taxa(ps)
  mgs_summarized <- phyloseq2metagenomeSeq(ps_summarized2)
  pData(mgs_summarized@expSummary$expSummary)$normFactors <- nf

  abd1 <- MRcounts(mgs)
  abd2 <- MRcounts(mgs_summarized, norm = TRUE)
  expect_equal(abd1, abd2)
})

test_that("result of metagenomeSeq", {
  res <- run_deseq2(
    pediatric_ibd,
    "Class",
    "Control", "CD",
    p_adjust = "fdr"
  )

  expect_output_file(
    round_DF(marker_table(res)),
    test_path("out/test-metagenomeSeq.txt"),
    print = TRUE
  )
})
