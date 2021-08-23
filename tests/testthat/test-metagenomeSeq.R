context("metagenomeSeq")

test_that("check the norm factors in metagenomeSeq", {
  # skip_on_bioc()
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
  ps <- phyloseq::subset_samples(
    cid_ying,
    Consistency %in% c("formed stool", "liquid", "semi-formed")
  )

  expect_output_file(
    round_DF(marker_table(mm_mgs)),
    test_path("out/test-metagenomeSeq.txt"),
    print = TRUE
  )

  # expect_error(
  #   run_metagenomeseq(pediatric_ibd, "Class"),
  #   "`contrast` is required"
  # )

  expect_error(
    run_metagenomeseq(ps, "Consistency"),
    "ZILN method do not allows"
  )

  expect_error(
    run_metagenomeseq(ps, "Consistency", contrast = c("Control", "CD")),
    "ZILN method do not allows"
  )

  # run_metagenomeseq(
  #   ps,
  #   "Consistency",
  #   contrast = c("Control", "CD"),
  #   p_adjust = "fdr",
  #   method = "ZIG"
  # )

  # run_metagenomeseq(ps, "Consistency", method = "ZIG")
})

test_that("get enrich group of a featrue of multiple groups comparison", {
  group_pairs <- list(
    c("a", "b"),
    c("a", "c"),
    c("b", "c")
  )
  logFC_pairs <- c(1, 1, 1)

  expect_identical(get_mgs_enrich_group(group_pairs, logFC_pairs), "a")
})
