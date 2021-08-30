context("test lefse output")

set.seed(2020)
# lefse - lda
data(kostic_crc)
mm_lefse <- run_lefse(
  kostic_crc,
  wilcoxon_cutoff = 0.01,
  group = "DIAGNOSIS",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4,
)

test_that("lefse output of oxygen", {
  expect_known_output(
    mm_lefse,
    test_path("out/test-lefse-out.txt"),
    print = TRUE
  )

  expect_known_output(
    marker_table(mm_lefse),
    test_path("out/test-lefse-out_marker.txt"),
    print = TRUE
  )

})
