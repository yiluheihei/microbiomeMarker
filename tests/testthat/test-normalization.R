context("taxa abundance normalization")

test_that("ensure the results are the same for phylosq, otu_table, data.frame, matrix object ", {
  ot <- otu_table(enterotypes_arumugam)
  df <- as.data.frame(ot)
  mat <- as.matrix(df)

  # tss
  tss_ps <- norm_tss(enterotypes_arumugam)
  tss_ot <- norm_tss(ot)
  tss_df <- normalize(df, "TSS")
  tss_mat <- normalize(mat, "TSS")
  expect_equal(tss_ot, otu_table(tss_ps))
  expect_equal(as.data.frame(tss_ot), tss_df)
  expect_equal(as.matrix(tss_df), tss_mat)



  # rarefy
  rarefy_ps <- norm_rarefy(enterotypes_arumugam, rng_seed = 2020)
  rarefy_ot <- norm_rarefy(ot, rng_seed = 2020)
  rarefy_df <- normalize(df, "rarefy", rng_seed = 2020)
  rarefy_mat <- normalize(mat, "rarefy", rng_seed = 2020)
  expect_equal(otu_table(rarefy_ps), rarefy_ot)
  expect_equal(as.data.frame(rarefy_ot), rarefy_df)
  expect_equal(as.matrix(rarefy_df), rarefy_mat)


  # css
  css_ps <- norm_css(enterotypes_arumugam)
  css_ot <- norm_css(ot)
  attr(css_ot, "metagenomeSeq_norm_factor") <- NULL
  css_df <- normalize(df, "CSS")
  css_mat <- normalize(mat, "CSS")
  expect_identical(otu_table(css_ps), css_ot)
  expect_identical(as.data.frame(css_ot), css_df)
  expect_identical(as.matrix(css_df), css_mat)

  # rle
  rle_ps <- norm_rle(enterotypes_arumugam)
  rle_ot <- norm_rle(ot)
  rle_df <- normalize(df, "RLE")
  rle_mat <- normalize(mat, "RLE")
  expect_equal(otu_table(rle_ps), rle_ot)
  expect_equal(as.data.frame(rle_ot), rle_df)
  expect_equal(as.matrix(rle_df), rle_mat)

  # tmm
  tmm_ps <- norm_tmm(enterotypes_arumugam)
  tmm_ot <- norm_tmm(ot)
  tmm_df <- normalize(df, "TMM")
  tmm_mat <- normalize(mat, "TMM")
  expect_equal(otu_table(tmm_ps), tmm_ot)
  expect_equal(as.data.frame(tmm_ot), tmm_df)
  expect_equal(as.matrix(tmm_df), tmm_mat)

  # clr
  clr_ps <- norm_clr(enterotypes_arumugam)
  clr_ot <- norm_clr(ot)
  clr_df <- normalize(df, "CLR")
  clr_mat <- normalize(mat, "CLR")
  expect_equal(otu_table(clr_ps), clr_ot)
  expect_equal(as.data.frame(clr_ot), clr_df)
  expect_equal(as.matrix(clr_df), clr_mat)


  # value
  value_ps <- norm_value(enterotypes_arumugam, 1e6)
  value_ot <- norm_value(ot, 1e6)
  value_df <- normalize(df, 1e6)
  value_mat <- normalize(mat, 1e6)
  expect_equal(otu_table(value_ps), value_ot)
  expect_equal(as.data.frame(value_ot), value_df)
  expect_equal(as.matrix(value_df), value_mat)
})
