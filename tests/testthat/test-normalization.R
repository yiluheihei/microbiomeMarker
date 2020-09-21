context("taxa abundance normalization")

test_that("ensure the results are the same for phylosq and otu_table object ", {
  # tss
  tss_ps <- norm_tss(enterotypes_arumugam)
  ot <- otu_table(enterotypes_arumugam)
  tss_ot <- norm_tss(ot)
  expect_equal(tss_ot, otu_table(tss_ps))

  # rarefy
  rarefy_ps <- norm_rarefy(enterotypes_arumugam, rng_seed = 2020)
  rarefy_ot <- norm_rarefy(ot, rng_seed = 2020)
  expect_equal(otu_table(rarefy_ps), rarefy_ot)

  # css
  css_ps <- norm_css(enterotypes_arumugam)
  css_ot <- norm_css(ot)
  expect_equal(otu_table(css_ps), css_ot)

  # rle
  rle_ps <- norm_rle(enterotypes_arumugam)
  rle_ot <- norm_rle(ot)
  expect_equal(otu_table(rle_ps), rle_ot)

  # tmm
  tmm_ps <- norm_tmm(enterotypes_arumugam)
  tmm_ot <- norm_tmm(ot)
  expect_equal(otu_table(tmm_ps), tmm_ot)

  # clr
  clr_ps <- norm_clr(enterotypes_arumugam)
  clr_ot <- norm_clr(ot)
  expect_equal(otu_table(clr_ps), clr_ot)

})
