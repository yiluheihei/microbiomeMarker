context("test abundance transformation")

ps_t <- transform_abundances(oxygen, "log10")
otutable_t <- transform_abundances(otu_table(oxygen), "log10")
microbiomeMarker_t <- transform_abundances(lefse_out, "log10")

test_that("return a object: the same class with the argument `object`", {
  expect_true(inherits(ps_t, "phyloseq"))
  expect_true(inherits(otutable_t, "otu_table"))
  expect_true(inherits(microbiomeMarker_t, "microbiomeMarker"))

  expect_equal(otu_table(ps_t), otutable_t)
})

test_that("transformation", {
  # transformed using log10(1 + x) if data contains zero
  t_log10 <- transform_abundances(enterotypes_arumugam, "log10")
  t_log10p <- transform_abundances(enterotypes_arumugam, "log10p")
  expect_warning(transform_abundances(enterotypes_arumugam, "log10"), "contains zeroes")
  expect_identical(t_log10, t_log10p)
})

