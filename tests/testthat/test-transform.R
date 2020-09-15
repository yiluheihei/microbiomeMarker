context("test abundance transformation")

ps_t <- transform_abundances(oxygen)
otutable_t <- transform_abundances(otu_table(oxygen))
microbiomeMarker_t <- transform_abundances(lefse_out)

test_that("return a object: the same class with the argument `object`", {
  expect_true(inherits(ps_t, "phyloseq"))
  expect_true(inherits(otutable_t, "otu_table"))
  expect_true(inherits(microbiomeMarker_t, "microbiomeMarker"))
})

test_that("transformation", {
  # transformed using log10(1 + x) if data contains zero
  t_log10 <- transform_abundances(enterotypes_arumugam, "log10")
  t_log10p <- transform_abundances(enterotypes_arumugam, "log10p")
  expect_warning(transform_abundances(enterotypes_arumugam, "log10"), "contains zeroes")
  expect_identical(t_log10, t_log10p)
})

