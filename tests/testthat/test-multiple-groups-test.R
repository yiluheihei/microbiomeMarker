context("test multiple groups test")

test_that("etaseq effect size", {
  etasq <- calc_etasq(c(1, 2, 1.2, 3, 4, 1.4), c("a", "b", "c", "a", "b", "c"))
  expect_equal(signif(etasq, 3), 0.421)
})

test_that("test multiple groups result", {
  skip_on_cran()
  skip_on_bioc()

  ps <- phyloseq::subset_samples(
    enterotypes_arumugam,
    Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
  )

  # error group
  expect_error(
    test_multiple_groups(ps, "Entertype"),
    regexp = "`group` must in the field of sample meta data",
    fixed = TRUE
  )

  res_anova <- test_multiple_groups(ps, "Enterotype", effect_size_cutoff = 0.7)
  expect_known_output(
    res_anova,
    test_path("out/test-multiple-group-anova.txt"),
    print = TRUE
  )

  res_kruk <- test_multiple_groups(ps, "Enterotype", method = "kruskal")
  expect_known_output(
    res_kruk,
    test_path("out/test-multiple-group-kruk.txt"),
    print = TRUE
  )

})
