context("test multiple group enterotype test")

data(enterotypes_arumugam)
enterotype <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)
# mm_anova <- run_test_multiple_groups(
#   enterotype,
#   group = "Enterotype",
#   method = "anova",
#   effect_size_cutoff = 0.7
# )

tukey_res <- run_posthoc_test(enterotype, "Enterotype", method = "tukey")

test_that("etaseq effect size", {
  etasq <- calc_etasq(c(1, 2, 1.2, 3, 4, 1.4), c("a", "b", "c", "a", "b", "c"))
  expect_equal(signif(etasq, 3), 0.421)
})

test_that("test multiple group enterotype result", {

  # error group
  expect_error(
    run_test_multiple_groups(enterotype, "Entertype"),
    regexp = "`group` must in the field of sample meta data",
    fixed = TRUE
  )

  # expect_known_output(
  #   print(marker_table(mm_anova), digits = 5),
  #   test_path("out/test-multiple-group-anova.txt"),
  #   print = TRUE
  # )

})

test_that("test post hoc test result", {
  expect_known_output(
    print(
      tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]],
      digits = 5),
    test_path("out/posthoc-turkey.txt"),
    print = TRUE
  )
})

test_that("test visualization of post hoc test, p value significance level ", {
  expect_equal(
    pvalue2siglevel(c(0.05, 0.01, 0.0001, 0.06)),
    c("*", "**", "***", "NS.")
  )
})

test_that(
  "test visualization of posthoc test, data of signicance level annotation", {
  # single feature
  abd <- tukey_res@abundance
  group <- abd$group
  pht_df <- data.frame(tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]])
  annotation_single <- get_sig_annotation_single(
    abd[["p__Bacteroidetes|g__Bacteroides"]],
    pht_df,
    group
  )
  annotation_single$y_position <- formatC(
    annotation_single$y_position,
    format = "g",
    digits = 5
  )
  expect_known_output(
    annotation_single,
    test_path("out/test-posthoc-vis-sig_annotation_single.txt"),
    print = TRUE
  )

  # all features
  annotation_all <- get_sig_annotation(tukey_res)
  annotation_all$y_position <- formatC(
    annotation_all$y_position,
    format = "g",
    digits = 5
  )
  expect_known_output(
    head(annotation_all),
    test_path("out/test-posthoc-vis-sig_annotation.txt"),
    print = TRUE
  )
})
