context("bar plot of effect size")

test_that("feature label in bar plot", {
  feature <- "Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae"
  short_feature <- "Bacteria|Bacteroidetes"

  expect_equal(get_feature_label(feature, 1), "Bacteroidaceae")
  expect_equal(get_feature_label(feature, 2), "Bacteroidales|Bacteroidaceae")
  expect_equal(
    get_feature_label(feature, 7),
    "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
  )
  expect_equal(
    get_feature_label(feature, 0),
    "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
  )
  expect_equal(get_feature_label(short_feature, 0), "Bacteria|Bacteroidetes")

  # replace "Unknown" in the species level as "sp."
  expect_equal(
    replace_unknown_species("s__Leucobacter_Unknown"),
    "s__Leucobacter_sp."
  )
  expect_equal(
    replace_unknown_species("g__abcd|s__Leucobacter_Unknown"),
    "g__abcd|s__Leucobacter_sp."
  )
  expect_equal(
    replace_unknown_species("g__abcd"),  "g__abcd"
  )
  expect_equal(
    replace_unknown_species("g__abcd|s__Leucobacter_sp."),
    "g__abcd|s__Leucobacter_sp."
  )
})

p_lda <- plot_ef_bar(mm_lefse)

test_that("label of x", {
  expect_identical(p_lda$labels$x, "LDA score (log10)")

  p_diff_mean <- plot_ef_bar(mm_welch)
  expect_identical(p_diff_mean$labels$x, "Differential means")

  p_eta_squared <- plot_ef_bar(mm_anova)
  expect_identical(p_eta_squared$labels$x, "Eta squared")

  p_logfc <- plot_ef_bar(mm_des)
  expect_identical(p_logfc$labels$x, "log2 Fold Change")
})

test_that("effect size in descending order", {
  ef <- p_lda$data$effect_size
  ef1 <- ef[p_lda$data$enrich_group == "Healthy"]
  ef2 <- ef[p_lda$data$enrich_group == "Tumor"]
  expect_true(all(diff(ef1) >= 0))
  expect_true(all(diff(ef2) >= 0))

})
