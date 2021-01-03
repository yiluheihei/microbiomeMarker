context("test multiple group enterotype test")

tukey_res <- posthoc_test(enterotype, "Enterotype", method = "tukey")

# round_DF <- function(DF) {
#   round2 <- function(x) {
#     ifelse(
#       x <= 1e-5,
#       as.numeric(formatC(x, format = "g", digits = 5)),
#       as.numeric(formatC(x, format = "f", digits = 5))
#     )
#   }
#   purrr::map_if(as.data.frame(DF), is.numeric, round2) %>%
#     dplyr::bind_cols() %>%
#     as.data.frame()
# }

test_that("etaseq effect size", {
  etasq <- calc_etasq(c(1, 2, 1.2, 3, 4, 1.4), c("a", "b", "c", "a", "b", "c"))
  expect_equal(signif(etasq, 3), 0.421)
})

test_that("test multiple group enterotype result", {
  skip_on_cran()
  skip_on_bioc()

  # error group
  expect_error(
    test_multiple_groups(enterotype, "Entertype"),
    regexp = "`group` must in the field of sample meta data",
    fixed = TRUE
  )

  expect_known_output(
    round_DF(marker_table(mm_anova)),
    test_path("out/test-multiple-group-anova.txt"),
    print = TRUE
  )

  expect_known_output(
    round_DF(marker_table(mm_kruskal)),
    test_path("out/test-multiple-group-kruk.txt"),
    print = TRUE
  )

})

test_that("test post hoc test result", {
  skip_on_cran()
  skip_on_bioc()

  expect_known_output(
   round_DF(tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]]),
    test_path("out/test-post-hoc-tukey.txt"),
    print = TRUE
  )

  games_res <- posthoc_test(enterotype, "Enterotype", method = "games_howell")
  expect_known_output(
    round_DF(games_res@result[["p__Bacteroidetes|g__Bacteroides"]]),
    test_path("out/test-post-hoc-games.txt"),
    print = TRUE
  )

  scheffe_res <- posthoc_test(enterotype, "Enterotype", method = "scheffe")
  expect_known_output(
    round_DF(scheffe_res@result[["p__Bacteroidetes|g__Bacteroides"]]),
    test_path("out/test-post-hoc-scheffe.txt"),
    print = TRUE
  )

  welch_res <- posthoc_test(enterotype, "Enterotype" , method = "welch_uncorrected")
  expect_known_output(
    round_DF(welch_res@result[["p__Bacteroidetes|g__Bacteroides"]]),
    test_path("out/test-post-hoc-welch.txt"),
    print = TRUE
  )
})

test_that("test visualization of post hoc test, p value significance level ", {
  expect_equal(
    pvalue2siglevel(c(0.05, 0.01, 0.0001, 0.06)),
    c("*", "**", "***", "NS.")
  )
})

test_that("test visualization of post hoc test, data of signicance level annotation", {
  # single feature
  abd <- tukey_res@abundance
  group <- abd$group
  pht_df <- as.data.frame(tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]])
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
  annotation_all$y_position <- formatC(annotation_all$y_position, format = "g", digits = 5)
  expect_known_output(
    head(annotation_all),
    test_path("out/test-posthoc-vis-sig_annotation.txt"),
    print = TRUE
  )
})
