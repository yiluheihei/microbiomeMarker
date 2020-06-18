context("test multiple groups test")

ps <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)

round_DF <- function(DF) {
  round2 <- function(x) {
    ifelse(
      x <= 1e-5,
      as.numeric(formatC(x, format = "g", digits = 5)),
      as.numeric(formatC(x, format = "f", digits = 5))
    )
  }
  purrr::map_if(as.data.frame(DF), is.numeric, round2) %>%
    dplyr::bind_cols() %>%
    as.data.frame()
}

test_that("etaseq effect size", {
  etasq <- calc_etasq(c(1, 2, 1.2, 3, 4, 1.4), c("a", "b", "c", "a", "b", "c"))
  expect_equal(signif(etasq, 3), 0.421)
})

test_that("test multiple groups result", {
  skip_on_cran()
  skip_on_bioc()

  # error group
  expect_error(
    test_multiple_groups(ps, "Entertype", rank_name = "Genus"),
    regexp = "`group` must in the field of sample meta data",
    fixed = TRUE
  )

  res_anova <- test_multiple_groups(
    ps,
    "Enterotype",
    rank_name = "Genus",
    effect_size_cutoff = 0.7
  )
  expect_known_output(
    round_DF(marker_table(res_anova)),
    test_path("out/test-multiple-group-anova.txt"),
    print = TRUE
  )

  res_kruk <- test_multiple_groups(
    ps,
    "Enterotype",
    rank_name = "Genus",
    method = "kruskal"
  )
  expect_known_output(
    round_DF(marker_table(res_kruk)),
    test_path("out/test-multiple-group-kruk.txt"),
    print = TRUE
  )

})

test_that("test post hoc test result", {
  skip_on_cran()
  skip_on_bioc()

  tukey_res <- posthoc_test(ps, "Enterotype", rank_name = "Genus", method = "tukey")
  expect_known_output(
   round_DF(tukey_res@result[[1]]),
    test_path("out/test-post-hoc-tukey.txt"),
    print = TRUE
  )

  games_res <- posthoc_test(ps, "Enterotype", rank_name = "Genus", method = "games_howell")
  expect_known_output(
    round_DF(games_res@result[[1]]),
    test_path("out/test-post-hoc-games.txt"),
    print = TRUE
  )

  scheffe_res <- posthoc_test(ps, "Enterotype", rank_name = "Genus", method = "scheffe")
  expect_known_output(
    round_DF(scheffe_res@result[[1]]),
    test_path("out/test-post-hoc-scheffe.txt"),
    print = TRUE
  )

  welch_res <- posthoc_test(ps, "Enterotype", rank_name = "Genus", method = "welch_uncorrected")
  expect_known_output(
    round_DF(welch_res@result[[1]]),
    test_path("out/test-post-hoc-welch.txt"),
    print = TRUE
  )
})
