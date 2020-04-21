#' Statistical test between two groups
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param groups character, the variable to set the group
#' @param rank_name character, taxonomic names of [`phyloseq::phyloseq-class`]
#' to compare
#' @param method test method
#' @param p_adjust_method method for multiple test correction
#' @param conf_level numeric, confidence level of interval
#' @importFrom phyloseq rank_names tax_glom
#' @importFrom dplyr select everything
#' @export
#' @author Yang Cao
test_two_groups <- function(ps,
                           groups,
                           rank_name,
                           method = "welch.test",
                           p_adjust_method = "fdr",
                           conf_level = 0.95) {
  stopifnot(inherits(ps, "phyloseq"))
  ranks <- rank_names(ps)
  if (!rank_name %in% ranks) {
    stop("`rank_name` must be one of availabel taxonomic ranks of `ps`")
  }

  # agglomerate tax in the same rank_name
  if (rank_name != ranks[length(ranks)]) {
    ps <- tax_glom(ps, taxrank = rank_name)
  }
  abd <- otu_table(ps) %>%
    t() %>%
    as.data.frame()
  # relative abundance
  abd_sum <- rowSums(abd)
  abd <- sweep(abd, 1, abd_sum, "/")

  sample_meta <- sample_data(ps)
  groups <- sample_meta[[groups]]
  abd_group <- split(abd, groups)

  test_res <- run_welch_test(abd_group, conf_level = conf_level)
  feature <- tax_table(ps)[, rank_name] %>% unclass()
  test_res$feature <- feature[, 1]

  # ratio proportion
  rp <- purrr::pmap_dbl(abd_group, ~ calc_ratio_proportion(.x, .y))
  test_res$ratio_proportion <- rp

  # set the ci and ratio proportion to 0, if both of the mean is 0
  test_res <- mutate(
    test_res,
    ci_lower = ifelse(.data$pvalue == 1, 0, .data$ci_lower),
    ci_upper = ifelse(.data$pvalue == 1, 0, .data$ci_upper)) %>%
    select(.data$feature, everything())

  microbiomeMarker(
    marker_table(test_res),
    otu_table(t(abd), taxa_are_rows = TRUE),
    tax_table(ps)
  )
}

#' welch test
#'
#' @param abd_group a two length list, each element represents the feature
#' abundance of a group
#' @param conf_level numeric, confidence level of the interval, default 0.95
run_welch_test <- function(abd_group, conf_level = 0.95) {
  if (length(abd_group) != 2) {
    stop("welch test requires test between two groups")
  }
  # welch test
  t_res <- purrr::pmap(
    abd_group,
    ~ t.test(.x, .y, conf.level = conf_level)
  )

  # p value
  p <- purrr::map_dbl(t_res, ~ .x$p.value)
  # set the p value to 1 is the result is NA
  p[is.na(p)] <- 1

  # mean proportion of each group
  # different between means
  t_estimate <- purrr::map(t_res, ~ .x$estimate)
  mean_g1 <- purrr::map_dbl(t_estimate, ~ .x[1])
  mean_g2 <- purrr::map_dbl(t_estimate, ~ .x[2])
  diff_means <- mean_g1 - mean_g2

  # confidence interval
  ci <- purrr::map(t_res, ~ .x$conf.int)
  ci_lower <- purrr::map_dbl(ci, ~ .x[1])
  ci_upper <- purrr::map_dbl(ci, ~ .x[2])

  group_names <- names(abd_group)
  mean_names <- paste(group_names, "mean_rel_freq", sep = "_")
  res <- data.frame(
    p,
    mean_g1*100,
    mean_g2*100,
    diff_means*100,
    ci_lower*100,
    ci_upper*100
  )
  names(res) <- c("pvalue", mean_names, "diff_mean", "ci_lower", "ci_upper")

  res
}

#' ratio proportion used for effect size
#'
#' @param abd1,abd2 numeric vector, abundance of a given feature of the group1
#' and group2
#' @param pseducount numeric, pseducount for unobserved data
#'
#' @return numeric ratio proportion for a feature
#' @noRd
calc_ratio_proportion <- function(abd1, abd2, pseudocount = 0.5) {
  n1 <- length(abd1)
  n2 <- length(abd2)

  mean_g1 <- sum(abd1)/n1
  mean_g2 <- sum(abd2)/n2

  if (mean_g1 == 0 || mean_g2 == 0) {
    pseudocount <- pseudocount/(mean_g1 + mean_g2)
    mean_g1 <- mean_g1 + pseudocount
    mean_g2 <- mean_g2 + pseudocount
  }

  res <- mean_g1/mean_g2
  if (is.na(res)) {
    res <- 0
  }

  res
}
