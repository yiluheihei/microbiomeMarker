#' Statistical test for multiple groups
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
#' @param method test method, must be one of "anova" or "kruskal"
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param p_value_cutoff numeric, p value cutoff, default 0.05
#' @param effect_size_cutoff numeric, cutoff of effect size default `NULL`
#' which means no effect size filter
#' @importFrom dplyr mutate bind_cols filter
#' @importFrom stats p.adjust
test_multiple_groups <- function(ps,
                                 group,
                                 method = c("anova", "kruskal"),
                                 p_adjust = c("none", "fdr", "bonferroni", "holm",
                                              "hochberg", "hommel", "BH", "BY"),
                                 p_value_cutoff = 0.05,
                                 effect_size_cutoff = NULL) {
  stopifnot(inherits(ps, "phyloseq"))

  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
  )
  method <- match.arg(method, c("anova", "kruskal"))

  # relative abundance
  abd <- otu_table(ps) %>%
    t() %>%
    as.data.frame()
  abd_sum <- rowSums(abd)
  abd_prop <- sweep(abd, 1, abd_sum, "/")

  sample_meta <- sample_data(ps)
  if (!group %in% names(sample_meta)) {
    stop("`group` must in the field of sample meta data")
  }
  groups <- sample_meta[[group]]


  if (method == "anova") {
    aov_df <- mutate(abd_prop, groups = groups)
    features <- names(abd_prop)
    formula_char <- paste(features, "~", "groups")
    pvalue <- purrr::map(
      formula_char,
      ~ aov(as.formula(.x), abd_prop) %>% summary(.)) %>%
      purrr::map_dbl(~.x[[1]][["Pr(>F)"]][1])
  } else {
    pvalue <- purrr::map_dbl(abd_prop, ~ kruskal.test(.x, groups)$p.value)
  }
  pvalue[is.na(pvalue)] <- 1

  # p value correction for multiple comparisons
  pvalue_corrected <- p.adjust(pvalue, method = p_adjust)

  ef <- purrr::map_dbl(abd_prop, calc_etasq, groups)

  # freq means prop
  abd_prop_groups <- split(abd_prop, groups)
  freq_means_prop <- purrr::map(abd_prop_groups, ~ colMeans(.x)*100)
  mean_names <- paste(
    names(abd_prop_groups),
    "mean_rel_freq_percent",
    sep = "_"
  )

  res <- bind_cols(
    data.frame(
      pvalue = pvalue,
      pvalue_corrected = pvalue_corrected,
      effect_size = ef
    ),
    freq_means_prop
  )
  names(res) <- c("pvalue", "pvalue_corrected", "effect_size", mean_names)

  # filter: pvalue and effect size
  res_filtered <- filter(res, .data$pvalue_corrected <= p_value_cutoff)
  if (!is.null(effect_size_cutoff)) {
    res_filtered <- filter(
      res_filtered,
      .data$effect_size >= effect_size_cutoff
    )
  }

  res_filtered
}


#' calculate eta-squared measurement of effect size commonly used in multiple
#' group statistical analysis
#' @param feature numeric vector, relative of abundance of a given feature
#' @param group vector in the same length with argument `feature`, groups of the
#' feature
#' @noRd
calc_etasq <- function(feature, group) {
  group_n <- table(group)

  if (any(group_n < 1)) {
    return(-1)
  }

  total_sum <- sum(feature)
  n <- length(feature)
  grand_mean <- total_sum/n

  total_ss <- sum((feature - grand_mean)^2)
  feature_split <- split(feature, group)
  between_group_ss <- purrr::map_dbl(feature_split, ~ sum(.x)*sum(.x)/length(.x))
  between_group_ss <- sum(between_group_ss) - total_sum*total_sum/n

  etasq <- ifelse(total_ss == 0, -1, between_group_ss/total_ss)

  etasq
}
