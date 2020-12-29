#' Statistical test for multiple groups
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation (default).
#'   * "log10", the transformation is `log10(object)`, and if the data contains
#'     zeros the transformation is `log10(1 + object)`.
#'   * "log10p", the transformation is `log10(1 + object)`.
#' @param norm the methods used to normalize the microbial abundance data. See
#'   [`normalize()`] for more details.
#'   Options include:
#'   * a integer, e.g. 1e6 (default), indicating pre-sample normalization of
#'     the sum of the values to 1e6.
#'   * "none": do not normalize.
#'   * "rarefy": random subsampling counts to the smallest library size in the
#'     data set.
#'   * "TSS": total sum scaling, also referred to as "relative abundance", the
#'     abundances were normalized by dividing the corresponding sample library
#'     size.
#'   * "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
#'     The scaling factor is then derived using a weighted trimmed mean over the
#'     differences of the log-transformed gene-count fold-change between the
#'     sample and the reference.
#'   * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#'   * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#'   * "CLR": centered log-ratio normalization.
#' @param norm_para arguments passed to specific normalization methods
#' @param method test method, must be one of "anova" or "kruskal"
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param p_value_cutoff numeric, p value cutoff, default 0.05.
#' @param effect_size_cutoff numeric, cutoff of effect size default `NULL`
#' which means no effect size filter
#' @importFrom dplyr mutate bind_cols filter select
#' @importFrom stats p.adjust
#' @seealso [posthoc_test()]
#' @export
test_multiple_groups <- function(ps,
                                 group,
                                 transform = c("identity", "log10", "log10p"),
                                 norm = "TSS",
                                 norm_para = list(),
                                 method = c("anova", "kruskal"),
                                 p_adjust = c("none", "fdr", "bonferroni", "holm",
                                              "hochberg", "hommel", "BH", "BY"),
                                 p_value_cutoff = 0.05,
                                 effect_size_cutoff = NULL) {
  stopifnot(inherits(ps, "phyloseq"))

  if (!check_rank_names(ps)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(available_ranks, collapse = ", ")
    )
  }

  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
  )
  method <- match.arg(method, c("anova", "kruskal"))

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)

  ps_summarized <- summarize_taxa(ps)
  # otus <-  otu_table(ps_summarized)
  otus <- abundances(ps_summarized, norm = FALSE)
  # normalize
  norm_para <- c(norm_para, method = norm, object = list(ps_summarized))
  ps_normed <- do.call(normalize, norm_para)

  feature <- tax_table(ps_summarized)@.Data[, 1]
  abd <- transpose_and_2df(otus)
  abd_norm <- abundances(ps_normed, norm = TRUE) %>%
    transpose_and_2df()

  sample_meta <- sample_data(ps_summarized)
  if (!group %in% names(sample_meta)) {
    stop("`group` must in the field of sample meta data")
  }
  groups <- sample_meta[[group]]

  if (method == "anova") {
    aov_df <- mutate(abd_norm, groups = groups)

    # separator "|" and some strings (such as "/", "-", "+") have a special
    # meaning in formula
    # replace this strings with ___(three underscores) before aov (new_feature),
    # and reset the names as `feature`
    names(aov_df) <- gsub("[-|+*//]", "___", names(aov_df))
    new_features <- setdiff(names(aov_df), "groups")

    formula_char <- paste(new_features, "~", "groups")
    pvalue <- purrr::map(
      formula_char,
      ~ aov(as.formula(.x), aov_df) %>% summary(.)) %>%
      purrr::map_dbl(~.x[[1]][["Pr(>F)"]][1])
  } else {
    pvalue <- purrr::map_dbl(abd_norm, ~ kruskal.test(.x, groups)$p.value)
  }
  pvalue[is.na(pvalue)] <- 1

  # p value correction for multiple comparisons
  pvalue_corrected <- p.adjust(pvalue, method = p_adjust)

  ef <- purrr::map_dbl(abd_norm, calc_etasq, groups)

  # mean abundances
  abd_means <- calc_mean(abd_norm, groups)
  row.names(abd_means) <- feature

  # enriched group
  group_enriched_idx <- apply(abd_means, 1, which.max)
  groups_uniq <- unique(groups)
  group_nms <- groups_uniq[charmatch(groups_uniq, names(abd_means))]
  group_enriched <- group_nms[group_enriched_idx]

  res <- bind_cols(
    data.frame(
      enrich_group = group_enriched,
      pvalue = pvalue,
      pvalue_corrected = pvalue_corrected,
      effect_size = ef
    ),
    abd_means
  )

  # append feature
  res <- mutate(res,feature = feature) %>%
    select(.data$feature, .data$enrich_group, everything())
  # row.names(res) <- feature[match(res$feature, feature)]
  row.names(res) <- paste0("feature", seq_len(nrow(res)))

  # filter: pvalue and effect size
  res_filtered <- filter(res, .data$pvalue_corrected <= p_value_cutoff)

  if (!is.null(effect_size_cutoff)) {
    res_filtered <- filter(
      res_filtered,
      .data$effect_size >= effect_size_cutoff
    )
  }

  # summarized tax table
  tax <- matrix(feature) %>%
    tax_table()
  row.names(tax) <- row.names(otus)

  if (nrow(res_filtered) == 0) {
    warning("No significant features were found, return all the features")
    marker <- microbiomeMarker(
      marker_table = marker_table(res),
      tax_table_orig = tax_table(ps),
      otu_table(t(abd), taxa_are_rows = TRUE),
      tax
    )
  } else {
    row.names(res_filtered) <- paste0("marker", seq_len(nrow(res_filtered)))
    marker <- microbiomeMarker(
      marker_table = marker_table(res_filtered),
      tax_table_orig = tax_table(ps),
      otu_table(t(abd), taxa_are_rows = TRUE),
      tax
    )
  }

  marker
}

# calculate mean abundance of each feature in each group
#' @importFrom dplyr bind_cols
#' @noRd
calc_mean <- function(abd_norm, groups) {
  abd_norm_groups <- split(abd_norm, groups)
  abd_means <- purrr::map(abd_norm_groups, ~ colMeans(.x)) %>%
    bind_cols() %>%
    as.data.frame()
  row.names(abd_means) <- names(abd_norm)
  names(abd_means) <- paste(
    names(abd_norm_groups),
    "mean_abundance",
    sep = ":"
  )

  abd_means
}

#' calculate eta-squared measurement of effect size commonly used in multiple
#' group statistical analysis
#' @param feature numeric vector, abundance of a given feature
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
