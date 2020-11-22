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
#' @param p_value_cutoff numeric, p value cutoff, default 0.05
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

  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
  )
  method <- match.arg(method, c("anova", "kruskal"))

  ranks <- rank_names(ps)
  diff_rank <- setdiff(ranks, available_ranks)
  if (length(diff_rank)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(available_ranks, collapse = ", ")
    )
  }

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)

  otus <- summarize_taxa(ps)
  # normalize
  norm_para <- c(norm_para, method = norm, object = list(otus))
  otus_norm <- do.call(normalize, norm_para)

  feature <- taxa_names(otus)
  abd <- transpose_and_2df(otus)
  abd_norm <- transpose_and_2df(otus_norm)

  sample_meta <- sample_data(ps)
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

  # freq means prop
  freq_means_prop <- calc_mean_prop(abd_norm, groups)
  row.names(freq_means_prop) <- feature

  # enriched group
  group_enriched_idx <- apply(freq_means_prop, 1, which.max)
  groups_uniq <- unique(groups)
  group_nms <- groups_uniq[charmatch(groups_uniq, names(freq_means_prop))]
  group_enriched <- group_nms[group_enriched_idx]

  res <- bind_cols(
    data.frame(
      enrich_group = group_enriched,
      pvalue = pvalue,
      pvalue_corrected = pvalue_corrected,
      effect_size = ef
    ),
    freq_means_prop
  )

  # append feature
  res <- mutate(res,feature = feature) %>%
    select(.data$feature, .data$enrich_group, everything())
  row.names(res) <- feature[match(res$feature, feature)]

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
  row.names(tax) <- feature

  if (nrow(res_filtered) == 0) {
    warning("No significant features were found, return all the features")
    marker <- microbiomeMarker(
      marker_table = marker_table(res),
      tax_table_orig = tax_table(ps),
      otu_table(t(abd), taxa_are_rows = TRUE),
      tax
    )
  } else {
    marker <- microbiomeMarker(
      marker_table = marker_table(res_filtered),
      tax_table_orig = tax_table(ps),
      otu_table(t(abd), taxa_are_rows = TRUE),
      tax
    )
  }

  marker
}

# calculate mean proportion of each feature in each group
#' @importFrom dplyr bind_cols
#' @noRd
calc_mean_prop <- function(abd_norm, groups) {
  abd_norm_groups <- split(abd_norm, groups)
  freq_means_prop <- purrr::map(abd_norm_groups, ~ colMeans(.x)*100) %>%
    bind_cols() %>%
    as.data.frame()
  row.names(freq_means_prop) <- names(abd_norm)
  names(freq_means_prop) <- paste(
    names(abd_norm_groups),
    "mean_rel_freq_percent",
    sep = ":"
  )

  freq_means_prop
}

# effect  size ------------------------------------------------------------

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

# post hoc test -----------------------------------------------------------

#' Post hoc pairwise comparisons for multiple groups test.
#'
#' Multiple group test, such as anova and Kruskal-Wallis rank sum test, can be
#' used to uncover the significant feature among all groups. Post hoc tests are
#' used to uncover specific mean differences between pair of groups.
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
#' @param rank_name character, taxonomic names of [`phyloseq::phyloseq-class`]
#' to compare
#' @param conf_level confidence level, default 0.95
#' @param method one of "tukey", "games_howell", "scheffe", "welch_uncorrected",
#' defining the method for the pairwise comparisons. See details for more
#' information.
#' @return a [postHocTest-class] object
#' @seealso [postHocTest-class], [test_multiple_groups()]
#' @importFrom IRanges DataFrameList
#' @importFrom dplyr mutate
#' @export
posthoc_test <- function(ps,
                         group,
                         rank_name,
                         conf_level = 0.95,
                         method = c("tukey", "games_howell", "scheffe", "welch_uncorrected")) {
  stopifnot(inherits(ps, "phyloseq"))

  method <- match.arg(
    method,
    c("tukey", "games_howell", "scheffe", "welch_uncorrected")
  )

  ranks <- rank_names(ps)
  if (!rank_name %in% ranks) {
    stop("`rank_name` must be one of available taxonomic ranks of `ps`")
  }
  # agglomerate tax in the same rank_name
  if (rank_name != ranks[length(ranks)]) {
    ps <- tax_glom(ps, taxrank = rank_name)
  }

  # relative abundance
  abd <- otu_table(ps) %>%
    t() %>%
    as.data.frame()
  abd_sum <- rowSums(abd)
  abd_norm <- sweep(abd, 1, abd_sum, "/")
  names(abd_norm) <- tax_table(ps)@.Data[, rank_name]

  groups <- sample_data(ps)[[group]]

  # mean proportion of each feature in each group
  # freq_means_prop <- calc_mean_prop(abd_norm, groups)
  # row.names(freq_means_prop) <- tax_table(ps)@.Data[, rank_name]

  result = switch(method,
    tukey = purrr::map(abd_norm,
      calc_tukey_test,
      groups,
      conf_level
    ),
    games_howell = purrr::map(abd_norm,
      calc_games_howell_test,
      groups,
      conf_level
    ),
    scheffe = purrr::map(abd_norm,
      calc_scheffe_test,
      groups,
      conf_level
    ),
    welch_uncorrected = purrr::map(abd_norm,
      calc_welch_uncorrected_test,
      groups, conf_level
    )
  )

  # convert to mean proportion
  result <- purrr::map(result, ~ mutate(
    .x,
    comparions = .data$comparisons,
    diff_mean_prop = .data$diff_means * 100,
    pvalue = .data$pvalue,
    ci_lower_prop = .data$ci_lower * 100,
    ci_upper_prop = .data$ci_upper * 100,
    .keep = "none"
  ))

  abundance_proportion <- abd_norm*100
  abundance_proportion$group <- groups
  postHocTest(
    result = DataFrameList(result),
    abundance_proportion = abundance_proportion,
    conf_level = conf_level,
    method = method,
    method_str = paste("Posthoc multiple comparisons of means:", method)
  )
}

#' Tukey Post-hoc Tests
#' @param obs numeric vector, relative abundance of a feature
#' @param groups character vector, the same length of the argument `obs`
#' @param conf_level confidence level, default 0.95
#' @importFrom stats aov TukeyHSD
#' @importFrom dplyr mutate
#' @noRd
calc_tukey_test <- function(obs, groups, conf_level = 0.95) {
  df <- data.frame(obs = obs, groups = groups)
  fit <- aov(obs ~ groups, df)
  tukey_groups <- TukeyHSD(fit, "groups", conf.level = conf_level)$groups
  res <- data.frame(comparisons = row.names(tukey_groups), tukey_groups) %>%
    mutate(
      comparisons = .data$comparisons,
      diff_means = .data$diff,
      pvalue = .data$p.adj,
      ci_lower = .data$lwr,
      ci_upper = .data$upr,
      .keep = "none"
    )

  res
}

#' Games Howell Post-hoc Tests
#' @importFrom stats var pairwise.table ptukey qtukey
#' @inheritParams calc_tukey_test
#' @references https://github.com/kassambara/rstatix/blob/master/R/games_howell_test.R
#' @noRd
calc_games_howell_test <- function(obs, groups, conf_level = 0.95) {
  groups <- factor(groups)
  groups_n <- length(levels(groups))
  if(groups_n == 1) {
    stop("The number of groups at least 3")
  }

  # Statistics for games howell tests
  grp_sizes <- tapply(obs, groups, length)
  # nb.groups <- length(grp.sizes)
  grp_means <- tapply(obs, groups, mean)
  grp_vars <- tapply(obs, groups, var)

  # Helper functions
  get_mean_diff <- function(i, j){
    grp_means[i] - grp_means[j]
  }
  get_weltch_sd <- function(i, j){
    vn1 <- grp_vars[i]/grp_sizes[i]
    vn2 <- grp_vars[j]/grp_sizes[j]
    if (vn1 == 0) {
      vn1 = 1e-6
    }
    if (vn2 ==0) {
      vn2 = 1e-6
    }
    sqrt(vn1 + vn2)
  }
  get_degree_of_freedom <- function(i, j){
    vn1 <- grp_vars[i]/grp_sizes[i]
    vn2 <- grp_vars[j]/grp_sizes[j]
    if (vn1 == 0) {
      vn1 = 1e-6
    }
    if (vn2 ==0) {
      vn2 = 1e-6
    }
    A <- (vn1 + vn2)^2
    B <- (vn1^2)/(grp_sizes[i] - 1)
    C <- (vn2^2)/(grp_sizes[j] - 1)
    A/(B+C)
  }
  correct_pairwise_table <- function(table) {
    comparisons <- purrr::map(
      colnames(table),
      ~ paste(row.names(table), "-", .x, sep = "")) %>%
      unlist()
    res <- c(table[, 1], table[, 2])
    names(res) <- comparisons
    res <- res[!is.na(res)]

    res
  }

  mean_diff_table <- pairwise.table(
    get_mean_diff, levels(groups),
    p.adjust.method = "none")
  mean_diff <- correct_pairwise_table(mean_diff_table)

  weltch_sd <- pairwise.table(
    get_weltch_sd, levels(groups),
    p.adjust.method = "none") %>%
    correct_pairwise_table()

  df <- pairwise.table(
    get_degree_of_freedom, levels(groups),
    p.adjust.method = "none") %>%
    correct_pairwise_table()

  t <- abs(mean_diff)/weltch_sd
  p <- ptukey(t*sqrt(2), groups_n, df, lower.tail = FALSE)
  se <- weltch_sd*sqrt(0.5)

  q <- qtukey(p = conf_level, groups_n, df = df)
  conf_high <- mean_diff + q*se
  conf_low <- mean_diff - q*se

  res <- data.frame(
    comparisons = names(mean_diff),
    diff_means = mean_diff,
    pvalue = p,
    ci_lower = conf_low,
    ci_upper = conf_high
  )

  res
}

#' scheffe Post-hoc Tests
#' @inheritParams calc_tukey_test
#' @importFrom stats aov model.tables pf qf
#' @noRd
#' @references https://github.com/AndriSignorell/DescTools/blob/master/R/Tests.r#L3755:1
calc_scheffe_test<- function(obs, groups, conf_level=0.95){
  x <- aov(obs ~ groups)
  mm <- model.tables(x, "means")
  if (is.null(mm$n))
    stop("no factors in the fitted model")
  tab <- mm$tables[["groups"]]
  MSE <- sum(x$residuals^2) / x$df.residual

  group_means <- as.vector(tab)
  nms <- names(tab)
  group_size <- mm$n[["groups"]]

  if (length(group_size) < length(group_means)) {
    group_size <- rep.int(group_size, length(group_means))
  }

  contrasts <- Contrasts(nms)
  diff_means <- apply(contrasts * group_means, 2, sum)
  sscoeff <- apply(contrasts * contrasts / group_size, 2, sum)

  dferr <- x$df.residual
  dfgrp <- length(x$residuals) - dferr - 1

  pvalue <- pf(
    diff_means^2/(MSE * sscoeff * dfgrp),
    df1 = dfgrp,
    df2 = dferr,
    lower.tail = FALSE
  )

  critvalue <- dfgrp * qf(1 - conf_level, dfgrp, dferr, lower.tail=FALSE)
  ci_lower <- diff_means - sqrt(critvalue) * sqrt(MSE * sscoeff)
  ci_upper <- diff_means + sqrt(critvalue) * sqrt(MSE * sscoeff)

  res <- data.frame(
    comparisons = names(diff_means),
    diff_means = diff_means,
    pvalue = pvalue,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

  res
}

Contrasts <- function (levs) {
  k = length(levs)
  M = data.frame(levs = levs)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      con = rep(0, k)
      con[i] = -1
      con[j] = 1
      nm = paste(levs[j], levs[i], sep = "-")
      M[[nm]] = con
    }
  }
  row.names(M) = levs

  return(M[-1])
}

#' welch's uncorrected Post-hoc Tests
#' @inheritParams calc_tukey_test
#' @noRd
calc_welch_uncorrected_test <- function(obs, groups, conf_level = 0.95) {
  group_means <- tapply(obs, groups, mean)
  nms <- names(group_means)

  contrasts <- Contrasts(nms)
  diff_means <- apply(contrasts * group_means, 2, sum)

  obs_split <- split(obs, groups)
  comparisons <- names(diff_means)
  comparison_groups <- strsplit(comparisons, "-", fixed = TRUE)

  welch_res <- purrr::map(
    comparison_groups,
    ~ t.test(obs_split[[.x[1]]], obs_split[[.x[2]]], conf.level = conf_level)
  )

  pvalue <- purrr::map_dbl(welch_res, ~ .x$p.value)
  ci_lower <- purrr::map_dbl(welch_res, ~ .x$conf.int[1])
  ci_upper <- purrr::map_dbl(welch_res, ~ .x$conf.int[2])

  res <- data.frame(
    comparisons = comparisons,
    diff_means = diff_means,
    pvalue = pvalue,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )

  res
}
