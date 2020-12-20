# post hoc test -----------------------------------------------------------

#' Post hoc pairwise comparisons for multiple groups test.
#'
#' Multiple group test, such as anova and Kruskal-Wallis rank sum test, can be
#' used to uncover the significant feature among all groups. Post hoc tests are
#' used to uncover specific mean differences between pair of groups.
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
                         transform = c("identity", "log10", "log10p"),
                         norm = "TSS",
                         norm_para = list(),
                         conf_level = 0.95,
                         method = c("tukey", "games_howell",
                                    "scheffe", "welch_uncorrected")) {
  stopifnot(inherits(ps, "phyloseq"))

  if (!check_rank_names(ps)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(available_ranks, collapse = ", ")
    )
  }

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  method <- match.arg(
    method,
    c("tukey", "games_howell", "scheffe", "welch_uncorrected")
  )

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)

  # ps_summarized <- summarize_taxa(ps)
  # otus <- otu_table(ps_summarized)
  # # normalize
  # norm_para <- c(norm_para, method = norm, object = list(otus))
  # otus_norm <- do.call(normalize, norm_para)

  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_norm <- do.call(normalize, norm_para)
  ps_summarized <- summarize_taxa(ps_norm)

  # relative abundance
  # abd <- transpose_and_2df(otus)
  # abd_norm <- transpose_and_2df(otus_norm)
  abd_norm <- transpose_and_2df(otu_table(ps_summarized))
  feature <- tax_table(ps_summarized)@.Data[, 1]
  names(abd_norm) <- feature
  # names(abd) <- feature
  # names(abd_norm) <- tax_table(ps)@.Data[, rank_name]

  groups <- sample_data(ps_summarized)[[group]]

  result = switch(
    method,
    tukey = purrr::map(
      abd_norm,
      calc_tukey_test,
      groups,
      conf_level
    ),
    games_howell = purrr::map(
      abd_norm,
      calc_games_howell_test,
      groups,
      conf_level
    ),
    scheffe = purrr::map(
      abd_norm,
      calc_scheffe_test,
      groups,
      conf_level
    ),
    welch_uncorrected = purrr::map(
      abd_norm,
      calc_welch_uncorrected_test,
      groups,
      conf_level
    )
  )

  # diff_means to diff_mean
  result <- purrr::map(result, ~ mutate(
    .x,
    comparions = .data$comparisons,
    diff_mean = .data$diff_means,
    pvalue = .data$pvalue,
    ci_lower = .data$ci_lower,
    ci_upper = .data$ci_upper,
    .keep = "none"
  ))

  abundance <- abd_norm
  abundance$group <- groups
  postHocTest(
    result = DataFrameList(result),
    abundance = abundance,
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
