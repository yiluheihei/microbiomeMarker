#' Perform differential analysis using ANCOM
#'
#' Perform significant test by comparing the pairwise log ratios between all
#' features.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param group_var character, the variable to set the group.
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of `phyloseq::rank_names(phyloseq)`,
#'   or "all" means to summarize the taxa by the top taxa ranks
#'   (`summarize_taxa(ps, level = rank_names(ps)[1])`), or "none" means perform
#'   differential analysis on the original taxa (`taxa_names(phyloseq)`, e.g.,
#'   OTU or ASV).
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
#'   * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' @param norm_para  named `list`. other arguments passed to specific
#'   normalization methods.  Most users will not need to pass any additional
#'   arguments here.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param pvalue_cutoff significance level for each of the statistical tests,
#'   default 0.05.
#' @param W_cutoff lower bound for the proportion for the W-statistic, default 0.7.
#' @param test character, the test to dtermine the p value of log ratio,
#'   one of "aov", "wilcox.test",  "kruskal.test".
#' @param ... additional arguments passed to the test function.
#'
#' @details
#' In an experiment with only two treatments, this tests the following
#' hypothesis for feature \eqn{i}:
#'
#' \deqn{H_{0i}: E(log(\mu_i^1)) =  E(log(\mu_i^2))}
#'
#' where \eqn{\mu_i^1} and \eqn{\mu_i^2} are the mean abundances for feature
#' \eqn{i} in the two groups.
#'
#' The developers of this method recommend the following significance tests
#' if there are 2 groups, use non-parametric Wilcoxon rank sum test
#' [`stats::wilcox.test()`]. If there are more than 2 groups, use nonparametric
#' [`stats::kruskal.test()`] or one-way ANOVA [`stats::aov()`].
#'
#' @return a [microbiomeMarker-class] object, in which the `slot` of `marker_table`
#' contains four variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched.
#' *  `effect_size`, differential means for two groups, or F statistic for more
#'   than two groups
#' * `W`, the W-statistic, number of features that a single feature is tested to
#'   be significantly different against.
#'
#' @references Mandal et al. "Analysis of composition of microbiomes: a novel
#' method for studying microbial composition", Microbial Ecology in Health
#' & Disease, (2015), 26.
#'
#' @author Huang Lin, Yang Cao
#'
#' @export
run_ancom <- function(ps,
                      group_var,
                      taxa_rank = "all",
                      transform = c("identity", "log10", "log10p"),
                      norm = "TSS",
                      norm_para = list(),
                      p_adjust = c("none", "fdr", "bonferroni", "holm",
                                    "hochberg", "hommel", "BH", "BY"),
                      pvalue_cutoff = 0.05,
                      W_cutoff = 0.75,
                      test = c("aov", "wilcox.test", "kruskal.test"),
                      ...) {
  stopifnot(inherits(ps, "phyloseq"))
  test <- match.arg(test, c("aov", "wilcox.test", "kruskal.test"))

  # check whether group_var is valid, write a function
  sample_meta <- sample_data(ps)
  meta_nms <- names(sample_meta)
  if (!group_var %in% meta_nms) {
    stop(
      group_var, " are not contained in the `sample_data` of `ps`",
      call. = FALSE
    )
  }

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY")
  )

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)

  # normalize the data
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)

  # summarize data
  # create a function, extract_summarize?
  # check taxa_rank
  check_taxa_rank(ps, taxa_rank)
  if (taxa_rank == "all") {
    ps_summarized <- summarize_taxa(ps_normed)
  } else if (taxa_rank =="none") {
    ps_summarized <- extract_rank(ps_normed, taxa_rank)
  } else {
    ps_summarized <-aggregate_taxa(ps_normed, taxa_rank) %>%
      extract_rank(taxa_rank)
  }

  feature_table <- abundances(ps_summarized, norm = TRUE)
  meta_data <- data.frame(sample_data(ps_summarized))
  cls_full <- meta_data[[group_var]]
  cls_n <- length(unique(cls_full))

  # effect size: CLR mean_difference or aov f statistic
  feature_table_clr <- norm_clr(otu_table(feature_table, taxa_are_rows = TRUE))
  # ef <- apply(
  #   feature_table_clr,
  #   1,
  #   calc_ancom_ef, group = cls_full
  # )
  feature_table_clr <- data.frame(t(feature_table_clr))
  ef <- vapply(
    feature_table_clr,
    calc_ancom_ef,
    FUN.VALUE = 0.0,
    group = cls_full
  )
  # names(ef) <- row.names(feature_table)
  # enrich_group
  group_enriched <- vapply(
    feature_table_clr,
    get_ancom_enrich_group,
    FUN.VALUE = character(1),
    group = cls_full
  )

  # ANCOM requires log transformation
  feature_table <- log(as.matrix(feature_table) + 1)
  n_taxa <- nrow(feature_table)
  taxa_id <- row.names(feature_table)
  n_samp <- ncol(feature_table)

  # Calculate the p-value for each pairwise comparison of taxa.
  p <- calc_ancom_pmat(feature_table, cls_full, test, ...)

  # Multiple comparisons correction.
  p_adjusted <- vapply(
    data.frame(p),
    p.adjust,
    FUN.VALUE = numeric(n_taxa),
    method = p_adjust
  )

  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < pvalue_cutoff.
  W <- apply(p_adjusted, 2, function(x) sum(x < pvalue_cutoff))

  # Organize outputs
  out_comp <- data.frame(
    feature = taxa_id,
    enrich_group = group_enriched,
    ef = ef,
    W = W,
    row.names = NULL,
    check.names = FALSE
  )
  # Declare a taxon to be differentially abundant based on the quantile of W
  # statistic. We perform (n_taxa - 1) hypothesis testings on each taxon, so
  # the maximum number of rejections is (n_taxa - 1).
  sig_out <- out_comp[out_comp$W > W_cutoff * (n_taxa -1), ]
  if (cls_n == 2) {
    names(sig_out)[3] <- "ef_CLR_diff_mean"
  } else {
    names(sig_out)[3] <- "ef_CLR_F_statistic"
  }

  marker <- return_marker(sig_out, out_comp)
  tax <- matrix(taxa_id) %>%
    tax_table()
  row.names(tax) <- row.names(feature_table)

  mm <- microbiomeMarker(
    marker_table = marker,
    norm_method = get_norm_method(norm),
    diff_method = "ANCOM",
    # tax_table = tax_table(ps),
    otu_table = otu_table(feature_table, taxa_are_rows = TRUE), # normalized feature table
    sam_data = sample_data(ps_normed),
    tax_table = tax
  )

  mm
}

# https://github.com/biocore/scikit-bio/blob/master/skbio/stats/composition.py#L811
#' Calculates pairwise pvalues between all features
#' @param feature_table matrix-like, logged feature table.
#' @param classes character vector, the same length with `log_ratio`.
#' @param test  character, the test to dtermine the p value of log ratio,
#'   one of "aov", "wilcox.test",  "kruskal.test".
#' @param ... extra arguments passed to the test.
#' @noRd
calc_ancom_pmat <- function(feature_table,
                            classes,
                            test = c("aov", "wilcox.test", "kruskal.test"),
                            ...) {
  test_method <- match.arg(test, c("aov", "wilcox.test", "kruskal.test"))
  # test <- match.fun(test_method)

  taxas <- row.names(feature_table)
  feature_table <- data.frame(t(feature_table))
  taxa_n <- ncol(feature_table)
  p <- matrix(NA, nrow = taxa_n, ncol = taxa_n)
  row.names(p) <- taxas
  colnames(p) <- taxas


  for (i in 1:(taxa_n -1)) {
    new_table <- -(feature_table[(i + 1):taxa_n] - feature_table[[i]])
    p[-(1:i), i] <- vapply(
      new_table,
      calc_ancom_p,
      FUN.VALUE = numeric(1),
      classes = classes, test = test_method, ...
    )
  }

  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p[upper.tri(p)] <- t(p)[upper.tri(p)]
  diag(p) <- 1 # let p-values on diagonal equal to 1
  p[is.na(p)] <- 1 # let p-values of NA equal to 1

  p
}

#' calculate the p value of a pair-wise log ratio
#' @param log_ratio  a numeric vector, a pair-wise log ratio.
#' @param classes character vector, the same length with `log_ratio`.
#' @param test  character, the test to dtermine the p value of log ratio,
#'   one of "aov", "wilcox.test",  "kruskal.test".
#' @param ... extra arguments passed to the test.
#' @noRd
calc_ancom_p <- function(log_ratio, classes, test, ...) {
  test_method <- match.arg(test, c("aov", "wilcox.test", "kruskal.test"))
  test_fun <- match.fun(test_method)

  fml <- formula(paste0("x ~ ", "class"))
  dat <- data.frame(x = log_ratio, class = classes)
  test_res <- test_fun(fml, data = dat, ...)
  pvalue <- ifelse(
    test_method == "aov",
    summary(test_res)[[1]]["class", "Pr(>F)"],
    test_res$p.value
  )

  pvalue
}


#' Identify structural zeros
#'
#' @author Huang Lin, Yang Cao
#' @references \url{https://github.com/FrederickHuangLin/ANCOMBC/blob/master/R/get_struc_zero.R}
#' @noRd
get_struc_zero <- function(ps, group_var, neg_lb) {
  stopifnot(inherits(ps, "phyloseq"))
  stopifnot(is.logical(neg_lb))
  stopifnot(length(group_var) == 1 & is.character(group_var))

  meta_tab <- sample_data(ps)
  check_var_in_meta(group_var, meta_tab)
  group <- factor(meta_tab[[group_var]])

  feature_tab <- as(otu_table(ps), "matrix")
  present_tab <- feature_tab
  present_tab[is.na(present_tab)] <- 0
  present_tab[present_tab != 0] <- 1
  n_taxa <- nrow(feature_tab)
  n_group <- nlevels(group)

  p_hat <- matrix(NA, nrow = n_taxa, ncol = n_group)
  rownames(p_hat) <- rownames(feature_tab)
  colnames(p_hat) <- levels(group)
  samp_size <- p_hat

  for (i in seq_len(n_taxa)) {
    p_hat[i, ] <- tapply(
      present_tab[i, ],
      group,
      function(x) mean(x, na.rm = TRUE)
    )
    samp_size[i, ] <- tapply(
      feature_tab[i, ],
      group,
      function(x) length(x[!is.na(x)])
    )
  }

  p_hat_lo <- p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

  zero_ind <- p_hat == 0

  if (neg_lb) {
    zero_ind[p_hat_lo <= 0] <- TRUE
  }
  colnames(zero_ind) <- paste0(
    "structural_zero (", group_var, " = ",colnames(zero_ind), ")"
  )

  data.frame(zero_ind)
}


#' calculate effect size, mean differences for two groups, and f statistic for
#' three groups
#' https://github1s.com/qiime2/q2-composition/blob/HEAD/q2_composition/_ancom.py
#' @noRd
#' @importFrom stats lm aov
calc_ancom_ef <- function(feature_abd, group) {
  group_n <- length(unique(group))
  if (group_n == 2) {
    ef <- abs(lm(feature_abd ~ group)$coefficients[2])
  } else if (group_n > 2) {
    # f statistic from aov
    ef <- summary(aov(feature_abd ~ group))[[1]]$`F value`[1]
  } else {
    stop("The number of group must be greater than 2", call. = FALSE)
  }

  ef
}

#' enrich group for ancom, rewrite this function in the later
#' split get_feature_enrich_group into two funcitons: enrich_group and log max mean
#' @noRd
get_ancom_enrich_group <- function(feature_abd, group) {
  abd_split <- split(feature_abd, group)
  abd_mean_group <- sapply(abd_split, mean)
  enrich_group <- names(abd_split)[which.max(abd_mean_group)]

  enrich_group
}

#' preprocess feature data using methods of ANCOM-II
#' @noRd
#' @importFrom stats dnorm lm na.omit quantile residuals sd
preprocess_ancom <- function(feature_table,
                             meta_data,
                             sample_var,
                             group_var = NULL,
                             out_cut = 0.05,
                             zero_cut = 0.90,
                             lib_cut,
                             neg_lb) {
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1)
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = rep(0, length(f)); e[!is.na(group)] = residuals(f_fit)
    y = t(t(z) - e)

    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n

        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break

        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }

      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
    out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
      unlist(tapply(i, group, function(j) outlier_check(j)))))

    feature_table[out_ind] = NA
  }

  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }

  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }

  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1

    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)

    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }

  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
  return(res)
}


