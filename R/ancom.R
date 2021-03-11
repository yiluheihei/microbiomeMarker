#' Perform differential analysis using ANCOM
#'
#' Perform significant test by comparing the pairwise log ratios (based on
#' relative abundances) between all features.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param group_var character, the variable to set the group.
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
#' @param norm_para  named `list`. other arguments passed to specific
#'   normalization methods.  Most users will not need to pass any additional
#'   arguments here.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param alpha significance level for each of the statistical tests,
#'   default 0.05.
#' @param theta lower bound for the proportion for the W-statistic, default 0.7.
#' @param adj_formula character string representing the formula for adjustment.
#' @param rand_formula character string representing the formula for random effects in lme
#' @param ... additional arguments passed to the test.
#' @return a [microbiomeMarker-class] object, in which the `slot` of `marker_table`
#' contains four variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched.
#' *  `effect_size`, differential means for two groups, or F statistic for more
#'   than two groups
#' * `W`, the W-statistic, number of features that a single feature is tested to
#'   be significantly different against.
#' @author Huang Lin, Yang Cao
#' @export
run_ancom <- function(ps,
                      group_var,
                      transform = c("identity", "log10", "log10p"),
                      norm = "TSS",
                      norm_para = list(),
                      p_adjust = c("none", "fdr", "bonferroni", "holm",
                                   "hochberg", "hommel", "BH", "BY"),
                      alpha = 0.05,
                      theta = 0.75,
                      adj_formula = NULL,
                      rand_formula = NULL,
                      ...) {
  stopifnot(inherits(ps, "phyloseq"))

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

  # ANOCOM requires relative abundance
  ps_normed <- norm_tss(ps_normed)

  # summarize data
  ps_summarized <- summarize_taxa(ps_normed)
  # ps_summarized = ps_normed

  feature_table <- abundances(ps_summarized, norm = TRUE)
  meta_data <- data.frame(sample_data(ps_summarized))
  cls_full <- meta_data[[group_var]]
  cls_n <- length(unique(cls_full))

  # effect size: CLR mean_difference or aov f statistic
  feature_table_clr <- norm_clr(otu_table(feature_table, taxa_are_rows = TRUE))
  ef <- apply(
    feature_table_clr,
    1,
    calc_ancom_ef, group = cls_full
  )
  # names(ef) <- row.names(feature_table)
  # enrich_group
  group_enriched <- apply(
    feature_table_clr,
    1,
    get_ancom_enrich_group, group = cls_full
  )

  # ANCOM requires log transformation
  feature_table <- log(as.matrix(feature_table) + 1)
  n_taxa <- nrow(feature_table)
  taxa_id <- row.names(feature_table)
  n_samp <- ncol(feature_table)

  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    if (cls_n == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun <- stats::wilcox.test
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun <- stats::kruskal.test
    }
    # Formula
    tformula <- formula(paste("x ~", group_var, sep = " "))
  } else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun <- stats::aov
    # Formula
    tformula <- formula(paste("x ~", group_var, "+", adj_formula, sep = " "))
  } else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun <- nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula <- formula(paste("x ~", group_var))
    }else {
      # Random coefficients/slope model
      tformula <- formula(paste("x ~", group_var, "+", adj_formula))
    }
  }

  # Calculate the p-value for each pairwise comparison of taxa.
  p_data <- matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) <-  taxa_id
  rownames(p_data) <- taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using
    # taxon i as the reference. e.g. the first alr matrix will be the log
    # abundance data (feature_table) recursively subtracted by the log abundance
    # of 1st taxon (1st column) column-wisely,
    # and remove the first i columns since: the first (i - 1) columns were
    # calculated by previous iterations, and the i^th column contains all zeros.
    # if feature_table is m x n, alr_data is n x m
    alr_data <- apply(feature_table, 1, function(x) x - feature_table[i, ])

    # apply(...) allows crossing the data in a number of ways and avoid explicit
    # use of loop constructs. Here, we basically want to iteratively subtract
    # each column of the feature_table by its i^th column.
    alr_data <- alr_data[, - (1:i), drop = FALSE]
    n_lr <- ncol(alr_data) # number of log-ratios (lr)
    alr_data <- cbind(alr_data, meta_data) # merge with the metadata

    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] <- apply(
        alr_data[, 1:n_lr, drop = FALSE],
        2,
        function(x){
          suppressWarnings(
            tfun(
              tformula,
              data = data.frame(x, alr_data, check.names = FALSE))$p.value
          )
        }
      )
    } else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] <- apply(
        alr_data[, 1:n_lr, drop = FALSE],
        2,
        function(x){
          fit = tfun(
            tformula,
            data = data.frame(x, alr_data, check.names = FALSE),
            na.action = na.omit
          )
          summary(fit)[[1]][group_var, "Pr(>F)"]
        }
      )
    } else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] <- apply(
        alr_data[, 1:n_lr, drop = FALSE],
        2,
        function(x){
          fit = tfun(
            fixed = tformula,
            data = data.frame(x, alr_data, check.names = FALSE),
            random = formula(rand_formula),
            na.action = na.omit,
            ...
          )
          stats::anova(fit)[group_var, "p-value"]
        }
      )
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] <- t(p_data)[upper.tri(p_data)]
  diag(p_data) <- 1 # let p-values on diagonal equal to 1
  p_data[is.na(p_data)] <- 1 # let p-values of NA equal to 1

  # Multiple comparisons correction.
  q_data <- apply(p_data, 2, function(x) p.adjust(x, method = p_adjust))

  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W <- apply(q_data, 2, function(x) sum(x < alpha))

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
  marker <- out_comp[out_comp$W > theta * (n_taxa -1), ]
  if (cls_n == 2) {
    names(marker)[3] <- "CLR_diff_mean"
  } else {
    names(marker)[3] <- "CLR_F_statistic"
  }
  marker <- marker_table(marker)
  tax <- matrix(taxa_id) %>%
    tax_table()
  row.names(tax) <- row.names(feature_table)

  mm <- microbiomeMarker(
    marker_table = marker,
    tax_table_orig = tax_table(ps),
    otu_table(feature_table, taxa_are_rows = TRUE), # normalized feature table
    sample_data(ps),
    tax
  )

  mm
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


