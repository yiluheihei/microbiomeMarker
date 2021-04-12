#' Differential analysis of compositions of microbiomes with bias correction
#' (ANCOM-BC).
#'
#' Differential abundance analysis for microbial absolute abundance data. This
#' function is a wrapper of [`ANCOMBC::ancombc()`].
#'
#' @param ps  a [`phyloseq::phyloseq-class`] object, which consists of a feature
#'   table, a sample metadata and a taxonomy table.
#' @param formula the character string expresses how the microbial absolute
#'   abundances for each taxon depend on the variables in metadata.
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
#' @param p_adjust method to adjust p-values by. Default is "holm".
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none". See [`stats::p.adjust()`] for more details.
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of
#'   zeroes greater than `zero_cut` will be excluded in the analysis. Default
#'   is 0.90.
#' @param lib_cut a numerical threshold for filtering samples based on library
#'   sizes. Samples with library sizes less than `lib_cut` will be excluded
#'   in the analysis. Default is 0, i.e. do not filter any sample.
#' @param group_var the name of the group variable in metadata. Specifying
#'   `group_var` is required for detecting structural zeros and performing
#'   global test.
#' @param struc_zero whether to detect structural zeros. Default is FALSE.
#' @param neg_lb whether to classify a taxon as a structural zero in the
#'   corresponding study group using its asymptotic lower bound.
#'   Default is FALSE.
#' @param tol the iteration convergence tolerance for the E-M algorithm.
#'   Default is 1e-05.
#' @param max_iter the maximum number of iterations for the E-M algorithm.
#'   Default is 100.
#' @param conserve whether to use a conservative variance estimate of
#'   the test statistic. It is recommended if the sample size is small and/or
#'   the number of differentially abundant taxa is believed to be large.
#'   Default is FALSE.
#' @param alpha level of significance. Default is 0.05.
#'
#' @references
#' Lin, Huang, and Shyamal Das Peddada. "Analysis of compositions of microbiomes
#' with bias correction." Nature communications 11.1 (2020): 1-11.
#'
#' @seealso [`ANCOMBC::ancombc`]
#'
#' @importFrom ANCOMBC ancombc
#' @export
run_ancombc <- function(ps,
                        formula,
                        transform = c("identity", "log10", "log10p"),
                        norm = "none",
                        norm_para = list(),
                        p_adjust = "holm",
                        zero_cut = 0.9,
                        lib_cut = 0,
                        group_var = NULL,
                        struc_zero = FALSE,
                        neg_lb = FALSE,
                        tol = 1e-05,
                        max_iter = 100,
                        conserve = FALSE,
                        alpha = 0.05) {
  stopifnot(inherits(ps, "phyloseq"))

  # check whether group_var is valid, write a function
  sample_meta <- sample_data(ps)
  meta_nms <- names(sample_meta)
  if (!is.null(group_var)) {
    if (!group_var %in% meta_nms) {
      stop(
        group_var, " are not contained in the `sample_data` of `ps`",
        call. = FALSE
      )
    }
  }

  # if it contains missing values for any
  # variable specified in the formula, the corresponding sampling fraction
  # estimate for this sample will return NA since the sampling fraction is
  # not estimable with the presence of missing values.
  # remove this samples
  vars_formula <- all.vars(stats::as.formula(paste("~", formula)))
  for (var in vars_formula) {
    ps <- remove_na_samples(ps, var)
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
  ps_normed <- summarize_taxa(ps_normed)

  # ancombc differential abundance analysis
  ancombc_out <- ANCOMBC::ancombc(
    ps_normed,
    formula = formula,
    p_adj_method = p_adjust,
    zero_cut = zero_cut,
    lib_cut = lib_cut,
    group = group_var,
    struc_zero = struc_zero,
    neg_lb = neg_lb,
    tol = tol,
    max_iter = max_iter,
    conserve = conserve,
    alpha = alpha,
    global = TRUE
  )

  groups <- sample_data(ps_normed)[[group_var]]
  n_group <- length(unique(groups))

  # multiple-group comparison will be performed while the group
  # variable has > 2 levels
  if (n_group > 2) {
    # ANCOM-BC global test to determine taxa that are differentially abundant
    # between three or more groups of multiple samples.
    # global result to marker_table
    mtab <- ancombc_out$res_global
  } else {
    ancombc_out_res <- ancombc_out$res
    mtab <- do.call(
      cbind,
      ancombc_out_res[c("W", "p_val", "q_val", "diff_abn")]
    )
    names(mtab) <- c("W", "p_val", "q_val", "diff_abn")
  }

  # enriched group
  enrich_abd <- get_ancombc_enrich_group(ps_normed, ancombc_out, group_var)
  norm_abd <- enrich_abd$abd
  group_enrich <- enrich_abd$group_enrich
  idx <- match(group_enrich$feature, rownames(mtab))
  group_enrich <- group_enrich[idx, ]

  mtab <- cbind(mtab, group_enrich)
  mtab <- mtab[mtab$diff_abn, ]
  mtab <- mtab[c("feature", "enrich_group", "W", "p_val", "q_val")]
  names(mtab) <- c("feature", "enrich_group", "effect_size", "pvalue", "padj")

  mm <- microbiomeMarker(
    marker_table = marker_table(mtab),
    norm_method = get_norm_method(norm),
    diff_method = "ancombc",
    sam_data = sample_data(ps_normed),
    # tax_table = tax_table(ps),
    otu_table = otu_table(norm_abd, taxa_are_rows = TRUE),
    tax_table = tax_table(ps_normed)
  )

  mm
}

get_ancombc_enrich_group <- function(ps, ancombc_out, group_var) {
  samp_frac <- ancombc_out$samp_frac
  # As shown in the ancombc vignette: if it contains missing values for any
  # variable specified in the formula, the corresponding sampling fraction
  # estimate for this sample will return NA since the sampling fraction is
  # not estimable with the presence of missing values.
  # Replace NA with 0
  samp_frac[is.na(samp_frac)] <- 0

  # Add pesudo-count (1) to avoid taking the log of 0
  log_abd <- log(abundances(ps, norm = TRUE) + 1)
  # Adjust the log observed abundances
  log_abd_adj <- sweep(log_abd, 2, samp_frac)

  groups <- sample_data(ps)[[group_var]]
  # remove groups with NA
  na_idx <- is.na(groups)
  log_abd_adj <- log_abd_adj[, !na_idx]
  groups <- groups[!na_idx]

  # mean absolute abundance
  abd_mean <- by(t(log_abd_adj), groups, colMeans)
  abd_mean <- do.call(cbind, abd_mean)
  idx_enrich <- apply(abd_mean, 1, which.max)
  group_enrich <- colnames(abd_mean)[idx_enrich]
  group_enrich <- data.frame(
    feature = rownames(abd_mean),
    enrich_group = group_enrich
  )

  list(abd = exp(log_abd_adj), group_enrich = group_enrich)
}
