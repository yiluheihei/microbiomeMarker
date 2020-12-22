# In the vignette of DESeq2:
# The values in the matrix should be un-normalized counts or estimated counts
# of sequencing reads (for single-end RNA-seq) or fragments (for paired-end
# RNA-seq). The RNA-seq workflow describes multiple techniques for preparing
# such count matrices. It is important to provide count matrices as input for
# DESeq2’s statistical model (Love, Huber, and Anders 2014) to hold, as only
# the count values allow assessing the measurement precision correctly. The
# DESeq2 model internally corrects for library size, so transformed or
# normalized values such as counts scaled by library size should not be used
# as input.
#
# reference source code:
# https://github.com/biocore/qiime/blob/76d633c0389671e93febbe1338b5ded658eba31f/qiime/support_files/R/DESeq2_nbinom.r

#' Perform DESeq analysis
#'
#' Differential expression analysis based on the Negative Binomial distribution
#' using **DESeq2**.
#'
#' @param ps  ps a [`phyloseq::phyloseq-class`] object.
#' @param group_var  character, the variable to set the group, must be one of
#'   the var of the sample metadata
#' @param subgroup1,subgroup2 character, subgroups to be compared.
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation (default).
#'   * "log10", the transformation is `log10(object)`, and if the data contains
#'     zeros the transformation is `log10(1 + object)`.
#'   * "log10p", the transformation is `log10(1 + object)`.
#' @param test,fitType,sfType these three parameters are inherited form  `
#'   [`DESeq2::DESeq()`].
#'   - `test`, should be either "Wald" or "LRT", which will
#'     then use either Wald significance tests, or the likelihood ratio test on
#'     the difference in deviance between a full and reduced model formula.
#'   - `fitType`, either "parametric", "local", "mean", or "glmGamPoi" for the
#'     type of fitting of dispersions to the mean intensity.
#'   - `sfType`, either "ratio", "poscounts", or "iterate" for the type of size
#'     factor estimation.
#' @param p_adjust method for multiple test correction, default `none`, for
#'   more details see [stats::p.adjust].
#' @param p_value_cutoff p_value_cutoff numeric, p value cutoff, default 0.05.
#' @param ... extra parameters passed to [`DESeq2::results()`].
#' @export
#' @return a [`microbiomeMarker-class`] object.
#' @seealso [`DESeq2::results()`],[`DESeq2::DESeq()`]
#' @importFrom stats formula
run_deseq2 <- function(ps,
                      group_var,
                      subgroup1,
                      subgroup2,
                      transform = c("identity", "log10", "log10p"),
                      test = c("Wald", "LRT"),
                      fitType = c("parametric", "local", "mean", "glmGamPoi"),
                      sfType = "poscounts",
                      p_adjust = c("none", "fdr", "bonferroni", "holm",
                                   "hochberg", "hommel", "BH", "BY"),
                      p_value_cutoff = 0.05,
                      ...) {
  # groups
  groups <- sample_data(ps)[[group_var]]
  lvl <- unique(groups)
  n_lvl <- length(lvl)
  if (n_lvl < 2) {
    stop("Differential analysis requires at least two groups.")
  }

  if (!any(c(subgroup1, subgroup2) %in% lvl)) {
    stop(
      paste0(subgroup1, " and ", subgroup2,
             " must be included in the field of `group_var`.")
    )
  }

  test <- match.arg(test, c("Wald", "LRT"))
  fitType <- match.arg(fitType, c("parametric", "local", "mean", "glmGamPoi"))
  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY")
  )

  if (! sfType %in% c("ratio","poscounts","iterate")) {
    stop("`sfType` muste be one of poscounts, ratio, or iterate")
  }

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)
  # summarize data
  ps_summarized <- summarize_taxa(ps)

  design <- formula(paste("~ ", group_var))
  dds <- suppressWarnings(
    phyloseq::phyloseq_to_deseq2(ps_summarized, design = design)
  )

  dds <- DESeq2::DESeq(
    dds,
    test = test,
    fitType = fitType,
    sfType = sfType,
    ...
  )

  # By default, independent filtering is performed to select a set of genes
  # for multiple test correction which maximizes the number of adjusted p-values
  # less than a given critical value alpha (by default 0.1).
  # The adjusted p-values for the genes which do not pass the filter threshold
  # are set to NA.
  # By default, results assigns a p-value of NA to genes containing count
  # outliers, as identified using Cook's distance.
  res <- DESeq2::results(
    object = dds,
    contrast = c(group_var, subgroup1, subgroup2),
    alpha = p_value_cutoff,
    pAdjustMethod = p_adjust
  )

  res_ordered <- res[order(res$padj), ] %>%
    as.data.frame()
  padj <- res_ordered$padj
  res_filtered <- res_ordered[!is.na(padj) & padj < p_value_cutoff, ]

  if (nrow(res_filtered) == 0) {
    warning("No significant features were found, return all the features")
    sig_feature <- cbind(feature = row.names(res_ordered), res_ordered)
  } else {
    sig_feature <- cbind(feature = row.names(res_filtered), res_filtered)
  }

  row.names(sig_feature) <- paste0("marker", seq_len(nrow(sig_feature)))
  sig_feature$enrich_group <- ifelse(
    sig_feature$log2FoldChange > 0,
    subgroup2,
    subgroup1
  )

  # reorder columns: feature, enrich_group, other columns
  other_col <- setdiff(names(sig_feature), c("feature", "enrich_group"))
  sig_feature <- sig_feature[, c("feature", "enrich_group", other_col)]

  # normalized counts
  # https://bioinformatics.stackexchange.com/questions/193/how-can-i-extract-normalized-read-count-values-from-deseq2-results
  counts_normalized <- DESeq2::counts(dds, normalized = TRUE)

  marker <- microbiomeMarker(
    marker_table = marker_table(sig_feature),
    tax_table_orig = tax_table(ps),
    otu_table(counts_normalized, taxa_are_rows = TRUE),
    tax_table(ps_summarized)
  )

  marker
}
