#' Perform differential analysis using edgeR
#'
#' Differential expression analysis based on the Negative Binomial distribution
#' using **edgeR**.
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
#' @param norm_para arguments passed to specific normalization methods. Most
#'   users will not need to pass any additional arguments here.
#' @param disp_para additional arguments passed to [`edgeR::estimateDisp()`]
#'   used for dispersions estimation. Most users will not need to pass any
#'   additional arguments here.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param p_value_cutoff numeric, p value cutoff, default 0.05
#' @param ... extra arguments passed to the model. See [`edgeR::glmFit()`] for
#'   more details.
#' @return  a [`microbiomeMarker-class`] object.
#' @export
#' @seealso [`edgeR::glmFit()`],[`edgeR::estimateDisp()`],[`normalize()`]
#' @author Yang Cao
run_edger <- function(ps,
                      group_var,
                      subgroup1,
                      subgroup2,
                      transform = c("identity", "log10", "log10p"),
                      norm = "TMM",
                      norm_para = list(),
                      disp_para = list(),
                      p_adjust = c("none", "fdr", "bonferroni", "holm",
                                   "hochberg", "hommel", "BH", "BY"),
                      p_value_cutoff = 0.05,
                      ...) {
  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY")
  )

  # filter the samples in subgroup1 or subgroup2
  groups <- sample_data(ps)[[group_var]]
  levels(groups) <- c(subgroup1, subgroup2)
  ps <- phyloseq::prune_samples(groups %in% c(subgroup1, subgroup2), ps)

  # preprocess phyloseq object
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)

  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)

  # summarize data and  add norm.factors var to samples of DGEList
  ps_summarized <- summarize_taxa(ps_normed)
  dge_summarized <- phyloseq2edgeR(ps_summarized)

  nf <- get_norm_factors(ps_normed)
  if (!is.null(nf)) {
    dge_summarized$samples$norm.factors <- nf
  } else {
    # for TSS, CRL and rarefy (no norm factors is saved),
    # normalized the feature table using TMM method in edgeR
    # using the default arguments.
    dge_summarized <- edgeR::calcNormFactors(dge_summarized, method = "TMM")
  }

  # estimate dispersion
  design <- stats::model.matrix(~groups)
  disp_para <- c(disp_para, y = list(dge_summarized), design = list(design))
  dge_summarized <- do.call(edgeR::estimateDisp, disp_para)

  # differential expression
  fit <- edgeR::glmFit(
    dge_summarized,
    design,
    dispersion = NULL,
    ...
  )
  lrt <- edgeR::glmLRT(fit)
  res <- edgeR::topTags(
    lrt,
    n = ntaxa(ps_summarized),
    adjust.method = p_adjust,
    sort.by = "PValue"
  )
  res <- res$table
  if ("FDR" %in% names(res)) {
    res <- dplyr::rename(res, pvalue = .data$PValue, padj = .data$FDR)
  } else if ("FWER" %in% names(res)) {
    res <- dplyr::rename(res, pvalue = .data$PValue, padj = .data$FWER)
  } else {
    res <- dplyr::rename(res, pvalue = .data$PValue)
    res$padj <- res$pvalue
  }

  res$enrich_group <- ifelse(
    res$logFC > 0,
    subgroup2,
    subgroup1
  )

  res_filtered <- res[res$padj < p_value_cutoff & !is.na(res$padj), ]

  if (nrow(res_filtered) == 0) {
    warning("No significant features were found, return all the features")
    sig_feature <- cbind(feature = row.names(res), res)
  } else {
    sig_feature <- cbind(feature = row.names(res_filtered), res_filtered)
  }
  # first two columns: feature enrich_group (write a function)
  other_col <- setdiff(names(sig_feature), c("feature", "enrich_group"))
  sig_feature <- sig_feature[, c("feature", "enrich_group", other_col)]
  row.names(sig_feature) <- paste0("marker", seq_len(nrow(sig_feature)))

  ef_nf <- dge_summarized$samples$lib.size * dge_summarized$samples$norm.factors
  ref_nf <- mean(ef_nf)
  counts_normalized <-
    sweep(
      as(otu_table(ps_summarized), "matrix"),
      MARGIN = 2,
      ef_nf,
      "/"
    ) *
    ref_nf
  row.names(counts_normalized) <- row.names(tax_table(ps_summarized))

  marker <- microbiomeMarker(
    marker_table = marker_table(sig_feature),
    tax_table_orig = tax_table(ps),
    otu_table(counts_normalized, taxa_are_rows = TRUE),
    tax_table(ps_summarized)
  )

  marker
}

#' Convert phyloseq data to edgeR `DGEList` object
#'
#' This function convert [`phyloseq::phyloseq-class`] object to
#' [`edgeR::DGEList-class`] object, can then can be used to perform
#' differential analysis using the methods in **edgeR**.
#'
#' @param ps a [`phyloseq::phyloseq-class`] object.
#' @param ... optional, additional named arguments passed  to
#'   [`edgeR::DGEList()`]. Most users will not need to pass any additional
#'   arguments here.
#' @return A [`edgeR::DGEList-class`] object.
#' @export
phyloseq2edgeR <- function(ps, ...) {
  ps <- keep_taxa_in_rows(ps)
  abd <- as(otu_table(ps), "matrix")

  # tax_table: annotation information
  taxa <- tax_table(ps, FALSE)
  if (!is.null(taxa)) {
    taxa <- data.frame(as(taxa, "matrix"))
  }

  # sample_data: information on each sample
  samp <- sample_data(ps, FALSE)

  dge <- edgeR::DGEList(
    counts = abd,
    samples = samp,
    genes = taxa,
    ...
  )

  dge
}
