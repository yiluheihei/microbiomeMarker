#' Perform differential analysis using edgeR
#'
#' Differential expression analysis based on the Negative Binomial distribution
#' using **edgeR**.
#'
#' @param ps  ps a [`phyloseq::phyloseq-class`] object.
#' @param group_var  character, the variable to set the group, must be one of
#'   the var of the sample metadata
#' @param contrast  a two length vector used for two groups comparison. The
#'   order determines the direction of fold change, the first element is the
#'   numerator for the fold change, and the second element is used as baseline
#'   (denominator for fold change); or `NULL` for multiple groups comparison.
#' @param method character, used for differential analysis, please see details
#'   below for more info.
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation
#'     (default).
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
#'     The scaling factor is then derived using a weighted trimmed mean over
#'     the differences of the log-transformed gene-count fold-change between
#'     the sample and the reference.
#'   * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#'   * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#'   * "CLR": centered log-ratio normalization.
#'   * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' @param norm_para arguments passed to specific normalization methods. Most
#'   users will not need to pass any additional arguments here.
#' @param disp_para additional arguments passed to [`edgeR::estimateDisp()`]
#'   used for dispersions estimation. Most users will not need to pass any
#'   additional arguments here.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param pvalue_cutoff numeric, p value cutoff, default 0.05
#' @param ... extra arguments passed to the model. See [`edgeR::glmQLFit()`]
#'   and [`edgeR::glmFit()`] for more details.
#' @return  a [`microbiomeMarker-class`] object.
#'
#' @details
#' **Note** that edgeR is designed to work with actual counts. This means that
#' transformation is not required in any way before inputting them to edgeR.
#'
#' There are two test methods for differential analysis in **edgeR**,
#' likelihood ratio test (LRT) and quasi-likelihood F-test (QLFT). The QLFT
#' method is recommended as it allows stricter error rate control by
#' accounting for the uncertainty in dispersion estimation.
#'
#' The para `contrast` is used to specify contrast of the linear model
#' coefficients to be tested equal to zero. `NULL` means performing multiple
#' groups comparisons. This is done by creating a matrix of all
#' pairwise independent contrasts. In this manner, users can
#  perform a ANOVA-like test to find features that are differential between
#  any of the groups.
#'
#' @export
#' @seealso [`edgeR::glmFit()`],[`edgeR::glmQLFit()`],[`edgeR::estimateDisp()`]
#'   ,[`normalize()`]
#' @author Yang Cao
#' @references
#' Robinson, Mark D., and Alicia Oshlack. "A scaling normalization method for
#' differential expression analysis of RNA-seq data." Genome biology 11.3
#' (2010): 1-9.
#'
#' Robinson, Mark D., Davis J. McCarthy, and Gordon K. Smyth. "edgeR: a
#' Bioconductor package for differential expression analysis of digital
#' gene expression data." Bioinformatics 26.1 (2010): 139-140.
run_edger <- function(ps,
                      group_var,
                      contrast = NULL,
                      method = c("LRT", "QLFT"),
                      transform = c("identity", "log10", "log10p"),
                      norm = "TMM",
                      norm_para = list(),
                      disp_para = list(),
                      p_adjust = c("none", "fdr", "bonferroni", "holm",
                                   "hochberg", "hommel", "BH", "BY"),
                      pvalue_cutoff = 0.05,
                      ...) {
  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  method <- match.arg(method, c("LRT", "QLFT"))
  p_adjust <- match.arg(
    p_adjust,
    c("none", "fdr", "bonferroni", "holm",
      "hochberg", "hommel", "BH", "BY")
  )

  groups <- sample_data(ps)[[group_var]]
  groups <- factor(groups)
  lvl_n <- nlevels(groups)

  if (lvl_n < 2) {
    stop("Differential analysis requires at least two groups.")
  }

  # For two-groups comparison, contrast must be a two-length vector,
  # for multiple groups comparison, contrast must be NULL
  if (!is.null(contrast)) {
    if (length(contrast) != 2) {
      stop(
        "`contrast` must be length 2 or NULL (for multiple groups comparison).",
        call. = FALSE
      )
    }

    for (cont in contrast) {
      if (! cont %in% groups) {
        stop(
          "The element of `contrast` should be one of ",
          paste(unique(groups), collapse = ", "),
          call. = FALSE
        )
      }
    }
    two_groups <- contrast
    idx <- match(contrast, levels(groups))
    contrast <- rep(0, lvl_n)
    contrast[idx[1]] <- 1
    contrast[idx[2]] <- -1
  } else {
    if (lvl_n > 2) {
      contrast <- generate_anova_contrast(levels(groups))
    } else {
      stop(
        "`contrast` must be a two length vector for two groups comparison.",
        call. = FALSE
      )
    }
  }

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

  design <- stats::model.matrix(~0+groups)
  disp_para <- c(disp_para, y = list(dge_summarized), design = list(design))
  dge_summarized <- do.call(edgeR::estimateDisp, disp_para)

  # differential expression
  # quasi-likelihood (QL) F-test is used as it reflects the uncertainty in
  #  estimating the dispersion for each feature, and gives stricter error
  #  rate control
  fit_fun <- ifelse(method == "LRT", edgeR::glmFit, edgeR::glmQLFit)
  test_fun <- ifelse(method == "LRT", edgeR::glmLRT, edgeR::glmQLFTest)
  fit <- fit_fun(dge_summarized, design, ...)
  lrt <- test_fun(fit, contrast = contrast)
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

  # normalized counts
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

  # enrich group
  if (!is.matrix(contrast)) {
    enrich_group <- ifelse(res$logFC > 0, two_groups[1], two_groups[2])
  } else {
    enrich_group <- get_sl_enrich_group(counts_normalized, groups)
  }
  res$enrich_group <- enrich_group

  res_filtered <- res[res$padj < pvalue_cutoff & !is.na(res$padj), ]
  # edgeR::decideTestsDGE(), dentify which genes are significantly
  # differentially expressed from an edgeR fit object containing p-values and
  # test statistics.

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
  # var of effect size: named as ef_<name of the actual effect size>,
  # e.g. ef_F, ef_LR
  ef_name <- ifelse(method == "LRT", "logFC", "F")


  # only keep five variables: feature, enrich_group, effect_size (LR for LRT
  # F for QLFT), pvalue, and padj, write a function? select_marker_var
  # (effect_size = "")
  keep_var <- c("feature", "enrich_group", ef_name, "pvalue", "padj")
  sig_feature <- sig_feature[keep_var]
  names(sig_feature)[3] <- paste0("ef_", ef_name)

  marker <- microbiomeMarker(
    marker_table = marker_table(sig_feature),
    norm_method = get_norm_method(norm),
    diff_method = paste("edgeR:", method),
    sam_data = sample_data(ps_normed),
    # tax_table = tax_table(ps),
    otu_table = otu_table(counts_normalized, taxa_are_rows = TRUE),
    tax_table = tax_table(ps_summarized)
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


# create all pair-wise comparisons (contrasts) for anova-like test
generate_anova_contrast <- function(levels) {
  n <- length(levels)
  design <- matrix(0, n, choose(n, 2))
  rownames(design) <- levels
  colnames(design) <- seq_len(choose(n, 2))
  k <- 0
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      k <- k + 1
      design[j, k] <- 1
      design[i, k] <- -1
      colnames(design)[k] <- paste(levels[j], "-", levels[i], sep = "")
    }
  }
  design
}
