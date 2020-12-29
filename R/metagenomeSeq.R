# We recommend fitFeatureModel over fitZig. MRcoefs, MRtable and MRfulltable
# are useful summary tables of the model outputs. We currently recommend using
# the zero-inflated log-normal model as implemented in fitFeatureModel.
#
# In metagenomeSeq, both fitZig and fitFeatureModel require the abundance
# profiles before normalization, and the normalized counts can be obtained
# according to the norm factors of MRexperiment.
#
# https://github.com/biocore/qiime/blob/master/qiime/support_files/R/fitZIG.r
# https://github.com/xia-lab/MicrobiomeAnalystR/blob/master/R/general_anal.R#L505

#' metagenomeSeq differential analysis
#'
#' Differential expression analysis based on the Zero-inflated Log-Normal
#' mixture model or Zero-inflated Gaussian mixture model using metagenomeSeq.
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
#' @param norm_para arguments passed to specific normalization methods.
#' @param model character, which model used for differential analysis,
#'   "fitFeatureModel" (Zero-inflated Log-Normal mixture model)" or "fitZig"
#'   (Zero-inflated Gaussian mixture model). As the authors of **metagenomeSeq**
#'   said, we currently recommend using the zero-inflated log-normal model.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param p_value_cutoff numeric, p value cutoff, default 0.05
#' @param ... extra arguments passed to the model. more details see
#'   [`metagenomeSeq::fitFeatureModel()`] and [`metagenomeSeq::fitZig()`]
#'   for more details.
#' @return  a [`microbiomeMarker-class`] object.
#' @export
#' @author Yang Cao
#' @importFrom stats model.matrix
#' @importFrom metagenomeSeq normFactors<- MRcounts
#' @importFrom Biobase pData<- pData
run_metagenomeseq <- function(ps,
                              group_var,
                              subgroup1,
                              subgroup2,
                              transform = c("identity", "log10", "log10p"),
                              norm = "CSS",
                              norm_para = list(),
                              model = c("fitFeatureModel", "fitZig"),
                              p_adjust = c("none", "fdr", "bonferroni", "holm",
                                           "hochberg", "hommel", "BH", "BY"),
                              p_value_cutoff = 0.05,
                              ...) {

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  model <- match.arg(model, c("fitFeatureModel", "fitZig"))
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

  # normalization, write a function here
  # fitZig fitFeatureModel
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)

  # summarize data
  ps_summarized <- summarize_taxa(ps_normed)
  mgs_summarized <- phyloseq2metagenomeSeq(ps_summarized)

  # extract norm factors and set the norm factors of MRexperiment
  nf <- get_norm_factors(ps_normed)
  if (!is.null(nf)) {
    pData(mgs_summarized@expSummary$expSummary)$normFactors <- nf
  }

  mod <- model.matrix(~groups)

  if (model == "fitFeatureModel") {
    tryCatch(
      fit <- metagenomeSeq::fitFeatureModel(mgs_summarized, mod, ...),
      error = function(e) {
         paste0(
           "fitFeatureModel model failed to fit to your data! ",
           "Consider fitZig model or further filtering your dataset!"
         )
      }
    )
  } else {
    tryCatch(
        fit <- metagenomeSeq::fitZig(mgs_summarized, mod, ...),
      error = function(e) {
        paste0(
          "fitZig model failed to fit to your data! ",
          "Consider fitFeatureModel model or further filtering your dataset!"
        )
      }
    )
  }

  res <- metagenomeSeq::MRfulltable(
    fit,
    number = ntaxa(ps_summarized),
    adjustMethod = p_adjust,
    group = 3
  )
  res <- dplyr::rename(res, pvalue = .data$pvalues, padj = .data$adjPvalues)

  # # enrich_group
  # if (model == "fitFeatureModel") {
  #   res$enrich_group <- ifelse(res$logFC > 0, subgroup2, subgroup1)
  # } else {
  #   res$enrich_group <- ifelse(res$oddsRatio > 1, subgroup2, subgroup1)
  # }
  res$enrich_group <- ifelse(
    res$`counts in group 0` > res$`counts in group 1`,
    subgroup1,
    subgroup2
  )


  res_filtered <- res[res$padj < p_value_cutoff & !is.na(res$padj), ]

  # write a function
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

  marker <- microbiomeMarker(
    marker_table = marker_table(sig_feature),
    tax_table_orig = tax_table(ps_normed),
    otu_table(ps_summarized),
    tax_table(ps_summarized)
  )

  marker
}



# This function is modified from `phyloseq::phyloseq_to_metagenomeSeq()`,
# There two changes: 1) do not coerce count data to vanilla matrix of integers;
# 2) do not normalize the count.
#
#
#' Convert phyloseq data to MetagenomeSeq `MRexperiment` object
#'
#' The phyloseq data is converted to the relevant
#' [`metagenomeSeq::MRexperiment-class`] object, which can then be tested in
#' the zero-inflated mixture model framework in the metagenomeSeq package.
#'
#' @param ps [`phyloseq::phyloseq-class`] object for
#'   `phyloseq2metagenomeSeq()`, or [`phyloseq::otu_table-class`] object
#'   for `otu_table2metagenomeseq()`.
#' @param ... optional, additional named arguments passed  to
#'   [`metagenomeSeq::newMRexperiment()`]. Most users will not need to pass
#'   any additional arguments here.
#' @return A [`metagenomeSeq::MRexperiment-class`] object.
#' @seealso [`metagenomeSeq::fitTimeSeries()`],[`metagenomeSeq::fitLogNormal()`],
#'   [`metagenomeSeq::fitZig()`],[`metagenomeSeq::MRtable()`],
#'   [`metagenomeSeq::MRfulltable()`]
#' @export
#' @importFrom Biobase AnnotatedDataFrame
#' @importMethodsFrom phyloseq t
phyloseq2metagenomeSeq <- function(ps, ...) {
  # Enforce orientation. Samples are columns
  if (!taxa_are_rows(ps) ) {
    ps <- t(ps)
  }

  count <- as(otu_table(ps), "matrix")
  # Create sample annotation if possible
  if (!is.null(sample_data(ps, FALSE))) {
    adf <- AnnotatedDataFrame(data.frame(sample_data(ps)))
  } else {
    adf <- NULL
  }

  # Create taxa annotation if possible
  if (!is.null(tax_table(ps, FALSE))) {
    tdf <- AnnotatedDataFrame(
      data.frame(
        OTUname = taxa_names(ps),
        data.frame(tax_table(ps)),
        row.names = taxa_names(ps)
      )
    )
  } else {
    tdf <- AnnotatedDataFrame(
      data.frame(
        OTUname = taxa_names(ps),
        row.names = taxa_names(ps)
      )
    )
  }

  # setting the norm factor, or the fitzig or
  # nf <- sample_data(ps)[["metagenomeSeq_norm_factor"]]

  # Create MRexperiment
  mr_obj = metagenomeSeq::newMRexperiment(
    counts = count,
    phenoData = adf,
    featureData = tdf,
    ...
  )

  mr_obj
}


#' @rdname phyloseq2metagenomeSeq
#' @export
otu_table2metagenomeSeq <- function(ps, ...) {
  stopifnot(inherits(ps, "otu_table"))
  # create a sample data with only one var "sample": sam1, sam2
  sdf <- sample_data(data.frame(sample = paste0("sam", 1:ncol(ps))))
  row.names(sdf) <- colnames(ps)

  ps <- phyloseq(
    ps,
    sdf
  )
  mgs <- phyloseq2metagenomeSeq(ps)

  mgs
}

