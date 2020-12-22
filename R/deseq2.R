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
# https://github.com/biocore/qiime/blob/master/qiime/support_files/R/DESeq2_nbinom.r

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
#' @param test,fitType,sfType,betaPrior,modelMatrixType,useT,minmu these seven
#'   parameters are inherited form [`DESeq2::DESeq()`].
#'   - `test`, should be either "Wald" or "LRT", which will
#'     then use either Wald significance tests, or the likelihood ratio test on
#'     the difference in deviance between a full and reduced model formula.
#'   - `fitType`, either "parametric", "local", "mean", or "glmGamPoi" for the
#'     type of fitting of dispersions to the mean intensity.
#'   - `sfType`, either "ratio", "poscounts", or "iterate" for the type of size
#'     factor estimation.
#'   - `betaPrior`, whether or not to put a zero-mean normal prior on the
#'     non-intercept coefficients.
#'   - `modelMatrixType`, either "standard" or "expanded", which describe how
#'     the model matrix,
#'   - `useT`, logical, where Wald statistics are assumed to follow a standard
#'     Normal.
#'   - `minmu`, lower bound on the estimated count for fitting gene-wise
#'     dispersion.
#'
#'   For more details, see [`DESeq2::DESeq()`].
#'
#' @param p_adjust method for multiple test correction, default `none`, for
#'   more details see [stats::p.adjust].
#' @param p_value_cutoff p_value_cutoff numeric, p value cutoff, default 0.05.
#' @param ... extra parameters passed to [`DESeq2::DESeq()`].
#' @export
#' @return a [`microbiomeMarker-class`] object.
#' @seealso [`DESeq2::results()`],[`DESeq2::DESeq()`]
#' @importFrom stats formula
#' @importFrom DESeq2 dispersions<-
#' @importMethodsFrom S4Vectors mcols
run_deseq2 <- function(ps,
                      group_var,
                      subgroup1,
                      subgroup2,
                      transform = c("identity", "log10", "log10p"),
                      test = c("Wald", "LRT"),
                      fitType = c("parametric", "local", "mean", "glmGamPoi"),
                      sfType = "poscounts",
                      betaPrior =FALSE,
                      modelMatrixType,
                      useT = FALSE,
                      minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
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

  dsg <- formula(paste("~ ", group_var))
  suppressWarnings(
    dds <- phyloseq::phyloseq_to_deseq2(ps_summarized, design = dsg)
  )

  # error: all gene-wise dispersion estimates are within 2 orders of magnitude
  # from the minimum value, which indicates that the count are not overdispersed
  #
  # If dispersion values are less than 1e-6  (minimal value is 1e-8),
  # it would be problematic to fit a dispersion trend in DESeq2.
  # The reason for a minimal value, is that for a given row of the count matrix,
  # the maximum likelihood estimate can tend to 0 (and so we have a rule to stop
  # after 1e-8)
  # https://support.bioconductor.org/p/63845/
  # https://support.bioconductor.org/p/122757/
  # https://github.com/biocore/qiime/blob/master/qiime/support_files/R/DESeq2_nbinom.r

  res_deseq <- try(
    DESeq2::DESeq(
      dds,
      test = test,
      fitType = fitType,
      sfType = sfType,
      quiet = TRUE,
      betaPrior = betaPrior,
      modelMatrixType = modelMatrixType,
      useT = useT,
      minmu = minmu,
      ...),
    silent = TRUE
  )

  if (inherits(res_deseq, "try-error") && fitType != "local") {
    warning("data is not overdispered, try `fitType = 'local'`")
    res_deseq <- try(
      DESeq2::DESeq(
        dds,
        test = test,
        fitType = "local",
        sfType = sfType,
        quiet = TRUE,
        betaPrior = betaPrior,
        modelMatrixType = modelMatrixType,
        useT = useT,
        minmu = minmu,
        ...),
      silent = TRUE
    )
  }
  if (inherits(res_deseq, "try-error") && fitType != "mean") {
    warning("data is not overdispered, try `fitType = 'mean'`")
    res_deseq <- try(
      DESeq2::DESeq(
        dds,
        test = test,
        fitType = "mean",
        sfType = sfType,
        quiet = TRUE,
        betaPrior = betaPrior,
        modelMatrixType = modelMatrixType,
        useT = useT,
        minmu = minmu,
        ...),
      silent = TRUE
    )
  }
  if (inherits(res_deseq, "try-error")) {
    warning("data is not overdispered, use gene-wise estimates as final estimates")
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- DESeq2::estimateDispersionsGeneEst(dds)
    DESeq2::dispersions(dds) <- mcols(dds)$dispGeneEst

    dds <- DESeq2::nbinomWaldTest(
      dds,
      betaPrior = betaPrior,
      quiet = TRUE,
      modelMatrixType = modelMatrixType,
      useT = useT,
      minmu = minmu
    )
  } else {
    dds <- res_deseq
  }

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


#' # used if data is not overdispered
#' # this function is modified from the DESeq2::DESeq()
#' # https://github.com/mikelove/DESeq2/blob/master/R/core.R#L271
#' #
#' #' @importMethodsFrom DESeq2 design
#' deseq2 <- function(object,
#'                    test = c("Wald", "LRT"),
#'                    sfType = "poscounts",
#'                    type = c("DESeq2", "glmGamPoi"),
#'                    betaPrior =FALSE,
#'                    full,
#'                    reduced,
#'                    minReplicatesForReplace = 7,
#'                    modelMatrixType,
#'                    useT = FALSE,
#'                    minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5) {
#'   stopifnot(is(object, "DESeqDataSet"))
#'   stopifnot(is.logical(betaPrior))
#'   test <- match.arg(test, choices=c("Wald","LRT"))
#'
#'   sfType <- match.arg(sfType, choices=c("ratio","poscounts","iterate"))
#'
#'   stopifnot(is.numeric(minReplicatesForReplace))
#'   modelAsFormula <- !is.matrix(full) & is(design(object), "formula")
#'
#'   if (test == "LRT") {
#'     if (missing(reduced)) {
#'       stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
#'     }
#'     if (betaPrior) {
#'       stop(
#'         paste(
#'           "test = 'LRT' does not support use of LFC shrinkage,",
#'           "use betaPrior=FALSE"
#'         )
#'       )
#'     }
#'     if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
#'       stop("test='LRT' does not support use of expanded model matrix")
#'     }
#'     if (is.matrix(full) | is.matrix(reduced)) {
#'       if (!(is.matrix(full) & is.matrix(reduced))) {
#'         stop(
#'           paste(
#'             "if one of 'full' and 'reduced' is a matrix,",
#'             "the other must be also a matrix")
#'         )
#'       }
#'     }
#'
#'     if (modelAsFormula) {
#'       DESeq2:::checkLRT(full, reduced)
#'     } else {
#'       DESeq2:::checkFullRank(full)
#'       DESeq2:::checkFullRank(reduced)
#'       if (ncol(full) <= ncol(reduced)) {
#'         stop(
#'           paste("the number of columns of 'full' should be more than",
#'                 "the number of columns of 'reduced'")
#'         )
#'       }
#'     }
#'   }
#'
#'   if (test == "Wald" & !missing(reduced)) {
#'     warning("'reduced' ignored when test='Wald'")
#'   }
#'
#'   if (modelAsFormula) {
#'     # run some tests common to DESeq, nbinomWaldTest, nbinomLRT
#'     DESeq2:::designAndArgChecker(object, betaPrior)
#'
#'     if (design(object) == formula(~1)) {
#'       warning("the design is ~ 1 (just an intercept). is this intended?")
#'     }
#'
#'     if (full != design(object)) {
#'       stop("'full' specified as formula should equal design(object)")
#'     }
#'     modelMatrix <- NULL
#'   } else {
#'     # model not as formula, so DESeq() is using supplied model matrix
#'     if (!quiet) message("using supplied model matrix")
#'     if (betaPrior == TRUE) {
#'       stop("betaPrior=TRUE is not supported for user-provided model matrices")
#'     }
#'     DESeq2:::checkFullRank(full)
#'     # this will be used for dispersion estimation and testing
#'     modelMatrix <- full
#'   }
#'
#'   attr(object, "betaPrior") <- betaPrior
#'
#'   if (test == "Wald") {
#'     object <- DESeq2::nbinomWaldTest(
#'       object,
#'       betaPrior=betaPrior,
#'       quiet = TRUE,
#'       modelMatrix = modelMatrix,
#'       modelMatrixType = modelMatrixType,
#'       useT=useT,
#'       minmu=minmu
#'     )
#'   } else if (test == "LRT") {
#'     object <- DESeq2::nbinomLRT(
#'       object,
#'       full = full,
#'       reduced = reduced,
#'       quiet = TRUE,
#'       minmu = minmu,
#'       type = type
#'     )
#'   }
#'
#'   sufficientReps <- any(DESeq2:::nOrMoreInCell(
#'     attr(object,"modelMatrix"),minReplicatesForReplace)
#'   )
#'   if (sufficientReps) {
#'     object <- DESeq2:::refitWithoutOutliers(
#'       object,
#'       test = test,
#'       betaPrior = betaPrior,
#'       full = full, reduced = reduced, quiet = TRUE,
#'       minReplicatesForReplace = minReplicatesForReplace,
#'       modelMatrix = modelMatrix,
#'       modelMatrixType = modelMatrixType
#'     )
#'   }
#'
#'   object
#' }
