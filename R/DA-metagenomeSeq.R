# We recommend fitFeatureModel over fitZig. MRcoefs, MRtable and MRfulltable
# are useful summary tables of the model outputs. We currently recommend using
# the zero-inflated log-normal model as implemented in fitFeatureModel.
#
# from biocore/qiime/blob/master/qiime/support_files/R/fitZIG.r
# /xia-lab/MicrobiomeAnalystR/blob/master/R/general_anal.R#L505
# https://support.bioconductor.org/p/78230/

# Difference between fitFeatureModel and fitZIG in metagenomeSeq,
# https://support.bioconductor.org/p/94138/.
#
# fitFeatureModel doesn't seem to allow for multiple comparisons.

#' metagenomeSeq differential analysis
#'
#' Differential expression analysis based on the Zero-inflated Log-Normal
#' mixture model or Zero-inflated Gaussian mixture model using metagenomeSeq.
#'
#' @param ps  ps a [`phyloseq::phyloseq-class`] object.
#' @param group  character, the variable to set the group, must be one of
#'   the var of the sample metadata.
#' @param confounders character vector, the confounding variables to be adjusted.
#'   default `character(0)`, indicating no confounding variable.
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of `phyloseq::rank_names(ps)`,
#'   or "all" means to summarize the taxa by the top taxa ranks
#'   (`summarize_taxa(ps, level = rank_names(ps)[1])`), or "none" means perform
#'   differential analysis on the original taxa (`taxa_names(ps)`, e.g.,
#'   OTU or ASV).
#' @param contrast this parameter only used for two groups comparison while
#'   there are multiple groups. For more please see the following details.
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
#' @param norm_para arguments passed to specific normalization methods.
#' @param method character, which model used for differential analysis,
#'   "ZILN" (Zero-inflated Log-Normal mixture model)" or "ZIG" (Zero-inflated
#'    Gaussian mixture model). And the zero-inflated log-normal model is
#'    preferred due to the high sensitivity and low FDR.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param pvalue_cutoff numeric, p value cutoff, default 0.05
#' @param ... extra arguments passed to the model. more details see
#'   [`metagenomeSeq::fitFeatureModel()`] and [`metagenomeSeq::fitZig()`],
#'   e.g. `control` (can be setted using [`metagenomeSeq::zigControl()`]) for
#'   [`metagenomeSeq::fitZig()`].
#'
#' @details
#' metagnomeSeq provides two differential analysis methods, zero-inflated
#' log-normal mixture model (implemented in
#' [`metagenomeSeq::fitFeatureModel()`]) and zero-inflated Gaussian mixture
#' model (implemented in [`metagenomeSeq::fitZig()`]). We recommend
#' fitFeatureModel over fitZig due to high sensitivity and low FDR. Both
#' [`metagenomeSeq::fitFeatureModel()`] and [`metagenomeSeq::fitZig()`] require
#' the abundance profiles before normalization.
#'
#' For [`metagenomeSeq::fitZig()`], the output column is the coefficient of
#' interest, and logFC column in the output of
#' [`metagenomeSeq::fitFeatureModel()`] is analogous to coefficient. Thus,
#' logFC is really just the estimate the coefficient of interest in
#' [`metagenomeSeq::fitFeatureModel()`]. For more details see
#' these question [Difference between fitFeatureModel and fitZIG
#' in metagenomeSeq](https://support.bioconductor.org/p/94138/).
#'
#' `contrast` must be a two length character or `NULL` (default). It is only
#' required to set manually for two groups comparison when there are multiple
#' groups. The order determines the direction of comparison, the first element
#' is used to specify the reference group (control). This means that, the first
#' element is the denominator for the fold change, and the second element is
#' used as baseline (numerator for fold change). Otherwise, users do required
#' to concern this paramerter (set as default `NULL`), and if there are
#' two groups, the first level of groups will set as the reference group; if
#' there are multiple groups, it will perform an ANOVA-like testing to find
#' markers which difference in any of the groups.
#'
#' Of note, [`metagenomeSeq::fitFeatureModel()`] is not allows for multiple
#' groups comparison.
#'
#' @return  a [`microbiomeMarker-class`] object.
#' @export
#' @author Yang Cao
#' @importFrom stats model.matrix
#' @importFrom metagenomeSeq normFactors<- MRcounts
#' @importFrom Biobase pData<- pData
#' @references
#' Paulson, Joseph N., et al. "Differential abundance analysis for microbial
#' marker-gene surveys." Nature methods 10.12 (2013): 1200-1202.
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")
#' )
#' run_metagenomeseq(ps, group = "Enterotype")
run_metagenomeseq <- function(ps,
    group,
    confounders =  character(0),
    contrast = NULL,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "CSS",
    norm_para = list(),
    method = c("ZILN", "ZIG"),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    ...) {
    ps <- check_rank_names(ps) %>% 
        check_taxa_rank(taxa_rank)
    
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    method <- match.arg(method, c("ZILN", "ZIG"))
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )
    
    if (length(confounders)) {
        confounders <- check_confounder(ps, group, confounders)
    }
    
    meta <- sample_data(ps)
    groups <- meta[[group]]
    groups <- make.names(groups)
    if (!is.null(contrast)) {
        contrast <- make.names(contrast)
    }
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    groups <- set_lvl(groups, contrast)
    lvl <- levels(groups)
    n_lvl <- length(lvl)

    if (n_lvl > 2 && method == "ZILN") {
        stop(
            "ZILN method do not allows for multiple groups comparison,\n",
            "please try set method = `ZIG`",
            call. = FALSE
        )
    }

    # contrast_new <- create_contrast(groups, contrast)
    # When running fitZig by default there is an additional covariate added to
    # the design matrix (scalingFactor), add var scalingFactor (set as zero)
    # if (n_lvl > 2) {
    #     old_contrast_nms <- row.names(contrast_new)
    #     contrast_new <- rbind(contrast_new, rep(0, ncol(contrast_new)))
    #     # row names of contrasts consistent with of coefficients
    #     # otherwise, warning: row names of contrasts don't match col names of
    #     # coefficients in the following `contrast.fit()`
    #     row.names(contrast_new) <- c(
    #         paste0("groups", old_contrast_nms),
    #         "scalingFactor"
    #     )
    # }

    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # normalization, write a function here
    # fitZig fitFeatureModel
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)

    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    mgs_summarized <- phyloseq2metagenomeSeq(ps_summarized)

    # extract norm factors and set the norm factors of MRexperiment
    nf <- get_norm_factors(ps_normed)
    if (!is.null(nf)) {
        pData(mgs_summarized@expSummary$expSummary)$normFactors <- nf
    } else {
        # for TSS, CRL and rarefy: normalized the feature table using CSS method
        ct <- metagenomeSeq::MRcounts(mgs_summarized, norm = FALSE)
        fun_p <- select_quantile_func(ct)
        mgs_summarized <- metagenomeSeq::cumNorm(
            mgs_summarized,
            p = fun_p(mgs_summarized)
        )
    }

    sl <- ifelse("sl" %in% names(norm_para), norm_para[["sl"]], 1000)
    counts_normalized <- metagenomeSeq::MRcounts(
        mgs_summarized,
        norm = TRUE,
        sl = sl
    )

    # mod <- model.matrix(~ 0 + groups)
    mod <- create_design(groups, meta, confounders)

    if (n_lvl == 2) {
        if (method == "ZILN") {
            tryCatch(
                fit <- metagenomeSeq::fitFeatureModel(mgs_summarized, mod, ...),
                error = function(e) {
                    paste0(
                        "fitFeatureModel model failed to fit to your data! ",
                        "Consider fitZig model or further ",
                        "filtering your dataset!"
                    )
                }
            )
        } else {
            tryCatch(
                fit <- metagenomeSeq::fitZig(mgs_summarized, mod, ...),
                error = function(e) {
                    paste0(
                        "fitZig model failed to fit to your data! ",
                        "Consider fitFeatureModel model or further ",
                        "filtering your dataset!"
                    )
                }
            )
        }

        # metagenomeSeq vignette: We recommend the user remove features based on
        # the number of estimated effective samples, please see
        # calculateEffectiveSamples. We recommend removing features with less 
        # than the average number of effective samples in all features. In 
        # essence, setting eff = .5 when using MRcoefs, MRfulltable, or MRtable.
        res <- metagenomeSeq::MRcoefs(
            fit,
            number = ntaxa(ps_summarized),
            adjustMethod = p_adjust,
            group = 3,
            eff = 0.5
        )
        res <- dplyr::rename(
            res, 
            pvalue = .data$pvalues, 
            padj = .data$adjPvalues
        )

        # For fitZig, the output var is the coefficient of interest (effect size
        # ), For fitFeaturemodel, logFC is anologous to coefficient of fitZig
        # (as logFC is really just the estimate the coefficient of interest).
        # Thus, we change the var of coefficent of interest to logFC for fitZig
        # https://support.bioconductor.org/p/94138/
        if (method == "ZIG") {
            names(res)[2] <- "logFC"
        }
        ef_var <- "logFC"
        res$enrich_group <- ifelse(res[[ef_var]] > 0, lvl[2], lvl[1])
    } else {
        fit <- metagenomeSeq::fitZig(mgs_summarized, mod, ...)
        zigfit <- slot(fit, "fit")
        # warning: row names of contrasts don't match col names of coefficients
        para_cf <- calc_coef(groups, mod, contrast)
        new_fit <- limma::contrasts.fit(zigfit, coefficients = para_cf)
        new_fit <- limma::eBayes(new_fit)
        res <- limma::topTable(
            new_fit,
            number = Inf,
            adjust.method = p_adjust,
        )
        res <- dplyr::filter(res, .data$adj.P.Val <= pvalue_cutoff) %>%
            dplyr::rename(pvalue = .data$P.Value, padj = .data$adj.P.Val)

        ef_var <- ifelse(is.null(contrast), "F", "logFC")

        # enrich group
        if (ef_var == "logFC") {
            # exp_lvl <- lvl[contrast_new == 1]
            # ref_lvl <- lvl[contrast_new == -1]
            # enrich_group <- ifelse(res$logFC > 0, exp_lvl, ref_lvl)
            enrich_group <- ifelse(res$logFC > 0, lvl[2], lvl[1])
        } else {
            cf <- zigfit$coefficients
            target_idx <- grepl("group", colnames(cf))
            cf <- cf[, target_idx]
            cf <- cbind(0, cf)
            enrich_group <- lvl[apply(cf, 1, which.max)]
            # sort the enrich_group according to the DE of topTags
            enrich_group <- enrich_group[match(row.names(res), row.names(cf))]
        }
        res$enrich_group <- enrich_group
    }
    res <- cbind(feature = row.names(res), res)
    res <- res[, c(
        "feature", "enrich_group",
        ef_var, "pvalue", "padj"
    )]
    row.names(res) <- paste0("marker", seq_len(nrow(res)))

    # rename the ef
    names(res)[3] <- ifelse(
        ef_var %in% c("logFC", "F"),
        paste0("ef_", ef_var),
        paste0("ef_", "coef")
    )
    sig_res <- res[res$padj < pvalue_cutoff & !is.na(res$padj), ]
    marker <- return_marker(sig_res, res)
    marker <- microbiomeMarker(
        marker_table = marker,
        norm_method = get_norm_method(norm),
        diff_method = paste0("metagenomeSeq: ", method),
        otu_table = otu_table(counts_normalized, taxa_are_rows = TRUE),
        sam_data = sample_data(ps_normed),
        tax_table = tax_table(ps_summarized)
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
#' @seealso [`metagenomeSeq::fitTimeSeries()`],
#'   [`metagenomeSeq::fitLogNormal()`],[`metagenomeSeq::fitZig()`],
#'   [`metagenomeSeq::MRtable()`],[`metagenomeSeq::MRfulltable()`]
#' @export
#' @importFrom Biobase AnnotatedDataFrame
#' @importMethodsFrom phyloseq t
#' @examples
#' data(caporaso)
#' phyloseq2metagenomeSeq(caporaso)
phyloseq2metagenomeSeq <- function(ps, ...) {
    # Enforce orientation. Samples are columns
    if (!taxa_are_rows(ps)) {
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

    # Create MRexperiment
    mr_obj <- metagenomeSeq::newMRexperiment(
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
    sdf <- sample_data(data.frame(sample = paste0("sam", seq_len(ncol(ps)))))
    row.names(sdf) <- colnames(ps)

    ps <- phyloseq(
        ps,
        sdf
    )
    mgs <- phyloseq2metagenomeSeq(ps)

    mgs
}


# get enrich group of a feature for multiple groups comparison
# group_pairs and logFC_pairs are the same length
get_mgs_enrich_group <- function(group_pairs, logFC_pairs) {
    all_groups <- unique(unlist(group_pairs))
    for (i in seq_along(group_pairs)) {
        group_low <- ifelse(
            logFC_pairs[i] > 0,
            group_pairs[[i]][2],
            group_pairs[[i]][1]
        )
        all_groups <- setdiff(all_groups, group_low)
    }

    all_groups
}
