#' Differential analysis using limma-voom
#'
#' @param ps  ps a [`phyloseq::phyloseq-class`] object.
#' @param group  character, the variable to set the group, must be one of
#'   the var of the sample metadata.
#' @param confounders character vector, the confounding variables to be adjusted.
#'   default `character(0)`, indicating no confounding variable.
#' @param contrast this parameter only used for two groups comparison while
#'   there are multiple groups. For more please see the following details.
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of
#'   `phyloseq::rank_names(phyloseq)`, or "all" means to summarize the taxa by
#'   the top taxa ranks (`summarize_taxa(ps, level = rank_names(ps)[1])`), or
#'   "none" means perform differential analysis on the original taxa
#'   (`taxa_names(phyloseq)`, e.g., OTU or ASV).
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
#' @param norm_para arguments passed to specific normalization methods. Most
#'   users will not need to pass any additional arguments here.
#' @param voom_span width of the smoothing window used for the lowess
#'   mean-variance trend for [`limma::voom()`]. Expressed as a proportion
#'   between 0 and 1.
#' @param p_adjust method for multiple test correction, default `none`,
#'   for more details see [stats::p.adjust].
#' @param pvalue_cutoff cutoff of p value, default 0.05.
#' @param ... extra arguments passed to [`limma::eBayes()`].
#'
#' @export
#' @return a [`microbiomeMarker-class`] object.
#' @references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014).
#'   voom: Precision weights unlock linear model analysis tools for RNA-seq read
#'   counts. Genome biology, 15(2), 1-17.
#'
#' @details
#' `contrast` must be a two length character or `NULL` (default). It is only
#' required to set manually for two groups comparison when there are multiple
#' groups. The order determines the direction of comparison, the first element
#' is used to specify the reference group (control). This means that, the first
#' element is the denominator for the fold change, and the second element is
#' used as baseline (numerator for fold change). Otherwise, users do required
#' to concern this parameter (set as default `NULL`), and if there are
#' two groups, the first level of groups will set as the reference group; if
#' there are multiple groups, it will perform an ANOVA-like testing to find
#' markers which difference in any of the groups.
#'
#' @examples
#' data(enterotypes_arumugam)
#' mm <- run_limma_voom(
#'     enterotypes_arumugam,
#'     "Enterotype",
#'     contrast = c("Enterotype 3", "Enterotype 2"),
#'     pvalue_cutoff = 0.01,
#'     p_adjust = "none"
#' )
#' mm
run_limma_voom <- function(ps,
    group,
    confounders =  character(0),
    contrast = NULL,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "none",
    norm_para = list(),
    voom_span = 0.5,
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    ...) {
    stopifnot(inherits(ps, "phyloseq"))
    ps <- check_rank_names(ps) %>% 
        check_taxa_rank(taxa_rank)
    
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
    # meta_nms <- names(meta)
    # if (!group %in% meta_nms) {
    #     stop(
    #         group, " are not contained in the `sample_data` of `ps`",
    #         call. = FALSE
    #     )
    # }
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

    transform <- match.arg(transform, c("identity", "log10", "log10p"))

    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # normalize the data
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    counts <- abundances(ps_summarized, norm = FALSE)
    # row.names(counts) <- tax_table(ps_summarized)[, 1]

    # design matrix
    design <- create_design(groups, meta, confounders)

    # library size
    nf <- get_norm_factors(ps_normed)
    lib_size <- phyloseq::sample_sums(ps)
    if (!is.null(nf)) {
        lib_size <- nf * lib_size
    }

    voom_out <- limma::voom(
        counts,
        design = design,
        lib.size = lib_size,
        span = voom_span
    )
    fit_out <- limma::lmFit(voom_out, design = design)
    
    para_cf <- calc_coef(groups, design, contrast)
    # fit_out <- limma::contrasts.fit(fit_out, coefficients = para_cf)
    
    # if (length(contrast_new) == n_lvl) {
    #     # warning: row names of contrasts don't match col names of coefficients
    #     fit_out <- limma::contrasts.fit(fit_out, contrast_new)
    # }
    test_out <- limma::eBayes(fit_out, ...)
    test_df <- limma::topTable(
        test_out,
        coef = para_cf,
        number = nrow(counts),
        adjust.method = p_adjust
    )

    counts_normed <- abundances(ps_summarized, norm = TRUE)
    if (n_lvl == 2 || !is.null(contrast)) {
        enrich_group <- ifelse(test_df$logFC > 0, lvl[2], lvl[1])
        
        ef <- test_df[["logFC"]]
        ef_name <- "ef_logFC"
    } else {
        cf <- fit_out$coefficients
        target_idx <- grepl("group", colnames(cf))
        cf <- cf[, target_idx]
        cf <- cbind(0, cf)
        enrich_group <- lvl[apply(cf, 1, which.max)]
        enrich_group <- enrich_group[match(row.names(test_df), row.names(cf))]
        
        ef <- test_df[["F"]]
        ef_name <- "ef_F_statistic"
    }
    # if (length(contrast_new) == n_lvl) {
    #     exp_lvl <- lvl[contrast_new == 1]
    #     ref_lvl <- lvl[contrast_new == -1]
    #     enrich_group <- ifelse(test_df$logFC > 0, exp_lvl, ref_lvl)
    # } else {
    #     cf <- fit_out$coefficients
    #     enrich_idx <- apply(cf, 1, which.max)
    #     enrich_group <- lvl[enrich_idx]
    #     enrich_group <- enrich_group[match(row.names(test_df), row.names(cf))]
    # }

    # if (length(contrast_new) == n_lvl) {
    #     ef <- test_df[["logFC"]]
    #     ef_name <- "ef_logFC"
    # } else {
    #     ef <- test_df[["F"]]
    #     ef_name <- "ef_F_statistic"
    # }

    marker <- data.frame(
        feature = row.names(test_df),
        enrich_group = enrich_group,
        ef = ef,
        pvalue = test_df$P.Value,
        padj = test_df$adj.P.Val
    )
    names(marker)[3] <- ef_name
    sig_marker <- dplyr::filter(marker, .data$padj <= pvalue_cutoff)
    marker <- return_marker(sig_marker, marker)
    mm <- microbiomeMarker(
        marker_table = marker,
        norm_method = get_norm_method(norm),
        diff_method = "limma_voom",
        sam_data = sample_data(ps_normed),
        otu_table = otu_table(counts_normed, taxa_are_rows = TRUE),
        tax_table = tax_table(ps_summarized)
    )

    mm
}
