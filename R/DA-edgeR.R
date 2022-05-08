#' Perform differential analysis using edgeR
#'
#' Differential expression analysis based on the Negative Binomial distribution
#' using **edgeR**.
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
#'   for more details see [stats::p.adjust].
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
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")
#' )
#' run_edger(ps, group = "Enterotype")
run_edger <- function(ps,
    group,
    confounders =  character(0),
    contrast = NULL,
    taxa_rank = "all",
    method = c("LRT", "QLFT"),
    transform = c("identity", "log10", "log10p"),
    norm = "TMM",
    norm_para = list(),
    disp_para = list(),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    ...) {
    ps <- check_rank_names(ps) %>% 
        check_taxa_rank( taxa_rank)
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    method <- match.arg(method, c("LRT", "QLFT"))
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

    meta <- data.frame(sample_data(ps))
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

    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)

    # summarize data and  add norm.factors var to samples of DGEList
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
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
    # if (!length(confounders)) {
    #     model_data <- data.frame(group = groups)
    #     design <- stats::model.matrix(~ 0 + group, data = model_data)
    # } else {
    #     model_data <- data.frame(group = groups)
    #     model_data[confounders] <- meta[confounders]
    #     design <- stats::model.matrix(
    #         formula(paste(
    #             "~ + ", 
    #             paste(c(confounders, "group"), collapse = " + "))),
    #         data = model_data
    #     )
    # }
    design <- create_design(groups, meta, confounders)
    disp_para <- c(disp_para, y = list(dge_summarized), design = list(design))
    dge_summarized <- do.call(edgeR::estimateDisp, disp_para)

    # differential expression
    # quasi-likelihood (QL) F-test is used as it reflects the uncertainty in
    #  estimating the dispersion for each feature, and gives stricter error
    #  rate control
    fit_fun <- ifelse(method == "LRT", edgeR::glmFit, edgeR::glmQLFit)
    test_fun <- ifelse(method == "LRT", edgeR::glmLRT, edgeR::glmQLFTest)
    fit <- fit_fun(dge_summarized, design, ...)
    para_cf <- calc_coef(groups, design, contrast)
    lrt <- test_fun(fit, coef = para_cf)
    # lrt <- test_fun(fit, contrast = contrast_new)
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
    ef_nf <- dge_summarized$samples$lib.size * 
        dge_summarized$samples$norm.factors
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

    if (n_lvl > 2) {
        if (is.null(contrast)) {
            cf <- fit$coefficients
            target_idx <- grepl("group", colnames(cf))
            cf <- cf[, target_idx]
            # the first coef is intercept, bind the the reference group as 0
            # (the first column)
            cf <- cbind(0, cf)
            enrich_group <- lvl[apply(cf, 1, which.max)]
            # sort the enrich_group according to the DE of topTags
            de_idx <- match(row.names(res), row.names(cf))
            enrich_group <- enrich_group[de_idx]
        } else {
            enrich_group <- ifelse(res$logFC > 0, contrast[2], contrast[1])
        }
    } else {
        enrich_group <- ifelse(res$logFC > 0, lvl[2], lvl[1])
    }
    res$enrich_group <- enrich_group

    # edgeR::decideTestsDGE(), identify which genes are significantly
    # differentially expressed from an edgeR fit object containing p-values and
    # test statistics.
    
    # first two columns: feature enrich_group (write a function)
    res <- cbind(feature = row.names(res), res)
    other_col <- setdiff(names(res), c("feature", "enrich_group"))
    res <- res[, c("feature", "enrich_group", other_col)]
    row.names(res) <- paste0("marker", seq_len(nrow(res)))

    # var of effect size: named as ef_<name> of the actual effect size,
    # two groups: logFC, multiple groups: F for QLFT method, LR for LFT method
    if (n_lvl == 2) { # two groups
        ef_name <- "logFC"
    } else { # multiple groups
        if (!is.null(contrast)) { # two groups comparison from multiple groups
            ef_name <- "logFC"
        } else {
            ef_name <- ifelse(method == "LRT", "LR", "F")
        }
    }

    # only keep five variables: feature, enrich_group, effect_size (LR for LRT
    # F for QLFT), pvalue, and padj, write a function? select_marker_var
    # (effect_size = "")
    keep_var <- c("feature", "enrich_group", ef_name, "pvalue", "padj")
    res <- res[keep_var]
    names(res)[3] <- paste0("ef_", ef_name)
    sig_res <- res[res$padj < pvalue_cutoff & !is.na(res$padj), ]
    marker <- return_marker(sig_res, res)

    marker <- microbiomeMarker(
        marker_table = marker,
        norm_method = ifelse(is.null(nf), "TMM", get_norm_method(norm)),
        diff_method = paste("edgeR:", method),
        sam_data = sample_data(ps_normed),
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
#' @examples
#' data(caporaso)
#' dge <- phyloseq2edgeR(caporaso)
phyloseq2edgeR <- function(ps, ...) {
    ps <- keep_taxa_in_rows(ps)
    abd <- as(otu_table(ps), "matrix")
    
    if (any(round(abd) != abd)) {
        warning(
            "Some counts are non-integers, they are rounded to integers.\n",
            "Raw count is recommended for reliable results for edger method.",
            call. = FALSE
        )
    }

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
