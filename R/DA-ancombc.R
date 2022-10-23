#' Differential analysis of compositions of microbiomes with bias correction
#' (ANCOM-BC).
#'
#' Differential abundance analysis for microbial absolute abundance data. This
#' function is a wrapper of [`ANCOMBC::ancombc()`].
#'
#' @param ps  a [`phyloseq::phyloseq-class`] object, which consists of a feature
#'   table, a sample metadata and a taxonomy table.
#' @param group the name of the group variable in metadata. Specifying
#'   `group` is required for detecting structural zeros and performing
#'   global test.
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
#' @param norm_para  named `list`. other arguments passed to specific
#'   normalization methods.  Most users will not need to pass any additional
#'   arguments here.
#' @param p_adjust method to adjust p-values by. Default is "holm".
#'   Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none". See [`stats::p.adjust()`] for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#'   less than `prv_cut` will be excluded in the analysis. Default
#'   is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#'   sizes. Samples with library sizes less than `lib_cut` will be excluded
#'   in the analysis. Default is 0, i.e. do not filter any sample.
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
#' @param pvalue_cutoff level of significance. Default is 0.05.
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
#' @references
#' Lin, Huang, and Shyamal Das Peddada. "Analysis of compositions of microbiomes
#' with bias correction." Nature communications 11.1 (2020): 1-11.
#'
#' @seealso [`ANCOMBC::ancombc`]
#'
#' @importFrom ANCOMBC ancombc
#' @importFrom stats relevel
#' @export
#' @return a [`microbiomeMarker-class`] object.
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")
#' )
#' run_ancombc(ps, group = "Enterotype")
run_ancombc <- function(ps,
    group,
    confounders = character(0),
    contrast = NULL,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "none",
    norm_para = list(),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    prv_cut = 0.1,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = FALSE,
    tol = 1e-05,
    max_iter = 100,
    conserve = FALSE,
    pvalue_cutoff = 0.05) {
    stopifnot(inherits(ps, "phyloseq"))
    ps <- check_rank_names(ps) %>%
        check_taxa_rank( taxa_rank)

    if (length(confounders)) {
        confounders <- check_confounder(ps, group, confounders)
    }

    # if it contains missing values for any
    # variable specified in the formula, the corresponding sampling fraction
    # estimate for this sample will return NA since the sampling fraction is
    # not estimable with the presence of missing values.
    # remove this samples
    fml_char <- ifelse(length(confounders),
                       paste(c(confounders, group), collapse = " + "),
                       group)
    # fml_char <- paste(c(confounders, group), collapse = " + ")
    # fml <- stats::as.formula(paste("~", fml_char))
    # vars_fml <- all.vars(fml)
    for (var in c(confounders, group)) {
        ps <- remove_na_samples(ps, var)
    }

    # check whether group is valid, write a function
    meta <- sample_data(ps)
    meta_nms <- names(meta)
    groups <- meta[[group]]
    groups <- make.names(groups)
    if (!is.null(contrast)) {
        contrast <- make.names(contrast)
    }
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    groups <- set_lvl(groups, contrast)
    sample_data(ps)[[group]] <- groups
    lvl <- levels(groups)
    n_lvl <- length(lvl)

    contrast <- check_contrast(contrast)

    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )

    # set the reference level for pair-wise comparison from mutliple groups
    # if (!is.null(contrast) && n_lvl > 2) {
    #     groups <- relevel(groups, ref = contrast[1])
    #     sample_data(ps)[[group]] <- groups
    # }


    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # normalize the data
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)

    global <- ifelse(n_lvl > 2, TRUE, FALSE)
    # ancombc differential abundance analysis
    
    if (taxa_rank == "all") {
        ancombc_taxa_rank <- rank_names(ps_summarized)[1]
    } else {
        ancombc_taxa_rank <- taxa_rank
    }
    
    ancombc_out <- ANCOMBC::ancombc(
        ps_summarized,
        tax_level = ancombc_taxa_rank,
        formula = fml_char,
        p_adj_method = p_adjust,
        prv_cut = prv_cut,
        lib_cut = lib_cut,
        group = group,
        struc_zero = struc_zero,
        neg_lb = neg_lb,
        tol = tol,
        max_iter = max_iter,
        conserve = conserve,
        alpha = pvalue_cutoff,
        global = global
    )

    # multiple-group comparison will be performed while the group
    # variable has > 2 levels
    keep_var <- c("W", "p_val", "q_val", "diff_abn")
    if (n_lvl > 2) {
        # ANCOM-BC global test to determine taxa that are differentially
        # abundant between three or more groups of multiple samples.
        # global result to marker_table
        if (is.null(contrast)) {
            mtab <- ancombc_out$res_global
        } else {
            exp_lvl <- paste0(group, contrast[2])
            ancombc_res <- ancombc_out$res
            mtab <- lapply(keep_var, function(x) ancombc_res[[x]][exp_lvl])
            mtab <- do.call(cbind, mtab)
        }
    } else {
        ancombc_out_res <- ancombc_out$res
        # drop intercept
        ancombc_out_res <- lapply(
            ancombc_out_res,
            function(x) x[-1])
        mtab <- do.call(
            cbind,
            ancombc_out_res[c("W", "p_val", "q_val", "diff_abn")]
        )
    }
    names(mtab) <- keep_var

    # determine enrich group based on coefficients
    # drop intercept
    cf <- ancombc_out$res$lfc[-1]
    if (n_lvl > 2) {
        if (!is.null(contrast)) {
            cf <- cf[exp_lvl]
            enrich_group <- ifelse(cf[[1]] > 0, contrast[2], contrast[1])
        } else {
            cf <- cbind(0, cf)
            enrich_group <- lvl[apply(cf, 1, which.max)]
        }
    } else {
        enrich_group <- ifelse(cf[[1]] > 0, lvl[2], lvl[1])
    }

    # # enriched group
    enrich_abd <- get_ancombc_enrich_group(ps_summarized, ancombc_out, group)
    norm_abd <- enrich_abd$abd

    mtab <- cbind(feature = row.names(mtab), mtab, enrich_group)
    mtab_sig <- mtab[mtab$diff_abn, ]
    mtab_sig <- mtab_sig[c("feature", "enrich_group", "W", "p_val", "q_val")]
    names(mtab_sig) <- c("feature", "enrich_group", "ef_W", "pvalue", "padj")
    marker <- return_marker(mtab_sig, mtab)
    marker <- microbiomeMarker(
        marker_table = marker,
        norm_method = get_norm_method(norm),
        diff_method = "ancombc",
        sam_data = sample_data(ps_summarized),
        otu_table = otu_table(norm_abd, taxa_are_rows = TRUE),
        tax_table = tax_table(ps_summarized)
    )

    marker
}

get_ancombc_enrich_group <- function(ps, ancombc_out, group) {
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

    groups <- sample_data(ps)[[group]]
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
