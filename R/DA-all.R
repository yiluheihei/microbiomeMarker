#' Find makers (differentially expressed metagenomic features)
#'
#' `run_marker` is a wrapper of all differential analysis functions.
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
#' @param da_method character to specify the differential analysis method. The
#'   options include:
#'   * "lefse", linear discriminant analysis (LDA) effect size (LEfSe) method,
#'     for more details see [`run_lefse()`].
#'   * "simple_t", "simple_welch", "simple_white", "simple_kruskal",
#'     and "simple_anova", simple statistic methods; "simple_t", "simple_welch"
#'     and "simple_white" for two groups comparison; "simple_kruskal", and
#'     "simple_anova" for multiple groups comparison. For more details see
#'     [`run_simple_stat()`].
#'   * "edger", see [`run_edger()`].
#'   * "deseq2", see [`run_deseq2()`].
#'   * "metagenomeseq", differential expression analysis based on the
#'     Zero-inflated Log-Normal mixture model or Zero-inflated Gaussian mixture
#'     model using metagenomeSeq, see [`run_metagenomeseq()`].
#'   * "ancom", see [`run_ancom()`].
#'   * "ancombc", differential analysis of compositions of microbiomes with
#'     bias correction, see [`run_ancombc()`].
#'   * "aldex", see [`run_aldex()`].
#'   * "limma_voom", see [`run_limma_voom()`].
#'   * "sl_lr", "sl_rf", and "sl_svm", there supervised leaning (SL) methods:
#'     logistic regression (lr), random forest (rf), or support vector machine
#'     (svm). For more details see [`run_sl()`].
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
#' @param norm_para arguments passed to specific normalization methods
#' @param p_adjust method for multiple test correction, default `none`,
#'   for more details see [stats::p.adjust].
#' @param pvalue_cutoff numeric, p value cutoff, default 0.05.
#' @param ... extra arguments passed to the corresponding differential analysis
#'   functions, e.g. [`run_lefse()`].
#' @return a [`microbiomeMarker-class`] object.
#' @details This function is only a wrapper of all differential analysis
#'   functions, We recommend to use the corresponding function, since it has a
#'   better default arguments setting.
#' @export
#' @seealso [`run_lefse()`],[`run_simple_stat()`],[`run_test_two_groups()`],
#'   [`run_test_multiple_groups()`],[`run_edger()`],[`run_deseq2()`],
#'   [`run_metagenomeseq`],[`run_ancom()`],[`run_ancombc()`],[`run_aldex()`],
#'   [`run_limma_voom()`],[`run_sl()`]
run_marker <- function(ps,
    group,
    da_method = c(
        "lefse", "simple_t", "simple_welch",
        "simple_white", "simple_kruskal",
        "simple_anova", "edger", "deseq2",
        "metagenomeseq", "ancom", "ancombc", "aldex",
        "limma_voom", "sl_lr", "sl_rf", "sl_svm"
    ),
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "none",
    norm_para = list(),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    ...) {
    stopifnot(inherits(ps, "phyloseq"))

    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    da_method <- match.arg(
        da_method,
        c(
            "lefse", "simple_t", "simple_welch",
            "simple_white", "simple_kruskal", "simple_anova",
            "edger", "deseq2", "metagenomeseq", "ancom",
            "ancombc", "aldex", "limma_voom",
            "sl_lr", "sl_rf", "sl_svm"
        )
    )
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )

    # group
    sample_meta <- sample_data(ps)
    if (!group %in% names(sample_meta)) {
        stop("`group` must in the field of sample meta data", call. = FALSE)
    }
    groups <- sample_meta[[group]]
    n_group <- length(unique(groups))
    if (n_group == 1) {
        stop("at least two groups required", call. = FALSE)
    }

    para <- c(
        list(...),
        ps = ps,
        group = group,
        taxa_rank = taxa_rank,
        transform = transform,
        norm = norm,
        norm_para = norm_para,
        p_adjust = p_adjust,
        pvalue_cutoff = pvalue_cutoff
    )

    test_fun <- switch(da_method,
        lefse = run_lefse,
        edger = run_edger,
        metagenomeseq = run_metagenomeseq,
        deseq2 = run_deseq2,
        ancom = run_ancom,
        ancombc = run_ancombc,
        aldex = run_aldex,
        limma_voom = run_limma_voom
    )

    if (da_method == "lefse") {
        para <- c(
            list(...),
            ps = ps,
            class = group,
            taxa_rank = taxa_rank,
            transform = transform,
            norm = norm,
            norm_para = norm_para
        )
    }
    if (da_method %in% c(
        "simple_t", "simple_welch",
        "simple_white", "simple_kruskal", "simple_anova"
    )) {
        test_method <- switch(da_method,
            simple_t = "t.test",
            simple_wilch = "welch.test",
            simple_white = "white.test",
            simple_anova = "anova",
            simple_kruskal = "kruskal"
        )
        para <- c(para, method = test_method)
        test_fun <- run_simple_stat
    }
    if (da_method %in% c("sl_lr", "sl_rf", "sl_svm")) { # sl methods
        sl_method <- gsub("sl_", "", da_method) %>%
            toupper()
        para <- c(
            list(...),
            ps = ps,
            group = group,
            taxa_rank = taxa_rank,
            transform = transform,
            norm = norm,
            norm_para = norm_para,
            method = sl_method
        )
        test_fun <- run_sl
    }
    mm <- do.call(test_fun, para)

    mm
}
