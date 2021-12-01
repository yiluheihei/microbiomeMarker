#' Simple statistical analysis of metagenomic profiles
#'
#' Perform simple statistical analysis of metagenomic profiles. This function
#' is a wrapper of `run_test_two_groups` and `run_test_multiple_groups`.
#'
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param group character, the variable to set the group
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
#' @param norm_para arguments passed to specific normalization methods
#' @param method test method, options include: "welch.test", "t.test" and
#'   "white.test" for two groups comparison,  "anova"and "kruskal" for multiple
#'   groups comparison.
#' @param p_adjust method for multiple test correction, default `none`,
#'   for more details see [stats::p.adjust].
#' @param pvalue_cutoff numeric, p value cutoff, default 0.05
#' @param diff_mean_cutoff,ratio_cutoff only used for two groups comparison,
#'   cutoff of different means and ratios, default `NULL` which means no effect
#'   size filter.
#' @param eta_squared_cutoff only used for multiple groups comparison, numeric,
#'   cutoff of effect size (eta squared) default `NULL` which means no effect
#'   size filter.
#' @param conf_level only used for two groups comparison, numeric, confidence
#'   level of interval.
#' @param nperm integer, only used for two groups comparison, number of
#'   permutations for white non parametric t test estimation
#' @param ... only used for two groups comparison, extra arguments passed to
#' [`t.test()`] or [`fisher.test()`].
#' @return a [`microbiomeMarker-class`] object.
#' @seealso [`run_test_two_groups()`],[`run_test_multiple_groups()`]
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")
#' )
#' run_simple_stat(ps, group = "Enterotype")
run_simple_stat <- function(ps,
    group,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "TSS",
    norm_para = list(),
    method = c(
        "welch.test", "t.test", "white.test",
        "anova", "kruskal"
    ),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    diff_mean_cutoff = NULL,
    ratio_cutoff = NULL,
    eta_squared_cutoff = NULL,
    conf_level = 0.95,
    nperm = 1000,
    ...) {
    stopifnot(inherits(ps, "phyloseq"))
    
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    method <- match.arg(
        method,
        c("welch.test", "t.test", "white.test", "anova", "kruskal")
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

    if (n_group == 2) {
        if (!method %in% c("welch.test", "t.test", "white.test")) {
            stop(
                "There are two groups here, please select welch.test, t.test, ",
                "or white.test for two groups comparison",
                call. = FALSE
            )
        }
        if (!missing(eta_squared_cutoff)) {
            warning(
                "`eta_squared_cutoff` is ignored since it is only used for ",
                "multiple groups comparison",
                call. = FALSE
            )
        }

        res <- run_test_two_groups(
            ps = ps,
            group = group,
            taxa_rank = taxa_rank,
            transform = transform,
            norm = norm,
            norm_para = norm_para,
            method = method,
            p_adjust = p_adjust,
            pvalue_cutoff = pvalue_cutoff,
            diff_mean_cutoff = diff_mean_cutoff,
            ratio_cutoff = ratio_cutoff,
            conf_level = conf_level,
            nperm = nperm,
            ...
        )
    } else {
        if (!method %in% c("anova", "kruskal")) {
            stop(
                "There are more than two groups, please select anova or ",
                "kruskal for multiple groups comparison",
                call. = FALSE
            )
        }

        if (!missing(diff_mean_cutoff)) {
            warning(
                "`diff_mean_cutoff` only worked for two groups comparison",
                call. = FALSE
            )
        }
        if (!missing(ratio_cutoff)) {
            warning(
                "`ratio_cutoff` only worked for two groups comparison",
                call. = FALSE
            )
        }
        if (!missing(nperm)) {
            warning(
                "`nperm` only worked for two groups comparison",
                call. = FALSE
            )
        }
        if (!missing(conf_level)) {
            warning(
                "`conf_level` only worked for two groups comparison",
                call. = FALSE
            )
        }

        res <- run_test_multiple_groups(
            ps = ps,
            group = group,
            taxa_rank = taxa_rank,
            transform = transform,
            norm = norm,
            norm_para = norm_para,
            method = method,
            p_adjust = p_adjust,
            pvalue_cutoff = pvalue_cutoff,
            effect_size_cutoff = eta_squared_cutoff
        )
    }

    res
}
