#' Liner discriminant analysis (LDA) effect size (LEFSe) analysis
#'
#' Perform Metagenomic LEFSe analysis based on phyloseq object.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object
#' @param group character, the column name to set the group
#' @param subgroup character, the column name to set the subgroup
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
#' @param kw_cutoff numeric, p value cutoff of kw test, default 0.05
#' @param wilcoxon_cutoff numeric, p value cutoff of wilcoxon test, default 0.05
#' @param lda_cutoff numeric, lda score cutoff, default 2
#' @param bootstrap_n integer, the number of bootstrap iteration for LDA,
#'   default 30
#' @param bootstrap_fraction numeric, the subsampling fraction value for each
#'   bootstrap iteration, default `2/3`
#' @param multigrp_strat logical, for multiple group tasks, whether the test is
#'   performed in a one-against one (more strict) or in a one-against all
#'   setting, default `FALSE`.
#' @param strict multiple testing options, 0 for no correction (default), 1 for
#'   independent comparisons, 2 for independent comparison.
#' @param sample_min integer, minimum number of samples per subclass for
#'   performing wilcoxon test, default 10
#' @param only_same_subgrp logical, whether perform the wilcoxon test only
#'   among the subgroups with the same name, default `FALSE`
#' @param curv logical, whether perform the wilcoxon test using the
#'   Curtis's approach, defalt `FALSE`
#' @importFrom  dplyr mutate filter arrange rowwise select
#' @importFrom  purrr map_dbl pmap_dbl pmap_chr
#' @importFrom stats p.adjust
#' @importFrom phyloseq rank_names tax_glom
#' @export
#' @return a [microbiomeMarker-class] object, in which the `slot` of
#' `marker_table`
#' contains four variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched.
#' * `lda`, logarithmic LDA score (effect size)
#' * `pvalue`, p value of kw test.
#' @author Yang Cao
#' @seealso [normalize]
#' @references Segata, Nicola, et al. Metagenomic biomarker discovery and
#' explanation. Genome biology 12.6 (2011): R60.
#' @examples
#' data(kostic_crc)
#' kostic_crc_small <- phyloseq::subset_taxa(
#'     kostic_crc,
#'     Phylum == "Firmicutes"
#' )
#' mm_lefse <- run_lefse(
#'     kostic_crc_small,
#'     wilcoxon_cutoff = 0.01,
#'     group = "DIAGNOSIS",
#'     kw_cutoff = 0.01,
#'     multigrp_strat = TRUE,
#'     lda_cutoff = 4
#' )
run_lefse <- function(ps,
    group,
    subgroup = NULL,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "CPM",
    norm_para = list(),
    kw_cutoff = 0.05,
    lda_cutoff = 2,
    bootstrap_n = 30,
    bootstrap_fraction = 2 / 3,
    wilcoxon_cutoff = 0.05,
    multigrp_strat = FALSE,
    strict = c("0", "1", "2"),
    sample_min = 10,
    only_same_subgrp = FALSE,
    curv = FALSE) {
    if (!inherits(ps, "phyloseq")) {
        stop("`ps` must be phyloseq object", call. = FALSE)
    }
    
    # check rank names and para taxa_rank
    ps <- check_rank_names(ps)
    ps <- check_taxa_rank(ps, taxa_rank)
    
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    strict <- match.arg(strict, c("0", "1", "2"))
    strict <- as.numeric(strict)

    # import input from the original lefse python script or galaxy,
    # will be dropped in the next release version
    summarized <- check_tax_summarize(ps)
    if (summarized && norm != "CPM") {
        stop(
            "`norm` must be a 'CPM' or 'none' while `ps` has been summarized",
            call. = FALSE
        )
    }

    # pre-processing, including: keep taxa in rows, filter taxa whose abundance
    # is zero, fix duplicated tax, transformation and normalization
    ps <- preprocess_ps(ps)
    # transformation
    ps <- transform_abundances(ps, transform = transform)
    # normalization
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)

    sample_meta <- sample_data(ps_normed)
    grp_info <- lefse_format_grp(sample_meta, group, subgroup = subgroup)
    grp <- grp_info$group
    subgrps <- grp_info$subgroup
    grp_hie <- grp_info$group_hie

    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    
    otus <- abundances(ps_summarized, norm = TRUE)
    # transform it for test
    otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
    feature <- tax_table(ps_summarized)@.Data[, 1]
    names(otus_test) <- feature
    
    # tax table
    tax <- matrix(feature) %>%
        tax_table()
    row.names(tax) <- row.names(otus)

    # kw rank sum test among classes
    kw_p <- purrr::map_dbl(otus_test, ~ kruskal.test(.x, grp)$p.value)

    # remove the taxa, while pvalue is na
    na_ind <- is.na(kw_p)
    if (sum(na_ind) >= 1) {
        otus_test <- otus_test[!na_ind]
        kw_p <- kw_p[!na_ind]
    }

    sig_ind <- kw_p <= kw_cutoff
    sig_otus <- otus_test[, sig_ind]

    # wilcox test is preformed for each class, if there is no subclass
    features_nms <- names(sig_otus)
    wilcoxon_p <- purrr::map2_lgl(
        sig_otus, features_nms,
        ~ test_rep_wilcoxon(
            subgrps, grp_hie,
            .x, .y,
            wilcoxon_cutoff = wilcoxon_cutoff,
            multicls_strat = multigrp_strat,
            strict = strict,
            sample_min = sample_min,
            only_same_subcls = only_same_subgrp,
            curv = curv
        )
    )
    sig_otus <- sig_otus[, wilcoxon_p, drop = FALSE]
    
    if (ncol(sig_otus) == 0) {
        warning("No marker was identified", call. = FALSE)
        mm <- microbiomeMarker(
            marker_table = NULL,
            norm_method = get_norm_method(norm),
            diff_method = "lefse",
            otu_table = otu_table(otus, taxa_are_rows = TRUE), # normalized
            # new var norm_factor (if it is calculated in normalize)
            sam_data = sample_data(ps_normed),
            tax_table = tax
        )
        
        return(mm)
    }

    # mean abundance in each group
    otus_enriched_group <- get_feature_enrich_group(grp, sig_otus)

    # bootsrap iteration of lda
    ldas <- bootstap_lda(
        sig_otus,
        boot_n = bootstrap_n,
        class = grp,
        sample_fract = bootstrap_fraction
    )

    lefse_res <- data.frame(
        feature = names(sig_otus),
        enrich_group = otus_enriched_group$group,
        ef_lda = ldas,
        pvalue = kw_p[sig_ind][wilcoxon_p],
        stringsAsFactors = FALSE
    )

    lefse_sig <- filter(lefse_res, .data$ef_lda >= lda_cutoff) %>%
        arrange(.data$enrich_group, desc(.data$ef_lda))

    lefse_out <- return_marker(lefse_sig, lefse_res)
    lefse_out$padj <- lefse_out$pvalue

    mm <- microbiomeMarker(
        marker_table = lefse_out,
        norm_method = get_norm_method(norm),
        diff_method = "lefse",
        otu_table = otu_table(otus, taxa_are_rows = TRUE), # normalized
        # new var norm_factor (if it is calculated in normalize)
        sam_data = sample_data(ps_normed),
        tax_table = tax
    )

    mm
}
