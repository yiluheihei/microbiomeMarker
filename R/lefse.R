#' Liner discriminant analysis (LDA) effect size (LEFSe) analysis
#'
#' Perform Metagenomic LEFSe analysis based on phyloseq object.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object
#' @param class character, the column name to set the class
#' @param subclass character, the varibale name to set the subclass
#' @param kw_cutoff numeric, p value cutoff of kw test, default 0.05
#' @param wilcoxon_cutoff numeric, p value cutoff of wilcoxon test, default 0.05
#' @param normalization norm set the normalization value
#' @param summarize logical, whether summarize the
#'   taxa (`TRUE` or `FALSE`), default `TRUE`
#' @param lda_cutoff numeric, lda score cutoff, default 2
#' @param bootstrap_n integer, the number of bootstrap iteration for LDA,
#'   default 30
#' @param bootstrap_fraction numeric, the subsampling fraction value for each
#'   bootstrap iteration, default `2/3`
#' @param multicls_strat logical, for multiple class tasks, whether the test is
#'   performed in a one-against one (more strict) or in a one-against all
#'   setting, default `FALSE`.
#' @param correct multiple testing options, 0 for no correction (default), 1 for
#'   independent comparisons, 2 for independent comparison
#' @param sample_min integer, minimum number of samples per subclass for
#'   performing wilcoxon test, default 10
#' @param only_same_subcls logical, whether perform the wilcoxon test only
#'   among the subclasses with the same name, default `FALSE`
#' @param curv logical, whether perform the wilcoxon test using the
#'   Curtis's approach, defalt `FALSE`
#'
#' @importFrom  dplyr mutate filter arrange rowwise select
#' @importFrom  purrr map_dbl pmap_dbl pmap_chr
#' @importFrom stats p.adjust
#' @importFrom phyloseq rank_names
#' @export
#' @return a data frame contains five variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched
#' * `log_max_mean`, the logarithm value of the highest mean among all the
#'   classes
#' * `lda`, logarithmic LDA score
#' * `p_value`, p value of kw test
#' @author Yang Cao \email{yiluheihei@gmail.com}
#' @references Segata, Nicola, et al. Metagenomic biomarker discovery and
#' explanation. Genome biology 12.6 (2011): R60.
lefse <- function(ps,
                  class,
                  subclass = NULL,
                  normalization = 1000000,
                  summarize = TRUE,
                  kw_cutoff = 0.05,
                  lda_cutoff = 2,
                  bootstrap_n = 30,
                  bootstrap_fraction = 2/3,
                  wilcoxon_cutoff = 0.05,
                  multicls_strat = FALSE,
                  correct = c("0", "1", "2"),
                  sample_min = 10,
                  only_same_subcls = FALSE,
                  curv = FALSE) {
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be phyloseq object", call. = FALSE)
  }

  correct <- match.arg(correct)
  correct <- as.numeric(correct)

  ranks <- rank_names(ps)
  diff_rank <- setdiff(ranks, availabel_ranks)
  if (length(diff_rank)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(availabel_ranks, collapse = ", ")
    )
  }

  # filter the taxa whose abundance is zero
  ps <- phyloseq_qc(ps)
  # fix duplicated tax
  ps <- fix_duplicate_tax(ps)

  sample_meta <- sample_data(ps)
  cls_info <- lefse_format_class(sample_meta, class, subcls = subclass)
  cls <- cls_info$cls
  subcls <- cls_info$subcls
  cls_hie <- cls_info$cls_hie

  # check whether the taxa is summarized first
  # e.g. phyloseq object is construct from lefse python script
  # just for experiment, may be dropped in the future
  if (check_tax_summarize(ps)) {
    otus <- otu_table(ps) %>%
      add_missing_levels()
  } else if (!summarize) {
    otus <- otu_table(ps)
  } else {
    otus <- summarize_taxa(ps)
  }
  otus_norm <- normalize_feature(otus, normalization = normalization)
  # transform it for test
  otus_test <- as.data.frame(t(otus_norm), stringsAsFactors = FALSE)

  # kw rank sum test among classes
  kw_p <- purrr::map_dbl(otus_test, ~ kruskal.test(.x, cls)$p.value)

  # remove the taxa, while pvalue is na
  na_ind <- is.na(kw_p)
  if (sum(na_ind) >= 1) {
    otus_test <- otus_test[-na_ind]
    kw_p <- kw_p[-na_ind]
  }

  sig_ind <- kw_p <= kw_cutoff
  sig_otus <- otus_test[, sig_ind]

  # wilcoxon rank sum test is preformed for each class, if there is no subclass
  features_nms <- names(sig_otus)
  wilcoxon_p <- purrr::map2_lgl(
    sig_otus, features_nms,
    ~ test_rep_wilcoxon(
      subcls, cls_hie,
      .x, .y,
      wilcoxon_cutoff = wilcoxon_cutoff,
      multicls_strat = multicls_strat,
      correct = correct,
      sample_min = sample_min,
      only_same_subcls = only_same_subcls,
      curv = curv
    )
  )
  sig_otus <- sig_otus[, wilcoxon_p]

  # mean abundance in each group
  otus_enriched_group <- get_feature_enrich_group(cls, sig_otus)

  # bootsrap iteration of lda
  ldas <- bootstap_lda(
    sig_otus,
    boot_n = bootstrap_n,
    class = cls,
    sample_fract = bootstrap_fraction
  )

  lefse_out <- data.frame(
    feature = names(sig_otus),
    enrich_group = otus_enriched_group$group,
    log_max_mean = otus_enriched_group$log_max_mean,
    lda = ldas,
    p_value = kw_p[sig_ind][wilcoxon_p],
    stringsAsFactors = FALSE) %>%
    filter(.data$lda >= lda_cutoff) %>%
    arrange(.data$enrich_group, desc(.data$lda)) %>%
    marker_table()

  if (summarize == "lefse" || summarize) {
    tax <- matrix(row.names(otus)) %>%
      tax_table()
    row.names(tax) <- row.names(otus)
  } else {
    tax <- tax_table(ps)
  }

  mm <- microbiomeMarker(lefse_out, otus_norm, tax)

  mm
}

# suppress the checking notes “no visible binding for global variable", which is
# caused by NSE
# use rlang:.data
# utils::globalVariables(c("lda"))
