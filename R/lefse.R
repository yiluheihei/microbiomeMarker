#' Liner discriminant analysis (LDA) effect size (LEFSe) analysis
#'
#' Perform Metagenomic LEFSe analysis based on phyloseq object.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object
#' @param class character, the column name to set the class
#' @param subclass character, the column name to set the subclass
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
#' @importFrom  dplyr mutate filter arrange rowwise select
#' @importFrom  purrr map_dbl pmap_dbl pmap_chr
#' @importFrom stats p.adjust
#' @importFrom phyloseq rank_names
#' @export
#' @return a [microbiomeMarker-class] object, in which the `slot` of `marker_table`
#' contains five variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched
#' * `log_max_mean`, the logarithm value of the highest mean among all the
#'   classes
#' * `lda`, logarithmic LDA score
#' * `p_value`, p value of kw test.
#' @author Yang Cao \email{yiluheihei@gmail.com}
#' @seealso [normalize]
#' @references Segata, Nicola, et al. Metagenomic biomarker discovery and
#' explanation. Genome biology 12.6 (2011): R60.
lefse <- function(ps,
                  class,
                  subclass = NULL,
                  transform = c("identity", "log10", "log10p"),
                  norm = 1000000,
                  norm_para = list(),
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

  if (!check_rank_names(ps)) {
    stop(
      "ranks of `ps` must be one of ",
      paste(available_ranks, collapse = ", ")
    )
  }

  transform <- match.arg(transform, c("identity", "log10", "log10p"))
  correct <- match.arg(correct, c("0", "1", "2"))
  correct <- as.numeric(correct)

  # import input from the original lefse python script or galaxy,
  # will be dropped in the next release version
  summarized <- check_tax_summarize(ps)
  if (summarized && !is.numeric(norm)) {
    stop(
      '`norm` must be a numeric or "none" while `ps` has been summarized',
      call. = FALSE
    )
  }

  # pre-processing, including: keep taxa in rows, filter taxa whose abundance is
  # zero, fix duplicated tax, transformation and normalization
  ps <- preprocess_ps(ps)
  # transformation
  ps <- transform_abundances(ps, transform = transform)
  # normalization
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)
  # ps <- normalize(ps, norm, ...)

  sample_meta <- sample_data(ps_normed)
  cls_info <- lefse_format_class(sample_meta, class, subcls = subclass)
  cls <- cls_info$cls
  subcls <- cls_info$subcls
  cls_hie <- cls_info$cls_hie

  ps_summarized <- summarize_taxa(ps_normed)
  # otus <- otu_table(ps_summarized)
  otus <- abundances(ps_summarized, norm = TRUE)
  # otus_norm <- normalize_feature(otus, normalization = normalization)
  # transform it for test
  otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
  feature <- tax_table(ps_summarized)@.Data[, 1]
  names(otus_test) <- feature

  # kw rank sum test among classes
  kw_p <- purrr::map_dbl(otus_test, ~ kruskal.test(.x, cls)$p.value)

  # remove the taxa, while pvalue is na
  na_ind <- is.na(kw_p)
  if (sum(na_ind) >= 1) {
    otus_test <- otus_test[!na_ind]
    kw_p <- kw_p[!na_ind]
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
  row.names(lefse_out) <- paste0("marker", seq_len(nrow(lefse_out)))

  # if (summarize == "lefse" || summarize) {
  #   tax <- matrix(row.names(otus)) %>%
  #     tax_table()
  # } else {
  #   tax <- tax_table(ps)
  # }
  # row.names(tax) <- row.names(otus)


  tax <- matrix(feature) %>%
      tax_table()
  row.names(tax) <- row.names(otus)

  mm <- microbiomeMarker(
    marker_table = lefse_out,
    tax_table_orig = tax_table(ps),
    otu_table(otus, taxa_are_rows = TRUE), # normalized feature table
    tax
  )

  mm
}

# suppress the checking notes “no visible binding for global variable", which is
# caused by NSE
# use rlang:.data
# utils::globalVariables(c("lda"))
