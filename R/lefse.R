#' Liner discriminant analysis (LDA) effect size (LEFSe) analysis
#'
#' Perform Metagenomic LEFSe analysis based on phyloseq object.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object
#' @param class character, the column name to specify the class
#' @param p_cutoff numeric, p value cutoff, default 0.05
#' @param norm norm set the normalization value
#' @param by_otu by otu or summarized taxa, default false
#' @param lda_cutoff numeric, lda score cutoff, default 2
#' @param correct use corrected p value (fdr) or not, default false
#'
#' @importFrom  dplyr mutate filter arrange rowwise select
#' @importFrom  purrr map_dbl pmap_dbl pmap_chr
#' @importFrom stats p.adjust
#' @export
#' @author Yang Cao \email{yiluheihei@gmail.com}
#' @references Segata, Nicola, et al. Metagenomic biomarker discovery and
#' explanation. Genome biology 12.6 (2011): R60.
lefse <- function(ps,
                  class,
                  norm = 1000000,
                  by_otu = FALSE,
                  p_cutoff = 0.05,
                  lda_cutoff = 2,
                  corrected = FALSE) {
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be phyloseq object")
  }

  # filter the taxa whose abundance is zero
  ps <- phyloseq_qc(ps)

  sample_meta <- sample_data(ps)
  class <- sample_meta[[class]]

  if (by_otu) {
    otus <- otu_table(ps)
    if (taxa_are_rows(otus)) {
      otus <- t(otus)
    }
    otus <- tibble::as_tibble(otus@.Data, rownames = NA)
  } else {
    otus <- summarize_taxa(ps, norm = norm) %>%
      t()
    taxa <- otus[1, ]
    otus <- tibble::as_tibble(otus, .name_repair = "unique")
    names(otus) <- taxa
    otus <- dplyr::slice(otus, -1) %>%
      purrr::map_df(as.numeric)
  }

  # kw rank sum test among classes
  kw_p <- map_dbl(otus, ~ kruskal.test(.x, class)$p.value)
  kw_fdr <- p.adjust(kw_p, method = "fdr")

  if (corrected) {
    sig_ind <- kw_fdr <= p_cutoff
  } else {
    sig_ind <- kw_p <= p_cutoff
  }
  otus <- otus[, sig_ind]


  # wilcoxon rank sum test is not preformed if there is no subclass

  # lda analysis
  lda_res <- MASS::lda(
    class ~ .,
    data = otus,
    tol = 1.0e-10
  )

  # dplyr verbs drop the row names automatically, otu_id is saved and added after
  # dplyr manipulation complete
  otu_id <- colnames(lda_res$means)
  lda_mean <- t(lda_res$means) %>%
    tibble::as_tibble()

  lda_max <- pmap_dbl(lda_mean, max)
  lda_min <- pmap_dbl(lda_mean, min)
  enriched <- pmap_chr(lda_mean, enrich_group)

  lda_mean <- mutate(lda_mean, max = lda_max, min = lda_min) %>%
    mutate(lda_score = signif(log10(1 + abs(max - min)/2), digits = 5)) %>%
    mutate(p_value = kw_p[sig_ind], fdr = kw_fdr[sig_ind], enrich_group = enriched) %>%
    mutate(otu = gsub("`", "", otu_id))# by default, as_tibble add ` around invalid names

  # significant feature,  order result by lda_score
  sig_feature <- filter(lda_mean, lda_score >= lda_cutoff) %>%
    arrange(desc(lda_score))

  sig_feature
}

# suppress the checking notes “no visible binding for global variable", which is
# caused by NSE
utils::globalVariables(c("fdr", "lda_score", "desc"))

#' helper function to extract enrich group
#' @noRd
enrich_group <- function(...) {
  lda_mean <- c(...)
  group_name <- names(lda_mean)

  group_name[which.max(lda_mean)]
}
