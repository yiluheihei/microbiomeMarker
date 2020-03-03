#' Liner discriminant analysis (LDA) effect size (LEFSe) analysis
#'
#' Perform Metagenomic LEFSe analysis based on phyloseq object.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object
#' @param tax_rank character, rank level
#' @param class character, the column name to specify the class
#' @param p_cutoff numeric, p value cutoff, default 0.05
#' @param lda_cutoff numeric, lda score cutoff, default 2
#' @importFrom  dplyr mutate filter arrange rowwise select
#' @importFrom  purrr map_dbl pmap_dbl pmap_chr
#' @importFrom stats p.adjust
#' @export
#' @author Yang Cao \email{yiluheihei@gmail.com}
#' @references Segata, Nicola, et al. Metagenomic biomarker discovery and
#' explanation. Genome biology 12.6 (2011): R60.
lefse <- function(ps, tax_rank, class, p_cutoff = 0.05, lda_cutoff = 2) {
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be phyloseq object")
  }

  sample_meta <- sample_data(ps)
  class <- sample_meta[[class]]

  otus <- otu_table(ps)
  if (taxa_are_rows(otus)) {
    otus <- t(otus)
  }
  otus <- tibble::as_tibble(otus@.Data, rownames = NA)

  # kw rank sum test among classes
  kw_p <- map_dbl(otus, ~ kruskal.test(.x, class)$p.value)
  kw_fdr <- p.adjust(kw_p, method = "fdr")

  # wilcoxon rank sum test is not preformed if there is no subclass

  # lda analysis
  lda_res <- MASS::lda(class~., data = otus)

  # dplyr verbs drop the row names automatically, otu_id is saved and added after
  # dplyr manipulation complete
  otu_id <- colnames(lda_res$means)
  lda_mean <- t(lda_res$means) %>%
    as_tibble()

  lda_max <- pmap_dbl(lda_mean, max)
  lda_min <- pmap_dbl(lda_mean, min)
  enriched <- pmap_chr(lda_mean, enrich_group)

  lda_mean <- mutate(lda_mean, max = lda_max, min = lda_min) %>%
    mutate(lda_score = signif(log10(1 + abs(max - min)/2), digits = 3)) %>%
    mutate(p_value = kw_p, fdr = kw_fdr, enrich_group = enriched) %>%
    mutate(otu = gsub("`", "", otu_id))# by default, as_tibble add ` around invalid names

  # significant feature,  order result by lda_score
  sig_feature <- filter(lda_mean, fdr < p_cutoff, lda_score > lda_cutoff) %>%
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
