#' Summarize taxa into a taxonomic level within each sample
#'
#' Provides summary information of the representation of a taxonomic levels within
#' each sample.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param level taxonomic level to summarize, default the top level rank of the
#'  `ps`.
#' @param absolute logical, whether return the absolute abundance or
#'   relative abundance, default `FALSE`
#' FALSE.
#' @param sep a character string to separate the taxonomic levels.
#' @return a [`phyloseq::otu_table-class`] object, where each row represents a
#'   taxa, and each col represents the taxa abundance of each sample.
#' @export
summarize_taxa <- function(ps,
                           level = rank_names(ps)[1],
                           absolute = TRUE,
                           sep = "|") {
  ps_ranks <- rank_names(ps)
  if (!level %in% ps_ranks) {
    stop("`level` must in the ranks of `ps` (rank_names(ps))")
  }

  # ranks <- setdiff(availabel_ranks, "Summarize")
  # level <- match(level, ranks)

  ind <- match(level, ps_ranks)
  levels <- ps_ranks[ind:length(ps_ranks)]
  res <- purrr::map(
    levels,
    ~.summarize_taxa_level(
      ps,
      rank = .x,
      absolute = absolute,
      sep = sep
    )
  )
  tax_nms <- purrr::map(res, row.names) %>% unlist()
  res <- bind_rows(res)
  row.names(res) <- tax_nms

  otu_table(res, taxa_are_rows = TRUE)
}

#' Summarize the taxa for the specific rank
#' @noRd
.summarize_taxa_level <- function(ps,
                                  rank_name,
                                  absolute = TRUE,
                                  sep = "|") {
  if (!absolute) {
    ps <- transform_sample_counts(ps, function(x)x/sum(x))
  }

  # norm the abundance data
  # if (norm > 0) {
  #   ps@otu_table <- ps@otu_table*norm
  # }

  otus <- otu_table(ps)
  otus_extend <- slot(otus, ".Data") %>%
    tibble::as_tibble()

  taxas <- tax_table(ps)@.Data %>%
    tibble::as_tibble()

  ranks <- setdiff(availabel_ranks, "Summarize")
  rank_level <- match(rank_name, ranks)
  select_ranks <- intersect(ranks[1:rank_level], rank_names(ps))

  consensus <- taxas[, select_ranks]  %>%
    purrr::pmap_chr(paste, sep = sep)
  otus_extend$consensus <- consensus

  taxa_summarized <- group_split(otus_extend, consensus) %>%
    purrr::map(.sum_consensus)
  taxa_summarized <- do.call(rbind, taxa_summarized)
  # filter taxa of which abundance is zero
  ind <- rowSums(taxa_summarized) != 0
  taxa_summarized <- taxa_summarized[ind, ]
    # tibble::rownames_to_column(var = "taxa") %>%
    # arrange(taxa)

  taxa_summarized
}

#' sum all otus which belongs to the same taxa
#' @noRd
.sum_consensus <- function(x) {
  consensus <- unique(x$consensus)
  if (length(consensus) != 1) {
    stop("consensus in the same group muste be the same")
  }

  x$consensus <- NULL
  res <- as.data.frame(t(colSums(x)))
  row.names(res) <- consensus

  return(res)
}

# suppress the checking notes “no visible binding for global variable", which is
# caused by NSE
# utils::globalVariables(c(".", "taxa", "as_tibble"))
