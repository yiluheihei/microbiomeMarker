#' Summarize taxa into a taxonomic level within each sample
#'
#' Provids summary information of the representation of a taxonomic levels within
#' each sample.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param level integer, taxonomic level to summarize by, default 7.
#' @param norm set the normalization value
#' @param absolute logical, whether return the absolute abundance or not, default
#' FALSE.
#' @param sep a character string to separate the taxonomic levels.
#'
#' @return a data frame, each row represnets a taxa, where each col represents
#' the taxa abunance of each sample
#' @export

summarize_taxa <- function(ps,
                           level = 7,
                           norm = 1000000,
                           absolute = FALSE,
                           sep = "|") {
  res <- purrr::map_dfr(
    1:level,
    ~.summarize_taxa_level(
      ps,
      rank = .x,
      norm = norm,
      absolute = absolute,
      sep = sep
    )
  )

  res
}

#' Summarize the taxa for the specific rank
#' @noRd
.summarize_taxa_level <- function(ps,
                                  rank = 6,
                                  norm = 1000000,
                                  absolute = FALSE,
                                  sep = "|") {
  if (!absolute) {
    ps <- microbiome::transform(ps, "compositional")
  }

  # norm the abundance data
  if (norm > 0) {
    ps@otu_table <- ps@otu_table*norm
  }

  otus <- otu_table(ps)
  otus_extend <- slot(otus, ".Data") %>%
    tibble::as_tibble()

  taxas <- tax_table(ps)@.Data %>%
    tibble::as_tibble()

  consensus <- taxas[, 1:rank]  %>%
    purrr::pmap_chr(paste, sep = sep)
  otus_extend$consensus <- consensus

  taxa_summarized <- dplyr::group_split(otus_extend, consensus) %>%
    purrr::map(.sum_consensus) %>%
    do.call(rbind, .)
  # filter taxa of which abundance is zero
  ind <- rowSums(taxa_summarized) != 0
  taxa_summarized <- taxa_summarized[ind, ] %>%
    tibble::rownames_to_column(var = "taxa") %>%
    dplyr::arrange(taxa)

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
utils::globalVariables(c(".", "taxa", "as_tibble"))
