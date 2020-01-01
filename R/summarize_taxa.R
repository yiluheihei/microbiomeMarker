#' Summarize taxa into a taxonomic level within each sample
#'
#' Provids summary information of the representation of a taxonomic levels within
#' each sample.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param level integer, taxonomic level to summarize by, default 6.
#' @param absolute logical, whether return the absolute abundance or not, default
#' FALSE.
#' @param sep a character string to separate the taxonomic levels.
#'
#' @return a data frame, each row represnets a taxa, where each col represents
#' the taxa abunance of each sample
#' @export
summarize_taxa <- function(ps, level = 6, absolute = FALSE, sep = "|") {
  if (!absolute) {
    ps <- microbiome::transform(ps, "compositional")
  }

  otus <- otu_table(ps)
  taxas <- tax_table(ps)

  consensus <- taxas[, 1:level] %>%
    slot(".Data") %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    purrr::map_df(~ ifelse(is.na(.x), "Other", .x)) %>%
    purrr::pmap_chr(paste, sep = sep)
  otus_extend <- slot(otus, ".Data") %>%
    as.data.frame(stringsAsFactors = FALSE)
  otus_extend$consensus <- consensus

  summarized_taxa <- dplyr::group_split(otus_extend, consensus) %>%
    purrr::map(sum_consensus) %>%
    do.call(rbind, .) %>%
    tibble::rownames_to_column(var = "taxa") %>%
    dplyr::arrange(taxa)
  summarized_taxa
}


#' sum all otus which belongs to the same taxa
#' @noRd
sum_consensus <- function(x) {
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
utils::globalVariables(c(".", "taxa"))
