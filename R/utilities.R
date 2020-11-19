#' phyloseq quality control, remove otu/asv of which abundance is zero
#' @noRd
phyloseq_qc <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

# only first letter in lower case
upper_firstletter <- function(x){
  paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
}

#' Transpose the phyloseq object to ensure taxa are in rows
#' @param ps a [phyloseq::phyloseq-class] object
#' @importMethodsFrom phyloseq t
#' @keywords internal
keep_taxa_in_rows <- function(ps) {
  if (!taxa_are_rows(ps)) {
    ps <- t(ps)
  }

  ps
}

# add ranks prefix, e.g k__, p__
add_prefix <- function(ps) {
  tax <- as(tax_table(ps), "matrix") %>%
    as.data.frame()
  lvl <- colnames(tax)

  diff_lvl <- setdiff(lvl, availabel_ranks)
  if (length(diff_lvl) != 0) {
    stop("rank names of `ps` must be one of Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }

  prefix <- substr(lvl, 1, 1) %>%
    tolower() %>%
    paste("__", sep = "")
  tax_new <- mapply(function(x, y) paste0(x, y), prefix, tax, SIMPLIFY = FALSE)
  tax_new <- do.call(cbind, tax_new)
  row.names(tax_new) <- row.names(tax)
  colnames(tax_new) <- lvl
  tax_table(ps) <- tax_new

  ps
}

#' add prefix of taxonomic ranks for summarized data construct from original
#' lefse (galaxy lefse or python app) input, p__, k__
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param ranks character vector, prefix of ranks to add, e.g. "p", "c"
#' @importFrom phyloseq taxa_names<-
#' @noRd
add_prefix_summarized <- function(ps, ranks_prefix, sep = "|") {
  tax <- tax_table(ps)@.Data[, 1]
  tax_split <- strsplit(tax, split = sep, fixed = TRUE)

  if (max(lengths(tax_split)) != length(ranks_prefix)) {
    stop("The length of `ranks_prefix` muste be equal to number of taxonomic ranks.")
  }

  # ensure the ranks_prefix is contained in availabel_ranks
  # and in descending order
  availabel_prefix <- substr(availabel_ranks, 1, 1) %>%
    tolower()
  if (!all(ranks_prefix %in% availabel_prefix)) {
    stop("all elements of ranks_prefix must be contained in availabel_ranks")
  }
  idx_descending <- sort(match(ranks_prefix, availabel_prefix))
  ranks_prefix <- availabel_prefix[idx_descending]

  tax_prefix <- purrr::map(
    tax_split,
    ~ paste0(ranks_prefix[1:length(.x)], "__", .x) %>%
      paste0(collapse = sep))
  tax_prefix <- do.call(rbind, tax_prefix)
  # row.names(tax_prefix) <- tax_prefix
  colnames(tax_prefix) <- paste0(availabel_ranks[idx_descending], collapse = sep)
  tax_table(ps) <- tax_table(tax_prefix)

  taxa_names(ps) <- tax_prefix[, 1]

  ps
}
