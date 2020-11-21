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

  # ensure the ranks_prefix is contained in available_ranks
  # and in descending order
  available_prefix <- get_available_prefix(available_ranks)
  if (!all(ranks_prefix %in% available_prefix)) {
    stop("all elements of ranks_prefix must be contained in available_ranks")
  }
  ranks_prefix <- keep_prefix_desc(ranks_prefix, type = "ranks_prefix")

  tax_prefix <- purrr::map(
    tax_split,
    ~ paste0(ranks_prefix[1:length(.x)], "__", .x) %>%
      paste0(collapse = sep))
  tax_prefix <- do.call(rbind, tax_prefix)
  # row.names(tax_prefix) <- tax_prefix
  colnames(tax_prefix) <- paste0(ranks_prefix, collapse = sep)
  tax_table(ps) <- tax_table(tax_prefix)

  taxa_names(ps) <- tax_prefix[, 1]

  ps
}

# extract the first letter of taxonomic ranks as the prefixes of the ranks
get_available_prefix <- function(ranks) {
  substr(ranks, 1, 1) %>%
    tolower()
}

# keep prefix in descending order: "k" "p" "c" "o" "f" "g" "s"
keep_prefix_desc <- function(ranks_prefix, type = c("ranks", "ranks_prefix")) {
  type <- match.arg(type, choices = c("ranks", "ranks_prefix"))
  available_prefix <- get_available_prefix(available_ranks)
  idx_desc <- sort(match(ranks_prefix, available_prefix))

  if (type == "ranks") {
    return(available_ranks[idx_desc])
  } else {
    return(available_prefix[idx_desc])
  }
}
