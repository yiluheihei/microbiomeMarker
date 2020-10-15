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
