#' phyloseq quality control, remove otu/asv of which abundance is zero
#' @noRd
phyloseq_qc <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

#' check whether tax abundance table is summarized or not
#' @noRd
check_tax_summarize <- function(ps) {
  is_summarize <- ifelse(
    ncol(tax_table(ps)) == 1,
    TRUE, FALSE
  )

  is_summarize
}


