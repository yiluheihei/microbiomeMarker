#' microbiomeMarker: A package for microbiome biomarker discovery
#'
#' The microboimeMarker package provides several methods to identify micribome
#' biomarker, such as lefse, deseq2.
#'
#' @docType package
#' @name microbiomeMarker-package
#' @aliases microbiomeMarker-package
#' @importFrom dplyr %>% group_by summarise filter select bind_rows
#' group_split arrange slice mutate desc group_modify ungroup
#' @importFrom phyloseq sample_data otu_table taxa_are_rows
#' transform_sample_counts tax_table taxa_sums prune_taxa phyloseq
#' phy_tree ntaxa nsamples refseq
#' @importClassesFrom phyloseq phyloseq
#' @importFrom rlang .data
#' @importFrom methods setClass setGeneric setMethod
#' @keywords internal
NULL
