#' microbiomeMarker: A package for microbiome biomarker discovery
#'
#' The microboimeMarker package provides several methods to identify micribome
#' biomarker, such as lefse and statistical analysis.
#'
#' @docType package
#' @name microbiomeMarker
#' @aliases microbiomeMarker-package
#' @importFrom dplyr %>% group_by summarise filter select bind_rows
#' group_split arrange slice mutate desc group_modify ungroup
#' @importFrom phyloseq sample_data otu_table taxa_are_rows
#' transform_sample_counts tax_table taxa_sums prune_taxa phyloseq
#' @importClassesFrom phyloseq phyloseq
#' @importFrom rlang .data
NULL
