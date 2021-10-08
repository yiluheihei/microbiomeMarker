# transpose otu_table and then convert it to data.frame
#' @importMethodsFrom phyloseq t
transpose_and_2df <- function(ot) {
    t(ot) %>% as.data.frame()
}
