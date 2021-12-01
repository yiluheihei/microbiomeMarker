#' Import function to read the output of picrust2 as phyloseq object
#' 
#' Import the output of picrust2 into phyloseq object
#' 
#' @param feature_tab character, file path of the prediction abundance table of 
#'   functional feature.
#' @param sam_tab character, file path of the sample meta data.
#' @param trait character, options are picrust2 function traits (including 
#'   "COG", "EC", "KO", "PFAM", "TIGRFAM", and "PHENO") and "PATHWAY".
#' @importFrom utils read.delim
#' @export
#' @return [`phyloseq::phyloseq-class`] object.
#' @details 
#' [PICRUSt2](https://huttenhower.sph.harvard.edu/picrust/) is a software for 
#' predicting abundances of functional profiles based on marker gene sequencing 
#' data. The functional profiles can be predicted from the taxonomic 
#' profiles using PICRUSt2. "Function" usually refers to gene families such as 
#' KEGG orthologs and Enzyme Classification numbers, but predictions can be 
#' made for any arbitrary trait.
#' 
#' In the `phyloseq object`, the predicted function abundance profile is stored 
#' in `otu_table` slot. And the functional trait is saved in `tax_table` slot, 
#' if the descriptions of function features is not added to the predicted table, 
#' `tax_table` will have only one rank `Picrust_trait` to represent the function 
#' feature id, or if the desciptions are added one more rank 
#' `Picrust_description` will be added to represent the description of 
#' function feature. 
#' @examples 
#' sam_tab <- system.file(
#'     "extdata", "picrust2_metadata.tsv",
#'     package = "microbiomeMarker")
#' feature_tab <- system.file(
#'     "extdata", "path_abun_unstrat_descrip.tsv.gz",
#'     package = "microbiomeMarker")
#' ps <- import_picrust2(feature_tab, sam_tab, trait = "PATHWAY")
#' ps
import_picrust2 <- function(feature_tab, 
    sam_tab = NULL, 
    trait = c("PATHWAY", "COG", "EC", "KO", "PFAM", "TIGRFAM", "PHENO")) {
    trait <- match.arg(
        trait, 
        c("COG", "EC", "KO", "PFAM", "TIGRFAM", "PHENO", "PATHWAY"))
    feature_tab <- read.delim(gzfile(feature_tab))
   
    # extract meta data of trait and stored it as taxonomy table in phyloseq
    feature_tab_vars <- names(feature_tab)
    tax_var <- feature_tab_vars[1]
    tax_rank <- "Picrust_trait"
    desp <- "description"
    if (desp %in% names(feature_tab)) {
        tax_var <- c(tax_var, desp)
        tax_rank <- c(tax_rank, paste0("Picrust_", desp))
    }
    tax_mat <- as.matrix(feature_tab[tax_var])
    name_prefix <- ifelse(trait == "PATHWAY", "path", "func")
    n_feature <- nrow(feature_tab)
    feature_nms <- paste0(name_prefix, seq_len(n_feature))
    row.names(tax_mat) <- tax_mat[, 1]
    colnames(tax_mat) <- tax_rank
    tax_tab <- tax_table(tax_mat)
    
    # sample data
    if (!is.null(sam_tab)) {
        sam_tab <- sample_data(read.delim(sam_tab))
    }
    
    # set the names of the feature
    feature_tab <- feature_tab[setdiff(feature_tab_vars, tax_var)]
    row.names(feature_tab) <- tax_mat[, 1]
    colnames(feature_tab) <- rownames(sam_tab)
    
    ps <- phyloseq(otu_table(as.matrix(feature_tab), taxa_are_rows = TRUE),
             tax_tab,
             sam_tab)
    ps
}
