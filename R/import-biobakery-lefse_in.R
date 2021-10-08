#' @title Import function to read the tab-delimited input file of biobakery
#' lefse
#'
#' @description For biobakey lefse, the input file must be a tab-delimited
#' text, consists of a list of numerical features, the class vector and
#' optionally the subclass and subject vectors. The features can be read counts
#' directly or abundance floating-point values more generally, and the first
#' field is the name of the feature. Class, subclass and subject vectors have a
#' name (the first field) and a list of non-numerical strings. This function
#' requires the features are organized in rows, although both column and row
#' feature organization is accepted in biobakery lefse.
#'
#' @param file the file path of tab-delimited input file of biobakery lefse
#' @param ranks_prefix character vector, prefix of taxonomic ranks to add,
#'   e.g. "p" for "Phylum", "g" for "Genus".
#' @param meta_rows integer vector, set which rows represent the meta data,
#'   such as class, subclass and subject, default `1`.
#' @param sep character, separator between different taxnomic ranks,
#'   default `|`.
#' @noRd
#' @return a [`phyloseq::phyloseq-class`] object.
#' @examples
#' # file <- system.file(
#' #   "extdata",
#' #   "hmp_small_aerobiosis.txt",
#' #   package = "microbiomeMarker"
#' # )
#' # six level of taxonomic ranks,
#' # meta data: row 1 represents class (oxygen_availability),
#' # row 2 represents subclass (body_site),
#' # row 3 represents subject (subject_id)
#' # ps <- import_biobakery_lefse_in(
#' #   file,
#' #   ranks_prefix = c("k", "p", "c", "o", "f", "g"),
#' #   meta_rows = 1:3,
#' # )
import_biobakery_lefse_in <- function(file,
    ranks_prefix,
    meta_rows = 1,
    sep = "|") {
    dat <- utils::read.delim(file, header = FALSE)

    # meta data of samples
    meta_nms <- dat[meta_rows, 1]
    sample_meta <- dat[meta_rows, -1] %>%
        t() %>%
        as.data.frame() %>%
        sample_data()
    colnames(sample_meta) <- meta_nms
    row.names(sample_meta) <- paste0("sa", seq_len(nrow(sample_meta)))

    # tax table
    tax <- dat[-meta_rows, 1, drop = FALSE] %>%
        as.matrix() %>%
        tax_table()
    # ensure the ranks_prefix is contained in available_ranks
    # and in descending order
    available_prefix <- get_available_prefix(available_ranks)
    if (!all(ranks_prefix %in% available_prefix)) {
        stop("all elements of ranks_prefix must be contained ",
            "in available_ranks"
        )
    }
    tax_nms <- keep_prefix_desc(ranks_prefix, type = "ranks") %>%
        paste0(collapse = sep)
    colnames(tax) <- tax_nms
    row.names(tax) <- paste0("feature", seq_len(nrow(tax)))

    # otu table
    otu <- dat[-meta_rows, -1, drop = FALSE] %>%
        apply(2, as.numeric) %>%
        otu_table(taxa_are_rows = TRUE)
    row.names(otu) <- taxa_names(tax)
    colnames(otu) <- sample_names(sample_meta)

    ps <- phyloseq(otu, tax, sample_meta) %>%
        add_prefix_summarized(ranks_prefix, sep)

    ps
}
