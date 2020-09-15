# This function is inspired from microbiome::abundances

#' Extract taxa abundances
#'
#' Extract taxa abundances from phyloseq objects.
#'
#' @param object [`otu_table-class`], [`phyloseq-class`], or [`microbiomeMarker-class`].
#' @param transform transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#' @return abundance matrix with taxa in rows and samples in columns.
#' @seealso [`otu_table-class`], [`phyloseq-class`], [`microbiomeMarker-class`],
#' [`transform_abundances`]
#' @docType methods
#' @rdname abundances-methods
#' @export
#' @examples
#' \dontrun{
#' data(oxygen)
#' abundances(oxygen)
#' }
setGeneric("abundances",
  function(object, transform = c("identity", "log10", "log10p")) {
    standardGeneric("abundances")
  }
)

# otu_table object
#' @aliases abundances, otu_table-method
#' @rdname abundances-methods
setMethod(abundances, "otu_table",
  function(object, transform = c("identity", "log10", "log10p")) {
    transform <- match.arg(transform)
    new_obj <- transform_abundances(object, transform = transform)

    as(otu_table(new_obj), "matrix")
    # otu <- as(object, "matrix")
    #
    # # ensure taxa are on the rows
    # if (!taxa_are_rows(object) && ntaxa(object) > 1 && nsamples(object) > 1) {
    #   otu <- t(otu)
    # }
    #
    # if (ntaxa(object) == 1) {
    #   otu <- matrix(otu, nrow=1)
    #   rownames(otu) <- taxa_names(object)
    #   colnames(otu) <- sample_names(object)
    # }
    #
    # if (nsamples(object) == 1) {
    #   otu <- matrix(otu, ncol=1)
    #   rownames(otu) <- taxa_names(object)
    #   colnames(otu) <- sample_names(object)
    # }
    #
    # otu <- transform_abundances(otu, transform = transform)
    # otu
  }
)

# phyloseq object
#' @aliases abundances,phyloseq-method
#' @rdname abundances-methods
setMethod(abundances, "phyloseq",
  function(object, transform = c("identity", "log10", "log10p")) {
    transform <- match.arg(transform)
    otu <- otu_table(object)
    otu <- abundances(otu, transform = transform)

    otu
  }
)

# microbiomeMarker object
#' @aliases abundances,microbiomeMarker-method
#' @rdname abundances-methods
setMethod(abundances, "microbiomeMarker",
  function(object, transform = c("identity", "log10", "log10p")) {
    transform <- match.arg(transform)
    otu <- otu_table(object)
    otu <- abundances(otu, transform = transform)

    otu
  }
)
