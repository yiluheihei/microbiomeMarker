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
#' @param norm logical, indicating whether or not to return the normalized
#'   taxa abundances.
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
  function(object,
           transform = c("identity", "log10", "log10p"),
           norm = FALSE) {
    standardGeneric("abundances")
  }
)

# otu_table object
#' @aliases abundances, otu_table-method
#' @rdname abundances-methods
setMethod(abundances, "otu_table",
  function(object,
           transform = c("identity", "log10", "log10p"),
           norm = FALSE) {
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    obj_transed <- transform_abundances(object, transform = transform)
    abd <- as(otu_table(obj_transed), "matrix")

    if (norm) {
      nf <- get_norm_factors(object)
      if (is.null(nf)) {
        # warning("return the original feature table.", call. = FALSE)
      } else {
        abd <- sweep(abd, 2, nf, "/")
      }
    }

    abd
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
  function(object,
           transform = c("identity", "log10", "log10p"),
           norm = FALSE) {
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    otu <- otu_table(object)
    if (norm) {
      nf <- get_norm_factors(object)
      if (!is.null(nf)) {
        attr(otu, "norm_factor") <- nf
      }
    }
    otu <- abundances(otu, transform = transform, norm = norm)

    otu
  }
)

# microbiomeMarker object
#' @aliases abundances,microbiomeMarker-method
#' @rdname abundances-methods
setMethod(abundances, "microbiomeMarker",
  function(object,
           transform = c("identity", "log10", "log10p")) {
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    otu <- otu_table(object)
    otu <- abundances(otu, transform = transform, norm = FALSE)

    otu
  }
)
