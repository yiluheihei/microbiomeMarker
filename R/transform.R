#' Transform the taxa abundances in `otu_table` sample by sample
#'
#' Transform the taxa abundances in `otu_table` sample by sample, which means
#' the counts of each sample will be transformed individually.
#'
#' @param object  [`otu_table-class`], [`phyloseq-class`], or
#'   [`microbiomeMarker-class`].
#' @param transform transformation to apply, the options inclulde:
#' * "identity", return the original data without any transformation.
#' * "log10", the transformation is `log10(object)`, and if the data contains
#'   zeros the transformation is `log10(1 + object)`.
#' * "log10p", the transformation is `log10(1 + object)`.
#' @importFrom phyloseq t otu_table<-
#' @return A object matches the class of argument `object` with the transformed
#'   `otu_table`.
#' @export
#' @seealso [`abundances()`]
#' @examples
#' data(oxygen)
#' x1 <- transform_abundances(oxygen)
#' head(otu_table(x1), 10)
#' x2 <- transform_abundances(oxygen, "log10")
#' head(otu_table(x2), 10)
#' x3 <- transform_abundances(oxygen, "log10p")
#' head(otu_table(x3), 10)
transform_abundances <- function(object,
    transform = c("identity", "log10", "log10p")) {
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    otu <- as(otu_table(object), "matrix")

    if (transform == "identity") {
        abd <- otu
    } else if (transform == "log10") {
        abd <- transform_log10(otu)
    } else {
        abd <- transform_log10p(otu)
    }

    otu_table(object) <- otu_table(abd, taxa_are_rows = taxa_are_rows(object))

    object
}

# the data is transformed using log10(1 + x) if the data contains zeroes
transform_log10 <- function(x) {
    if (min(x) == 0) {
        warning("OTU table contains zeroes. Using log10(1 + x) instead.")
        x_norm <- log10(1 + x)
    } else {
        x_norm <- log10(x)
    }

    x_norm
}

# the data is transformed using log10(1 + x)
transform_log10p <- function(x) {
    log10(1 + x)
}
