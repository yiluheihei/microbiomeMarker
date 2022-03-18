#' Subset microbiome markers
#'
#' Subset markers based on an expression related to the columns and values
#' within the `marker_table` slot of `mm`.
#'
#' @param mm a [`microbiomeMarker-class`] or [`marker_table-class`] object.
#' @param ... the subsetting expression passed to [`base::subset()`].
#' @return a subset object in the same class with `mm`.
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' mm <- run_limma_voom(
#'     enterotypes_arumugam,
#'     "Enterotype",
#'     contrast = c("Enterotype 3", "Enterotype 2"),
#'     pvalue_cutoff = 0.01,
#'     p_adjust = "none"
#' )
#' subset_marker(mm, pvalue < 0.005)
subset_marker <- function(mm, ...) {
    if (is.null(marker_table(mm))) {
        warning("No marker_table in `mm`")
        return(mm)
    }

    marker <- data.frame(marker_table(mm))
    marker_new <- subset(marker, ...)

    if (inherits(mm, "marker_table")) {
        return(marker_table(marker_new))
    } else {
        marker_table(mm) <- marker_new
        return(mm)
    }
}
