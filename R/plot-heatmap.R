#' Heatmap of microbiome marker
#'
#' Display the microbiome marker using heatmap, in which rows represents the
#' marker and columns represents the samples.
#'
#' @inheritParams plot_abundance
#' @param transform  transformation to apply, for more details see
#' [`transform_abundances()`]:
#'   * "identity", return the original data without any transformation.
#'   * "log10", the transformation is `log10(object)`, and if the data contains
#'     zeros the transformation is `log10(1 + object)`.
#'   * "log10p", the transformation is `log10(1 + object)`.
#' @param cluster_marker,cluster_sample logical, controls whether to perform
#'  clustering in markers (rows) and samples (cols), default `FALSE`.
#' @param sample_label logical, controls whether to show the sample labels in
#'   the heatmap, default `FALSE`.
#' @param ... extra arguments passed to [`ComplexHeatmap::Heatmap()`].
#' @export
#' @seealso [`transform_abundances`],[`ComplexHeatmap::Heatmap()`]
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @return a [`ComplexHeatmap::Heatmap-class`] object.
#' @examples
#' data(kostic_crc)
#' mm_lefse <- run_lefse(
#'   kostic_crc,
#'   wilcoxon_cutoff = 0.01,
#'   group = "DIAGNOSIS",
#'   kw_cutoff = 0.01,
#'   multigrp_strat = TRUE,
#'   lda_cutoff = 4
#' )
#' plot_heatmap(mm_lefse, group = "DIAGNOSIS")
plot_heatmap <- function(mm,
                         transform = c("log10", "log10p", "identity"),
                         cluster_marker = FALSE,
                         cluster_sample = FALSE,
                         markers = NULL,
                         label_level = 1,
                         max_label_len = 60,
                         sample_label = FALSE,
                         group,
                         ...) {
  stopifnot(inherits(mm, c("microbiomeMarker", "marker_table")))
  stopifnot(is.logical(sample_label))
  stopifnot(is.logical(cluster_marker))
  stopifnot(is.logical(cluster_sample))

  transform <- match.arg(transform, c("log10", "log10p", "identity"))
  mm <- transform_abundances(mm, transform = transform)

  marker <- marker_table(mm)
  # maker subset white a function here
  if (!is.null(markers)) {
    ind <- match(markers, marker$feature)
    ind_na <- is.na(ind)
    if (all(ind_na)) {
      stop("all the elements of `markers` should be a contained in marker_table")
    }

    if (any(ind_na)) {
      warning(
        paste(marker[ind_na], collapse = " "),
        "are not contained in the `marker_table`"
      )
    }

    marker <- marker[ind, ]
    # reorder to keep in increase order of effect size
    marker <- marker[order(marker$effect_size), ]
  }

  abd <- as(otu_table(mm), "matrix")
  marker_abd <- abd[match(marker$feature, row.names(abd)), ] %>%
    as.data.frame()

  # set the feature label
  labels <- get_features_labels(
    row.names(marker_abd),
    label_level,
    max_label_len
  )
  row.names(marker_abd) <- labels

  groups <- sample_data(mm)[[group]]
  group_lvl <- unique(groups)
  idx <- lapply(group_lvl, function(x)which(groups == x))
  marker_abd <- marker_abd[, unlist(idx)]
  column_nms <- rep(group_lvl, times = lengths(idx))

  if (!sample_label) {
    colnames(marker_abd) <- NULL
  }

  if (transform == "identity") {
    lgd_title <- "Abundance"
  } else {
    lgd_title <- "Log10 Abundance"
  }


  p <- Heatmap(
    as.matrix(marker_abd),
    cluster_rows = cluster_marker,
    cluster_columns = cluster_sample,
    top_annotation = HeatmapAnnotation(Group = column_nms),
    name = lgd_title,
    ...
  )

  p
}
