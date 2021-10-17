#' plot the abundances of markers
#'
#' @inheritParams plot_ef_bar
#' @param group character, the variable to set the group
#' @return a [`ggplot2::ggplot`] object.
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw element_text
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
#' plot_abundance(mm, group = "Enterotype")
plot_abundance <- function(mm,
    label_level = 1,
    max_label_len = 60,
    markers = NULL,
    group) {
    stopifnot(inherits(mm, c("microbiomeMarker", "marker_table")))

    sample_meta <- sample_data(mm)
    sample_meta_nms <- names(sample_meta)
    if (!group %in% sample_meta_nms) {
        stop(
            "`group_var` must be one of the sample-level variables",
            call. = FALSE
        )
    }

    marker <- marker_table(mm)

    # maker subset white a function here: replicate code in .plot_ef
    if (!is.null(markers)) {
        ind <- match(markers, marker$feature)
        ind_na <- is.na(ind)
        if (all(ind_na)) {
            stop(
                "all the elements of `markers` should",
                " be a contained in marker_table",
                call. = FALSE
            )
        }

        if (any(ind_na)) {
            warning(
                paste(marker[ind_na], collapse = " "),
                "are not contained in the `marker_table`"
            )
        }

        marker <- marker[ind, ]

        # reorder to keep in increase order of effect size
        marker <- marker[order(marker[[3]]), ]
    }

    abd <- as(otu_table(mm), "matrix")
    marker_abd <- abd[match(marker$feature, row.names(abd)), ] %>%
        as.data.frame()

    groups <- sample_meta[[group]]
    names(groups) <- names(marker_abd)

    marker_abd$feature <- row.names(marker_abd)
    marker_abd <- tidyr::pivot_longer(
        marker_abd,
        -.data$feature,
        names_to = "sample",
        values_to = "abd"
    )
    marker_abd[[group]] <- groups[match(marker_abd$sample, names(groups))]

    p <- ggplot(
        marker_abd,
        aes(x = .data$abd, y = .data$feature, fill = .data[[group]])
    ) +
        geom_boxplot() +
        labs(x = "Abundance", y = NULL) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(labels = function(x) {
            get_features_labels(x, label_level, max_label_len)
        }) +
        theme_bw()

    p
}
