#' dotplot of effect size of microbiomeMarker data
#'
#' dotplot of effect size microbiomeMarker data.
#'
#' @param mm a [microbiomeMarker-class] object
#' @param label_level integer, number of label levels to be displayed, default
#'   `1`, `0` means display the full name of the feature
#' @param max_label_len integer, maximum number of characters of feature labels,
#'   default `60`
#' @param direction the direction of bar, horizontal (`h`) or vertical (`v`),
#'   default `h`
#' @param n integer, number of markers to dispalay, default `10`. If the number
#'   of markers is less than 10, the actual number of markers will be displayed.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_x_continuous theme_bw
#'   scale_y_discrete coord_flip guide_axis
#' @return a ggplot project
#' @export
# plot_ef_dot <- function()
