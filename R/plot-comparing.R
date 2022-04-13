#' Plotting DA comparing result
#'
#' @param x an `compareDA` object, output from [`compare_DA()`].
#' @param sort character string specifying sort method. Possibilities are
#'   "score" which is calculated as \eqn{(auc - 0.5) * power - fdr}, "auc" for
#'   area under the ROC curve, "fpr" for false positive rate, "power" for
#'   empirical power.
#' @param ... extra arguments, just ignore it.
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_boxplot facet_wrap
#' labeller labs
#'
#' @return a [`ggplot2::ggplot`] object containing 4 subplots: "auc", "fdr",
#' "power"and "score" plot.
#'
#' @export
plot.compareDA <- function(x,
                           sort = c("score", "auc", "fpr", "power"),
                           ...) {
    sort <- match.arg(sort, c("score", "auc", "fpr", "power"))

    metrics <- x$metrics
    calls <- metrics$call
    new_metrics <- metrics[c('auc', 'fpr', 'power', 'fdr')]
    metric_med <- stats::aggregate(new_metrics,
                                   by = list(call = calls),
                                   FUN = median)
    metric_med$score <- (metric_med$auc - 0.5) * metric_med$power -
        metric_med$fdr
    metric_med$method <- metrics$method[match(metric_med$call, metrics$call)]
    metric_med <- metric_med[order(metric_med$score, decreasing = TRUE), ]

    sort_metric_method <- function(df, sort) {
        is_desc <- ifelse(sort == "fdr", FALSE, TRUE)
        df <- df[order(df[[sort]], decreasing = is_desc), ]
        method_lvl <- df[["method"]]

        method_lvl
    }

    method_lvl <- sort_metric_method(metric_med, sort)
    metrics$method <- factor(metrics$method, levels = method_lvl)
    metrics$score <- (metrics$auc - 0.5) * metrics$power - metrics$fdr

    metrics <- tidyr::pivot_longer(metrics,
                        cols = tidyr::one_of("score", "fdr", "auc", "power"))

    facet_labs <- c("Score", "Area Under the Curve",
                    "Power", "False Discovery Rate")
    names(facet_labs) <- c("score", "auc", "power", "fdr")
    p <- ggplot(metrics, aes_string("method", "value")) +
        geom_point() +
        geom_boxplot() +
        facet_wrap(c("name"),
                   ncol = 1,
                   labeller = labeller(.rows = facet_labs)) +
        labs(x = NULL, y = NULL)

    p
}
