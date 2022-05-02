#' `postHocTest` plot
#'
#' Visualize the result of post-hoc test using ggplot2
#'
#' @param pht a [`postHocTest-class`] object
#' @param feature character, to plot the post-toc test result of this feature
#' @param step_increase numeric vector with the increase in fraction of total
#' height for every additional comparison to minimize overlap, default `0.12`.
#' @name plot_postHocTest
#' @return a `ggplot` object
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes labs geom_errorbar geom_point geom_boxplot
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
#' ) %>%
#'     phyloseq::subset_taxa(Phylum == "Bacteroidetes")
#' pht <- run_posthoc_test(ps, group = "Enterotype")
#' plot_postHocTest(pht, feature = "p__Bacteroidetes|g__Alistipes")
plot_postHocTest <- function(pht,
    feature,
    step_increase = 0.12) {
    abd_long <- pht@abundance %>%
        tidyr::pivot_longer(-.data$group, names_to = "feat")

    if (!is.null(feature)) {
        abd_long <- filter(abd_long, .data$feat %in% feature)
    }

    annotation <- get_sig_annotation(pht, step_increase = step_increase)

    p_box <- suppressWarnings(
        ggplot(abd_long, aes(x = .data$group, y = .data$value)) +
            geom_boxplot() +
            ggsignif::geom_signif(
                data = annotation[annotation$feature %in% feature, ],
                aes(
                    xmin = .data$xmin, xmax = .data$xmax,
                    annotations = .data$annotation, 
                    y_position = .data$y_position),
                manual = TRUE, textsize = 3, vjust = 0.2
            ) +
            labs(x = NULL, y = "Abundance")
    )

    test_res <- as.data.frame(pht@result[[feature]])
    p_test <- ggplot(test_res, aes(x = .data$comparisons)) +
        geom_errorbar(
            aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
            width = 0.2
        ) +
        geom_point(aes(y = .data$diff_mean)) +
        labs(x = NULL, y = "95% confidence intervals")

    patchwork::wrap_plots(p_box + p_test)
}


#' get the annotation used for setting the location of the the significance bars
#'
#' @param pht a [postHocTest-class] object.
#' @param step_increase numeric vector with the increase in fraction of total
#' height for every additional comparison to minimize overlap, default `0.12`.
#' @noRd
get_sig_annotation <- function(pht, step_increase = 0.12) {
    abd <- pht@abundance
    group <- abd$group
    abd <- dplyr::mutate(abd, group = NULL)
    pht_res <- pht@result


    sig_annotation <- purrr::map2_df(
        abd, as.list(pht_res),
        ~ get_sig_annotation_single(.x, .y,
            group = group,
            step_increase = step_increase
        )
    )
    sig_annotation$feature <- rep(names(abd), each = nrow(pht_res[[1]]))

    sig_annotation
}


#' get the annotation data frame of a single feature
#'
#' @param abd numeric vector, abundance of a given feature
#' @param pht_df `data.frame` or `DFrame`, post hoc test result of a given
#' feature
#' @param group character vector the same length with `abd`, the group of the
#' samples
#' @param step_increase numeric vector with the increase in fraction of total
#' height for every additional comparison to minimize overlap, default `0.12`.
#' @return a data frame with four variables, `start`, `end`, `y_position`,
#' and `annotation`.
#' @noRd
get_sig_annotation_single <- function(abd,
    pht_df,
    group,
    step_increase = 0.12) {
    if (inherits(pht_df, "DFrame")) {
        pht_df <- as.data.frame(pht_df)
    }

    y_max <- split(abd, group) %>%
        purrr::map(range) %>%
        purrr::map_dbl(2)
    y_range <- max(abd) - min(abd)

    comps <- strsplit(pht_df$comparisons, "-", fixed = TRUE)
    start <- purrr::map_chr(comps, 1)
    end <- purrr::map_chr(comps, 2)
    y_max <- purrr::map2_dbl(start, end, ~ max(y_max[.x], y_max[.y]))
    y_pos <- purrr::map_dbl(
        seq_len(3),
        ~ y_max[.x] + y_range * step_increase * (.x - 1)
    )

    annotate_df <- data.frame(
        xmin = start,
        xmax = end,
        y_position = y_pos,
        annotation = pvalue2siglevel(pht_df$pvalue),
        stringsAsFactors = FALSE
    )

    annotate_df
}


#' covert p value to significance level
#'
#' <= 0.001 "\\*\\*\\*", <= 0.01 "\\*\\*", <=0.05 "\\*", > 0.05 "NS."
#' @noRd
#'
pvalue2siglevel <- function(p) {
    p[p <= 0.001] <- "***"
    p[p <= 0.01 & p > 0.001] <- "**"
    p[p > 0.01 & p <= 0.05] <- "*"
    p[p > 0.05] <- "NS."

    p
}
