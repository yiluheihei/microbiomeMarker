#' lefse bar plot
#'
#' @param lefse_out data frame, out of \code{\link{lefse}}
#' @param lda_threshold numberic, threshold of lda score, features whose lda score is less
#'   than, default is 2
#' then it will not be plotted
#' @import ggplot2
#' @importFrom dplyr arrange mutate desc
#' @return a ggplot project
#' @export
lefse_barplot <- function(lefse_out, lda_threshold = 2) {
  lefse_out_arranged <- arrange(
    lefse_out,
    .data$enrich_group,
    desc(.data$lda_score)
  ) %>%
    mutate(otu = factor(.data$otu, levels = .data$otu))

  p <-
    ggplot(
      lefse_out_arranged,
      aes(.data$otu, .data$lda_score, fill = .data$enrich_group)
    ) +
    geom_col() +
    labs(x = "Features", y = "LDA score", fill = NULL) +
    scale_x_discrete(limits = rev(lefse_out_arranged$otu)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y = "LDA SCORE (log10)") +
    coord_flip() +
    theme_bw()

  p
}

# suppress the checking notes “no visible binding for global variable", which is
# caused by NSE
# utils::globalVariables(c("otu", "lda_score"))
