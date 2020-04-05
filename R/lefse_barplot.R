# bar plot ----------------------------------------------------------------

#' lefse bar plot
#'
#' @param mm a [microbiomeMarker-class] object
#' @param label_level integer, number of label levels to be displayed, default
#'   `1`, `0` means display the full name of the feature
#' @param max_label_len integer, maximum number of characters of feature labels,
#' default `60`
#' @param direction the direction of bar, horizontal (`h`) or vertical (`v`),
#' default `h`
#' @importFrom ggplot2 ggplot aes geom_col labs scale_x_continuous theme_bw
#' scale_y_discrete coord_flip guide_axis
#' @return a ggplot project
#' @export
lefse_barplot <- function(mm,
                          label_level = 1,
                          max_label_len = 60,
                          direction = c("h", "v")) {
  lefse_out <- mm@marker_table
  lefse_out$feature <- factor(lefse_out$feature, rev(lefse_out$feature))
  nms_check <- any(c("feature", "lda", "enrich_group") %in% names(lefse_out))
  if (!nms_check) {
    stop("`lefse_out` must contains variable `feature`, `lda` and `enrich_group`")
  }
  direction <- match.arg(direction)

  p <-
    ggplot(
      lefse_out,
      aes(.data$lda, .data$feature, fill = .data$enrich_group)
    ) +
    geom_col() +
    labs(x = "LDA score (log10)", y = NULL, fill = NULL) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw()
  if (direction == "h") {
     p <- p +
       scale_y_discrete(
         labels = function(x)
           purrr::map_chr(x, ~ get_feature_label(.x, label_level, max_label_len))
       )
  } else {
    p <- p +
      scale_y_discrete(
        labels = function(x)
          purrr::map_chr(x, ~ get_feature_label(.x, label_level, max_label_len)),
        guide = guide_axis(angle = -90)
      ) +
      coord_flip()

  }
  p
}

#' get the feature label
#' @noRd
get_feature_label <- function(feature,
                              label_level = 1,
                              max_label_len = 60,
                              sep = "|") {
  if (length(feature) != 1) {
    stop("`feature` muste be a character vector of length 1")
  }
  if (label_level == 0) {
    feature <- feature
  } else {
    feature <- strsplit(feature, split = sep, fixed = TRUE) %>%
      unlist() %>%
      rev()
    feature_level <- length(feature)
    feature <- ifelse(
      label_level > feature_level,
      paste(rev(feature[1:feature_level]), collapse = sep),
      paste(rev(feature[1:label_level]), collapse = sep)
    )
  }

  feature_len <- nchar(feature)
  if (feature_len > max_label_len) {
    feature_letters <- unlist(strsplit(feature, ""))
    feature <- paste(
      paste(feature_letters[1:(max_label_len/2-2)], collapse = ""),
      "..",
      paste(
        feature_letters[(feature_len - max_label_len/2 + 3):feature_len],
        collapse = ""
      ),
      sep = ""
    )
  }

  feature
}

