#' bar and dot plot of effect size of microbiomeMarker data
#'
#' bar and dot plot of effect size microbiomeMarker data. This function returns
#' a `ggplot2` object that can be saved or further customized using **ggplot2**
#' package.
#'
#' @param mm a [microbiomeMarker-class] object
#' @param label_level integer, number of label levels to be displayed, default
#'   `1`, `0` means display the full name of the feature
#' @param max_label_len integer, maximum number of characters of feature label,
#'   default `60`
#' @param markers character vector, markers to display, default `NULL`,
#'   indicating plot all markers.
#' @importFrom ggplot2 ggplot aes geom_col labs scale_x_continuous theme_bw
#' scale_y_discrete guide_axis
#' @return a ggplot project
#' @export
#' @rdname effect_size-plot
#' @aliases ef-barplot,ef-dotplot
plot_ef_bar <- function(mm,
                        label_level = 1,
                        max_label_len = 60,
                        markers = NULL) {
  .plot_ef(mm, label_level, max_label_len, markers, "bar")
}

#' @rdname effect_size-plot
#' @export
plot_ef_dot <- function(mm,
                        label_level = 1,
                        max_label_len = 60,
                        markers = NULL) {
  .plot_ef(mm, label_level, max_label_len, markers, "dot")
}


# plot of effect size
.plot_ef <- function(mm,
                    label_level = 1,
                    max_label_len = 60,
                    markers = NULL,
                    type = c("bar", "dot")) {

  stopifnot(inherits(mm, c("microbiomeMarker", "marker_table")))
  type <- match.arg(type, c("bar", "dot"))

  marker <- marker_table(mm)

  # effect size names
  orig_ef_nm <- names(marker)[3]
  names(marker)[3] <- "effect_size"

  # labels of x
  # effect size: lda for lefse, diff_mean for classical test, logFC for
  # metagenomeSeq, DESeq2, edgeR
  if (orig_ef_nm == "lda") {
    label_x <- "LDA score (log10)"
  } else if (orig_ef_nm == "diff_mean") {
    label_x <- "Differential means"
  } else if (orig_ef_nm == "logFC") {
    label_x <- "log2 Fold Change"
  } else if (orig_ef_nm == "eta_squared") {
    label_x <- "Eta squared"
  } else {
    stop(
      "The name of third column muste be one of lda ",
      "diff_mean, eat_squared or logFC"
    )
  }

  # the levels of features: in increase order of effect size
  # marker <- marker[order(marker$effect_size), ]

  # maker subset
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
    # marker <- marker[order(marker$effect_size), ]
  }

  # increase order in each group
  marker <- dplyr::arrange(data.frame(marker), .data$enrich_group, .data$effect_size)
  feat <- marker$feature
  marker$feature <- factor(feat, levels = feat)

  nms_check <- any(c("feature", "enrich_group") %in% names(marker))
  if (!nms_check) {
    stop("`marker_table` must contains variable `feature` and `enrich_group`")

  }

  marker$logp <- -log10(marker$padj)

  if (type == "bar") {
    p <-
      ggplot(
        marker,
        aes(.data$effect_size, .data$feature, fill = .data$enrich_group)) +
      geom_col() +
      labs(x = label_x, y = NULL, fill = "Enriched group") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_discrete(labels = function(x) {
        get_features_labels(x, label_level, max_label_len)}) +
      theme_bw()
  } else {
    p <- ggplot(
      marker,
      aes(.data$effect_size, .data$feature,
          color = .data$enrich_group, size = .data$logp)) +
      geom_point() +
      labs(x = label_x, y = NULL, color = "Enriched group", size = "-log10(pvalue)") +
      scale_y_discrete(labels = function(x) {
        get_features_labels(x, label_level, max_label_len)}) +
      theme_bw()
  }

  p
}


#' get the labels of markers which will be used in the barplot
#' @noRd
get_features_labels <- function(features, label_level, max_label_len) {
  purrr::map_chr(features, ~ get_feature_label(.x, label_level, max_label_len))
}

#' get the label of a single feature
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
  # replace "Unknown" label in the species level as "sp."
  feature <- replace_unknown_species(feature)

  feature
}

# replace "Unknown" label in the species level as "sp."
replace_unknown_species <- function(feature, sep = "|") {
  species_hased <- grepl("s__", feature, fixed = TRUE)
  if (!species_hased) {
    return(feature)
  }

  taxa_lvl <- strsplit(feature, sep, fixed = TRUE)
  n_lvl <- length(taxa_lvl)
  sp <- taxa_lvl[[n_lvl]]
  sp <- gsub("Unknown", "sp.", feature, fixed = TRUE)
  taxa_lvl[[n_lvl]] <- sp
  feature <- paste(taxa_lvl, collapse = sep)

  feature
}
