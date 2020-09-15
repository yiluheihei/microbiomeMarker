# marker_table class ------------------------------------------------------

#' The S4 class for storing microbiome marker information
#'
#' This Class is inherit from `data.frame`. Rows represent the microbiome
#' markers and variables represents feature of the marker.
#'
#' @name marker_table-class
#' @aliases marker_table-class
#' @field names,row.names a character vector, inherited from the input
#' data.frame
#' @field .data a list, each element corresponding the each column of the input
#' data.frame
#' @field .S3Class character, the S3 class `marker_table` inherited from:
#' "`data.frame`"
#' @author Yang Cao
#' @exportClass marker_table
setClass("marker_table", contains = "data.frame")

# validator of marker_table
validity_marker_table <- function(object) {
  msg <- NULL
  if (!"feature" %in% names(object)) {
    msg <- c(
      msg,
      "marker table must contain variable `feature` to save the name of marker"
    )
  }
  if (any(dim(object) == 0)) {
    msg <- c(msg, "marker table must have non-zero dimensions")
  }

  if (length(msg)) {
    return(msg)
  } else {
    return(TRUE)
  }
}

setValidity("marker_table", validity_marker_table)

# microbiomeMarker class --------------------------------------------------

#' The main class for microbiomeMarker data
#'
#' `microbiomeMarker-class` is inherited from the [`phyloseq::phyloseq-class`]
#' by adding a custom slot `microbiome_marker` to save the differential analysis
#' results. And it provides a seamless interface with **phyloseq**, which makes
#' **microbiomeMarker** simple and easy to use. For more details on see the
#' document of [`phyloseq::phyloseq-class`].
#' @name microbiomeMarker-class
#' @aliases microbiomeMarker-class
#' @importClassesFrom phyloseq phyloseq
#' @slot marker_table a data.frame, the result of microbiome
#'  differential analysis
#' @seealso [`phyloseq::phyloseq-class`]
#' @exportClass microbiomeMarker
`microbiomeMarker-class` <- setClass("microbiomeMarker",
  slots = c(
    marker_table = "marker_table"
  ),
  contains = "phyloseq",
  prototype = list(marker_table = NULL)
)

#' Build microbiomeMarker-class objects
#'
#' This the constructor to build the [`microbiomeMarker-class`] object, don't use
#' the `new()` constructor.
#' @param marker_table a [marker_table-class] object
#' differtial analysis
#' @param ... arguments passed to [phyloseq::phyloseq()]
#' @seealso [phyloseq::phyloseq()],
#' @references [Is it bad practice to access S4 objects slots directly using @?](https://stackoverflow.com/questions/9900134/is-it-bad-practice-to-access-s4-objects-slots-directly-using/9900822#9900822)
#' @name microbiomeMarker
#' @export
microbiomeMarker <- function(marker_table, ...) {
  ps_slots <- list(...)
  msg <- "slot `otu_table` and `tax_table` are required"
  if (!length(ps_slots)) {
    stop(msg)
  }

  ps <- phyloseq(...)
  if (!inherits(ps, "phyloseq")) {
    stop(msg)
  }

  `microbiomeMarker-class`(marker_table = marker_table, ps)
}

# validity for microbiomeMarker, at least contains two slots: otu_table,
#  tax_table
validity_microbiomeMarker <- function(object) {
  msg <- NULL
  # if (is.null(object@microbiome_marker)) {
  #   msg <- c(msg, "microbiome_marker slot is required")
  # }
  otu <- object@otu_table
  tax <- object@tax_table
  marker <- object@marker_table

  # marker in marker_table must be contained in tax_table
  if (!all(marker$feature %in% tax)) {
    msg <- c(msg, "marker in marker_table must be contained in tax_table")
  }

  if (is.null(otu)) {
    msg <- c(msg, "slot `otu_table` is required")
  }
  if (is.null(tax)) {
    msg <- c(msg, "slot `tax_table` is required")
  }
  if (nrow(otu) != nrow(tax)) {
    msg <- c(
      msg,
      "nrow of `otu_table` should be equal to the nrow of `tax_table`"
    )
  }
  if (nrow(marker) > nrow(otu)) {
    msg <- c(
      msg,
      paste0(
        "The number of different feature must be smaller than the",
        " total number of feature"
      )
    )
  }

  if (length(msg)) {
    return(msg)
  } else {
    return(TRUE)
  }
}

setValidity("microbiomeMarker", validity_microbiomeMarker)

# postHocTest  class ------------------------------------------------------

#' The postHocTest Class, represents the result of post-hoc test result among
#' multiple groups
#'
#' @slot result  a [`IRanges::DataFrameList-class`], each `DataFrame` consists of
#' five variables:
#' * `comparisons`: character, specify which two groups to test (the group names
#'   are separated by "_)
#' * `diff_mean_prop`: numeric, difference in mean proportions
#' * `pvalue`: numeric, p values
#' * `ci_lower_prop` and `ci_upper_prop`: numeric, lower and upper confidence
#' interval of difference in mean proportions
#' @slot abundance_proportion abundance proportion of each feature in each group
#' @slot conf_level confidence level
#' @slot method method used for post-hoc test
#' @slot method_str method illustration
#' @name postHocTest-class
#' @aliases postHocTest-class
#' @author Yang Cao
#' @exportClass postHocTest
#' @importClassesFrom IRanges DataFrameList
setClass("postHocTest",
  slots = c(
    result = "DataFrameList",
    abundance_proportion = "data.frame",
    conf_level = "numeric",
    method = "character",
    method_str = "character"
  ),
  prototype = list(
    result = NULL,
    conf_level = NULL,
    method = NULL,
    method_str = NULL
  )
)

# validity for postHocTest
validity_postHocTest <- function(object) {
  msg <- NULL

  # result <- object@result
  # diff_var <- setdiff(
  #   c("comparisons", "diff_means", "pvalue", "ci_lower", "ci_upper"),
  #   names(result)
  # )
  # if (length(diff_var) > 0) {
  #   msg <- c(
  #     msg,
  #     paste0("`", diff_var, "`", collapse = ", ")
  #   )
  # }

  conf_level <- object@conf_level
  if (!is.numeric(conf_level) || conf_level < 0 || conf_level > 1) {
    msg <- c(
      msg,
      "conf_level must in the range of (0,1)"
    )
  }

  method <- object@method
  if (!method %in% c("tukey", "games_howell", "scheffe", "welch_uncorrected")) {
    msg <- c(
      msg,
      "method must be one of tukey, games_howell, scheffe or welch_uncorrected"
    )
  }

  if (length(msg)) {
    return(msg)
  } else {
    return(TRUE)
  }
}

setValidity("postHocTest", validity_postHocTest)


#' Build postHocTest
#' @noRd
postHocTest <- function(result,
                        abundance_proportion,
                        conf_level = 0.95,
                        method = "tukey",
                        method_str = paste("Posthoc multiple comparisons of means: ", method)) {
  new(
    "postHocTest",
    result = result,
    abundance_proportion = abundance_proportion,
    conf_level = conf_level,
    method = method,
    method_str = method_str
  )
}
