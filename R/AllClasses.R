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
#' @export
setClass("marker_table",contains = "data.frame")

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
#' @slot marker_table a  data.frame, the result of microbiome
#'  differntial analysis
#' @seealso [`phyloseq::phyloseq-class`]
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

  `microbiomeMarker-class`(ps, marker_table = marker_table)
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





