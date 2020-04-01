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
#' @slot microbiome_marker a  data.frame, the result of microbiome
#'  differntial analysis
#' @seealso [`phyloseq::phyloseq-class`]
`microbiomeMarker-class` <- setClass("microbiomeMarker",
  slots = c(
    microbiome_marker = "data.frame"
  ),
  contains = "phyloseq",
  prototype = list(microbiome_marker = NULL)
)

#' Build microbiomeMarker-class objects from their slots
#'
#' This the constructor to build the `microbiomeMarker-class` object, don't use
#' the `new()` constructor.
#' @param microbiome_marker a data frame, the result of microbiome
#' differtial analysis
#' @param ... arguments passed to [phyloseq::phyloseq()]
#' @seealso [phyloseq::phyloseq()],
#' @references [Is it bad practice to access S4 objects slots directly using @?](https://stackoverflow.com/questions/9900134/is-it-bad-practice-to-access-s4-objects-slots-directly-using/9900822#9900822)
#' @name microbiomeMarker
#' @export
microbiomeMarker <- function(microbiome_marker, ...) {
  ps_slots <- list(...)
  msg <- "slot `otu_table` and `tax_table` are required"
  if (!length(ps_slots)) {
    stop(msg)
  }

  ps <- phyloseq(...)
  if (!inherits(ps, "phyloseq")) {
    stop(msg)
  }

  `microbiomeMarker-class`(ps, microbiome_marker = microbiome_marker)
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
  marker <- object@microbiome_marker

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
