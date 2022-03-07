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
#' @field .data a list, each element corresponding the each column of the
#' input data.frame
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
            "marker table must contain variable `feature`: the name of marker"
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

################################################################################
# A class may be defined as the union of other classes; that is, as a virtual
# class defined as a superclass of several other classes. This is a way of
# dealing with the expected scenarios in which one or more of the slot is not
# available, in which case NULL will be used instead.
################################################################################
#' @importClassesFrom phyloseq taxonomyTable
#' @keywords internal
setClassUnion("marker_tableOrNULL", c("marker_table", "NULL"))
#' @keywords internal
setClassUnion("taxonomyTableOrNULL", c("taxonomyTable", "NULL"))
#' @keywords internal
setClassUnion("characterOrNULL", c("character", "NULL"))
#' @keywords internal
setClassUnion("numericOrNULL", c("numeric", "NULL"))

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
#' @slot marker_table a data.frame, a [`marker_table-class`] object.
#' @slot norm_method character, method used to normalize the input `phyloseq`
#'   object.
#' @slot diff_method character, method used for marker identification.
#' @seealso [`phyloseq::phyloseq-class`], [`marker_table-class`],
#' [summarize_taxa()]
#' @exportClass microbiomeMarker
#' @return a [`microbiomeMarker-class`] object.
`microbiomeMarker-class` <- setClass("microbiomeMarker",
    slots = c(
        marker_table = "marker_tableOrNULL",
        norm_method = "characterOrNULL",
        diff_method = "characterOrNULL"
    ),
    contains = "phyloseq",
    prototype = list(
        marker_table = NULL,
        norm_method = NULL,
        diff_method = NULL
    )
)

#' Build microbiomeMarker-class objects
#'
#' This the constructor to build the [`microbiomeMarker-class`] object, don't
#' use the `new()` constructor.
#' @param marker_table a [`marker_table-class`] object differtial analysis.
#' @param norm_method character, method used to normalize the input `phyloseq`
#'   object.
#' @param diff_method character, method used for microbiome marker
#'   identification.
#' @param ... arguments passed to [phyloseq::phyloseq()]
#' @seealso [phyloseq::phyloseq()]
#' @name microbiomeMarker
#' @export
#' @return  a [`microbiomeMarker-class`] object.
#' @examples
#' microbiomeMarker(
#'     marker_table = marker_table(data.frame(
#'         feature = c("speciesA", "speciesB"),
#'         enrich_group = c("groupA", "groupB"),
#'         ef_logFC = c(-2, 2),
#'         pvalue = c(0.01, 0.01),
#'         padj = c(0.01, 0.01),
#'         row.names = c("marker1", "marker2")
#'     )),
#'     norm_method = "TSS",
#'     diff_method = "DESeq2",
#'     otu_table = otu_table(matrix(
#'         c(4, 1, 1, 4),
#'         nrow = 2, byrow = TRUE,
#'         dimnames = list(c("speciesA", "speciesB"), c("sample1", "sample2"))
#'     ),
#'     taxa_are_rows = TRUE
#'     ),
#'     tax_table = tax_table(matrix(
#'         c("speciesA", "speciesB"),
#'         nrow = 2,
#'         dimnames = list(c("speciesA", "speciesB"), "Species")
#'     )),
#'     sam_data = sample_data(data.frame(
#'         group = c("groupA", "groupB"),
#'         row.names = c("sample1", "sample2")
#'     ))
#' )
microbiomeMarker <- function(marker_table = NULL,
    norm_method = NULL,
    diff_method = NULL,
    ...) {
    ps_slots <- list(...)
    ps_component_cls <- vapply(ps_slots, class, character(1))
    if (!"otu_table" %in% ps_component_cls) {
        stop("otu_table is required")
    }
    if (!"taxonomyTable" %in% ps_component_cls) {
        stop("tax_table is required")
    }

    # set the rownmaes of marker_table as "markern"
    if (!is.null(marker_table)) {
        rownames(marker_table) <- paste0("marker", seq_len(nrow(marker_table)))
    }

    new(
        "microbiomeMarker",
        marker_table = marker_table,
        norm_method = norm_method,
        diff_method = diff_method,
        ...
    )
}

# validity for microbiomeMarker, at least contains two slots: otu_table,
#  tax_table
#' @importMethodsFrom phyloseq taxa_names
validity_microbiomeMarker <- function(object) {
    msg <- NULL
    otu <- object@otu_table
    tax <- object@tax_table
    marker <- object@marker_table
    norm_method <- object@norm_method
    diff_method <- object@diff_method

    # summarized taxa
    if (is.null(tax)) {
        msg <- c(msg, "tax_table is required")
    }

    if (is.null(otu)) {
        msg <- c(msg, "otu_table is required")
    }

    # marker in marker_table must be contained in tax_table
    if (!is.null(marker) && !is.null(tax) &&
        !all(marker$feature %in% tax@.Data[, 1])) {
        msg <- c(msg, "marker in marker_table must be contained in tax")
    }

    if (!is.null(otu) && !is.null(tax) && nrow(otu) != nrow(tax)) {
        msg <- c(
            msg,
            "nrow of `otu_table` must be equal to the length of `tax_table()`"
        )
    }

    if (!is.null(tax) && !is.null(marker) && nrow(marker) > nrow(tax)) {
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
#' @slot result  a [`IRanges::DataFrameList-class`], each `DataFrame` consists
#' of five variables:
#' * `comparisons`: character, specify which two groups to test (the group names
#'   are separated by "_)
#' * `diff_mean`: numeric, difference in mean abundances
#' * `pvalue`: numeric, p values
#' * `ci_lower` and `ci_upper`: numeric, lower and upper confidence interval of
#'   difference in mean abundances
#' @slot abundance abundance of each feature in each group
#' @slot conf_level confidence level
#' @slot method method used for post-hoc test
#' @slot method_str method illustration
#' @name postHocTest-class
#' @aliases postHocTest-class
#' @author Yang Cao
#' @exportClass postHocTest
#' @importClassesFrom IRanges DataFrameList
#' @return a [`postHocTest-class`] object.
setClass("postHocTest",
    slots = c(
        result = "DataFrameList",
        abundance = "data.frame",
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
    conf_level <- object@conf_level
    if (!is.numeric(conf_level) || conf_level < 0 || conf_level > 1) {
        msg <- c(
            msg,
            "conf_level must in the range of (0,1)"
        )
    }

    method <- object@method
    if (!method %in% 
            c("tukey", "games_howell", "scheffe", "welch_uncorrected")) {
        msg <- c(
            msg,
            paste(
                "method must be one of tukey, games_howell, scheffe or",
                "welch_uncorrected"
            )
        )
    }

    if (length(msg)) {
        return(msg)
    } else {
        return(TRUE)
    }
}

setValidity("postHocTest", validity_postHocTest)


#' Build postHocTest object
#'
#' This function is used for create `postHocTest` object, and is only used for
#' developers.
#'
#' @param result a [`IRanges::SimpleDFrameList-class`] object.
#' @param abundance data.frame.
#' @param conf_level numeric, confidence level.
#' @param method character, method for posthoc test.
#' @param  method_str character, illustrates which method is used for posthoc
#'   test.
#' @return a [`postHocTest-class`] object.
#' @export
#' @examples
#' require(IRanges)
#' pht <- postHocTest(
#'     result = DataFrameList(
#'         featureA = DataFrame(
#'             comparisons = c("group2-group1", 
#'                 "group3-group1", 
#'                 "group3-group2"),
#'             diff_mean = runif(3),
#'             pvalue = rep(0.01, 3),
#'             ci_lower = rep(0.01, 3),
#'             ci_upper = rep(0.011, 3)
#'         ),
#'         featureB = DataFrame(
#'             comparisons = c("group2-group1", 
#'                 "group3-group1", 
#'                 "group3-group2"),
#'             diff_mean = runif(3),
#'             pvalue = rep(0.01, 3),
#'             ci_lower = rep(0.01, 3),
#'             ci_upper = rep(0.011, 3)
#'         )
#'     ),
#'     abundance = data.frame(
#'         featureA = runif(3),
#'         featureB = runif(3),
#'         group = c("group1", "group2", "grou3")
#'     )
#' )
#' pht
postHocTest <- function(result,
    abundance,
    conf_level = 0.95,
    method = "tukey",
    method_str =
        paste(
            "Posthoc multiple comparisons of means: ",
            method
        )) {
    new(
        "postHocTest",
        result = result,
        abundance = abundance,
        conf_level = conf_level,
        method = method,
        method_str = method_str
    )
}
