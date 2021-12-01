# marker_table class -----------------------------------------------------------

#' Build or access the marker_table
#'
#' This is the recommended function for both building and accessing microbiome
#' marker table ([`marker_table-class`]).
#' @param object an object among the set of classes defined by the
#' microbiomeMarker package that contain `marker_table`
#' @export
#' @rdname marker_table-methods
#' @return a [`marker_table-class`] object.
#' @examples
#' data(enterotypes_arumugam)
#' mm <- run_limma_voom(
#'     enterotypes_arumugam,
#'     "Enterotype",
#'     contrast = c("Enterotype 3", "Enterotype 2"),
#'     pvalue_cutoff = 0.05,
#'     p_adjust = "fdr"
#' )
#' marker_table(mm)
setGeneric(
    "marker_table",
    function(object) standardGeneric("marker_table")
)

# build marker_table from data.frame
#' @aliases marker_table,data.frame-method
#' @rdname marker_table-methods
setMethod("marker_table", "data.frame", function(object) {
    mt <- new("marker_table", object)
    row.names(mt) <- paste0("marker", seq_len(nrow(object)))

    mt
})


# access the marker_table of a microbiomeMarker-class object
#' @rdname marker_table-methods
#' @aliases marker_table,microbiomeMarker-method

setMethod("marker_table", "microbiomeMarker", function(object) {
    object@marker_table
})


# Assign marker_table -----------------------------------------------------
#' Assign marker_table to `object`
#'
#' This function replace the `marker_table` slot of `object` with `value`.
#'
#' @param object a [`microbiomeMarker-class`] object to modify.
#' @param value new value to replace the `marker_table` slot of `object`.
#'   Either a `marker_table-class`, a `data.frame` that can be coerced
#'   into `marker_table-class`.
#' @export
#' @rdname assign-marker_table
#' @aliases assign-marker_table marker_table<-
#' @return a [`microbiomeMarker-class`] object.
#' @examples
#' data(enterotypes_arumugam)
#' mm <- run_limma_voom(
#'     enterotypes_arumugam,
#'     "Enterotype",
#'     contrast = c("Enterotype 3", "Enterotype 2"),
#'     pvalue_cutoff = 0.1,
#'     p_adjust = "fdr"
#' )
#' mm_marker <- marker_table(mm)
#' mm_marker
#' marker_table(mm) <- mm_marker[1:2, ]
#' marker_table(mm)
"marker_table<-" <- function(object, value) {
    if (!inherits(value, "marker_table") && !is.null(value)) {
        value <- marker_table(value)
    }

    microbiomeMarker(
        marker_table = value,
        norm_method = object@norm_method,
        diff_method = object@diff_method,
        otu_table = object@otu_table,
        tax_table = object@tax_table,
        phy_tree = object@phy_tree,
        refseq = object@refseq
    )
}
# microbiomeMarker class ------------------------------------------------------

# modified from the show method of phyloseq
# https://github.com/joey711/phyloseq/blob/master/R/show-methods.R#L47-L82
#' @rdname microbiomeMarker-class
#' @param object a `microbiomeMarker-class` object
#' @export
setMethod("show", "microbiomeMarker", function(object) {
    cat("microbiomeMarker-class inherited from phyloseq-class", fill = TRUE)
    norm <- object@norm_method
    if (!is.null(norm)) {
        if (grepl("per-sample normalized", norm)) {
            norm <- gsub(".*to ", "", norm)
            cat(
                "normalization: per-sample to value [", norm, "]",
                fill = TRUE
            )
        } else {
            cat(
                "normalization method:              [", norm, "]",
                fill = TRUE
            )
        }
    }

    if (!is.null(object@diff_method)) {
        cat(
            "microbiome marker identity method: [",
            object@diff_method,
            "]",
            fill = TRUE
        )
    }
    
    if (!is.null(object@marker_table)) {
        cat(
            "marker_table() Marker Table:       [",
            nrow(object@marker_table), "microbiome markers with",
            ncol(object@marker_table), "variables ]",
            fill = TRUE
        )
    } else {
        cat(
            "marker_table() Marker Table:       [",
            "no microbiome markers were identified ]",
            fill = TRUE
        )
    }

    # print otu_table (always there).
    cat(
        "otu_table()    OTU Table:          [",
        ntaxa(otu_table(object)), "taxa and ",
        nsamples(otu_table(object)), "samples ]",
        fill = TRUE
    )

    # print Sample Data if there
    if (!is.null(sample_data(object, FALSE))) {
        cat(
            "sample_data()  Sample Data:        [", dim(sample_data(object))[1],
            "samples by ", dim(sample_data(object))[2],
            "sample variables ]",
            fill = TRUE
        )
    }

    # print tax Tab if there
    if (!is.null(tax_table(object, FALSE))) {
        cat(
            "tax_table()    Taxonomy Table:     [", dim(tax_table(object))[1],
            "taxa by", dim(tax_table(object))[2],
            "taxonomic ranks ]",
            fill = TRUE
        )
    }

    # print tree if there
    if (!is.null(phy_tree(object, FALSE))) {
        cat(
            "phy_tree()    Phylogenetic Tree:   [", ntaxa(phy_tree(object)),
            "tips and", phy_tree(object)$Nnode,
            "internal nodes ]",
            fill = TRUE
        )
    }

    # print refseq summary if there
    if (!is.null(refseq(object, FALSE))) {
        cat(
            "refseq()      ", class(refseq(object))[1],
            ":         [", ntaxa(refseq(object)),
            " reference sequences ]",
            fill = TRUE
        )
    }
})

# get the number of markers -----------------------------------------------

#' Get the number of microbiome markers
#' @param object a [`microbiomeMarker-class`] or [`marker_table-class`] object
#' @docType methods
#' @rdname nmarker-methods
#' @return an integer, the number of microbiome markers
#' @export
#' @examples
#' mt <- marker_table(data.frame(
#'     feature = c("speciesA", "speciesB"),
#'     enrich_group = c("groupA", "groupB"),
#'     ef_logFC = c(-2, 2),
#'     pvalue = c(0.01, 0.01),
#'     padj = c(0.01, 0.01),
#'     row.names = c("marker1", "marker2")
#' ))
#' nmarker(mt)
setGeneric("nmarker", function(object) standardGeneric("nmarker"))

#' @rdname nmarker-methods
#' @aliases nmarker,microbiomeMarker-method
setMethod("nmarker", "microbiomeMarker", function(object) {
    marker <- marker_table(object)
    ifelse(is.null(marker), 0L,  nrow(marker))
})

#' @rdname nmarker-methods
#' @aliases nmarker,marker_table-method
setMethod("nmarker", "marker_table", function(object) {
    ifelse(is.null(object), 0L,  nrow(object))
})

# postHocTest class -------------------------------------------------------
#' @rdname postHocTest-class
#' @aliases show, postHocTest-method
#' @param object a `postHocTest-class` object
#' @export
setMethod("show", "postHocTest", function(object) {
    cat("postHocTest-class object", fill = TRUE)
    result <- object@result
    var_mean <- c(
        "pair groups to test which separated by '-'",
        "difference in mean proportions",
        "post hoc test p values",
        "lower confidence interval",
        "upper confidence interval"
    )
    cat(
        "Pairwise test result of", length(result), " features, ",
        "DataFrameList object, each DataFrame has five variables:\n       ",
        paste0(
            names(result[[1]]),
            c("    : ", ": ", "        : ", " : ", " : "),
            var_mean,
            collapse = "        ",
            "\n"
        )
    )
    cat(
        "Posthoc multiple comparisons of means",
        " using ", object@method, " method",
        fill = TRUE
    )
})
