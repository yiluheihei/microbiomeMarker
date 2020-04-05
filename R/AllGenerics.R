# marker_table class -----------------------------------------------------------

#' Build or access the marker_table
#'
#' This is the recommended function for both building and accessing microbiome
#' marker table ([`marker_table-class`]).
#' @param object an object among the set of classes difined by the
#' microbiomeMarker package that contain `marker_table`
#' @export
#' @rdname marker_table-methods
setGeneric(
  "marker_table",
  function(object)standardGeneric("marker_table")
)

# build marker_table from data.frame
#' @aliases marker_table,data.frame-method
#' @rdname marker_table-methods
setMethod("marker_table", "data.frame", function(object) {
  mt <- new("marker_table", object)

  mt
})


# access the marker_table of a microbiomeMarker-class object
#' @rdname marker_table-methods
#' @aliases marker_table,microbiomeMarker-method

setMethod("marker_table", "microbiomeMarker", function(object) {
  object@marker_table

})

# microbiomeMarker class ------------------------------------------------------

# modified from the show method of phyloseq
# https://github.com/joey711/phyloseq/blob/master/R/show-methods.R#L47-L82
#' @rdname microbiomeMarker-class
#' @param object a `microbiomeMarker-class` object
#' @export
setMethod("show", "microbiomeMarker", function(object){
  cat("microbiomeMarker-class inherited from phyloseq-class", fill = TRUE)
  cat(
    "marker_table  Marker Table:      [",
    nrow(object@marker_table), "microbiome markers with",
    ncol(object@marker_table), "variables ]",
    fill = TRUE
  )

  # print otu_table (always there).
  cat(
    "otu_table()   OTU Table:         [",
    ntaxa(otu_table(object)), "taxa and ",
    nsamples(otu_table(object)), "samples ]",
    fill = TRUE
  )

  # print Sample Data if there
  if(!is.null(sample_data(object, FALSE))){
    cat(
      "sample_data() Sample Data:       [", dim(sample_data(object))[1],
      "samples by ", dim(sample_data(object))[2],
      "sample variables ]",
      fill = TRUE
    )
  }

  # print tax Tab if there
  if(!is.null(tax_table(object, FALSE))){
    cat(
      "tax_table()   Taxonomy Table:    [", dim(tax_table(object))[1],
      "taxa by", dim(tax_table(object))[2],
      "taxonomic ranks ]",
      fill = TRUE
    )
  }

  # print tree if there
  if(!is.null(phy_tree(object, FALSE))){
    cat(
      "phy_tree()    Phylogenetic Tree: [", ntaxa(phy_tree(object)),
      "tips and", phy_tree(object)$Nnode,
      "internal nodes ]",
      fill = TRUE
    )
  }

  # print refseq summary if there
  if(!is.null(refseq(object, FALSE))){
    cat(
      "refseq()      ", class(refseq(object))[1],
      ":      [", ntaxa(refseq(object)),
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
setGeneric("nmarker", function(object) standardGeneric("nmarker"))

#' @rdname nmarker-methods
#' @aliases nmarker,microbiomeMarker-method
setMethod("nmarker", "microbiomeMarker", function(object) {
  nrow(marker_table(object))
})

#' @rdname nmarker-methods
#' @aliases nmarker,marker_table-method
setMethod("nmarker", "marker_table", function(object) {
  nrow(object)
})


# get the number of classes for microbiome discorvery ---------------------

#' Get the number of classes for microbiome discovery
#' @inheritParams nmarker
#' @docType methods
#' @rdname nclass-methods
#' @return an integer
#' @export
setGeneric("nclass", function(object) standardGeneric("nclass"))

#' @rdname nclass-methods
#' @aliases nclass,microbiomeMarker-method
setMethod("nclass", "microbiomeMarker", function(object) {
  length(unique(marker_table(object)$enrich_group))
})

#' @rdname nclass-methods
#' @aliases nclass,marker_table-method
setMethod("nclass", "marker_table", function(object) {
  length(unique(object$enrich_group))
})




