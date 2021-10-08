#' Assign a new OTU table
#'
#' Assign a new OTU table in microbiomeMarker object
#' @param x [`microbiomeMarker-class`]
#' @param value [`otu_table-class`], [`phyloseq-class`],
#'   or [`microbiomeMarker-class`]
#' @docType methods
#' @name assign-otu_table
#' @importMethodsFrom phyloseq otu_table<-
#' @importClassesFrom phyloseq otu_table
#' @return a [`microbiomeMarker-class`] object.
NULL

#' @rdname assign-otu_table
#' @aliases otu_table<-,microbiomeMarker,otu_table-method
setMethod(
    "otu_table<-", c("microbiomeMarker", "otu_table"),
    function(x, value) {
        microbiomeMarker(
            x@marker_table,
            x@norm_method,
            x@diff_method,
            otu_table = value,
            sam_data = x@sam_data,
            phy_tree = x@phy_tree,
            refseq = x@refseq,
            tax_table = x@tax_table
        )
    }
)

#' @rdname assign-otu_table
#' @aliases otu_table<-,microbiomeMarker,phyloseq-method
setMethod(
    "otu_table<-", c("microbiomeMarker", "phyloseq"),
    function(x, value) {
        microbiomeMarker(
            x@marker_table,
            x@norm_method,
            x@diff_method,
            otu_table = otu_table(value),
            sam_data = x@sam_data,
            phy_tree = x@phy_tree,
            refseq = x@refseq,
            tax_table = x@tax_table
        )
    }
)

#' @rdname assign-otu_table
#' @aliases otu_table<-,microbiomeMarker,microbiomeMarker-method
setMethod(
    "otu_table<-", c("microbiomeMarker", "microbiomeMarker"),
    function(x, value) {
        microbiomeMarker(
            x@marker_table,
            x@norm_method,
            x@diff_method,
            otu_table = otu_table(value),
            sam_data = x@sam_data,
            phy_tree = x@phy_tree,
            refseq = x@refseq,
            tax_table = x@tax_table
        )
    }
)
