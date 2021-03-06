% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\name{microbiomeMarker}
\alias{microbiomeMarker}
\title{Build microbiomeMarker-class objects}
\usage{
microbiomeMarker(marker_table, tax_table_orig, ...)
}
\arguments{
\item{marker_table}{a \code{\linkS4class{marker_table}} object
differtial analysis}

\item{tax_table_orig}{a character vector, representing the summarized taxa}

\item{...}{arguments passed to \code{\link[phyloseq:phyloseq]{phyloseq::phyloseq()}}}
}
\description{
This the constructor to build the \code{\linkS4class{microbiomeMarker}} object, don't use
the \code{new()} constructor.
}
\references{
\href{https://stackoverflow.com/questions/9900134/is-it-bad-practice-to-access-s4-objects-slots-directly-using/9900822#9900822}{Is it bad practice to access S4 objects slots directly using @?}
}
\seealso{
\code{\link[phyloseq:phyloseq]{phyloseq::phyloseq()}},
}
