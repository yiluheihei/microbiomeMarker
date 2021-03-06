% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metagenomeSeq.R
\name{run_metagenomeseq}
\alias{run_metagenomeseq}
\title{metagenomeSeq differential analysis}
\usage{
run_metagenomeseq(
  ps,
  group_var,
  subgroup1,
  subgroup2,
  transform = c("identity", "log10", "log10p"),
  norm = "CSS",
  norm_para = list(),
  model = c("fitFeatureModel", "fitZig"),
  p_adjust = c("none", "fdr", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY"),
  pvalue_cutoff = 0.05,
  ...
)
}
\arguments{
\item{ps}{ps a \code{\link[phyloseq:phyloseq-class]{phyloseq::phyloseq}} object.}

\item{group_var}{character, the variable to set the group, must be one of
the var of the sample metadata}

\item{subgroup1, subgroup2}{character, subgroups to be compared.}

\item{transform}{character, the methods used to transform the microbial
abundance. See \code{\link[=transform_abundances]{transform_abundances()}} for more details. The
options include:
\itemize{
\item "identity", return the original data without any transformation (default).
\item "log10", the transformation is \code{log10(object)}, and if the data contains
zeros the transformation is \code{log10(1 + object)}.
\item "log10p", the transformation is \code{log10(1 + object)}.
}}

\item{norm}{the methods used to normalize the microbial abundance data. See
\code{\link[=normalize]{normalize()}} for more details.
Options include:
\itemize{
\item a integer, e.g. 1e6 (default), indicating pre-sample normalization of
the sum of the values to 1e6.
\item "none": do not normalize.
\item "rarefy": random subsampling counts to the smallest library size in the
data set.
\item "TSS": total sum scaling, also referred to as "relative abundance", the
abundances were normalized by dividing the corresponding sample library
size.
\item "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
The scaling factor is then derived using a weighted trimmed mean over the
differences of the log-transformed gene-count fold-change between the
sample and the reference.
\item "RLE", relative log expression, RLE uses a pseudo-reference calculated
using the geometric mean of the gene-specific abundances over all
samples. The scaling factors are then calculated as the median of the
gene counts ratios between the samples and the reference.
\item "CSS": cumulative sum scaling, calculates scaling factors as the
cumulative sum of gene abundances up to a data-derived threshold.
\item "CLR": centered log-ratio normalization.
}}

\item{norm_para}{arguments passed to specific normalization methods.}

\item{model}{character, which model used for differential analysis,
"fitFeatureModel" (Zero-inflated Log-Normal mixture model)" or "fitZig"
(Zero-inflated Gaussian mixture model). As the authors of \strong{metagenomeSeq}
said, we currently recommend using the zero-inflated log-normal model.}

\item{p_adjust}{method for multiple test correction, default \code{none},
for more details see \link[stats:p.adjust]{stats::p.adjust}.}

\item{pvalue_cutoff}{numeric, p value cutoff, default 0.05}

\item{...}{extra arguments passed to the model. more details see
\code{\link[metagenomeSeq:fitFeatureModel]{metagenomeSeq::fitFeatureModel()}} and \code{\link[metagenomeSeq:fitZig]{metagenomeSeq::fitZig()}}
for more details.}
}
\value{
a \code{\linkS4class{microbiomeMarker}} object.
}
\description{
Differential expression analysis based on the Zero-inflated Log-Normal
mixture model or Zero-inflated Gaussian mixture model using metagenomeSeq.
}
\author{
Yang Cao
}
