#' Normalize the microbial abundance data
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param method the methods used to normalize the microbial abundance data.
#'   Options includes:
#'   * "none": do not normalize.
#'   * "rarefy": random subsampling counts to the smallest library size in the
#'     data set.
#'   * "TSS": total sum scaling, also referred to as "relative abundance", the
#'     abundances were normalized by dividing the corresponding sample library
#'     size.
#'   * "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
#'     The scaling factor is then derived using a weighted trimmed mean over the
#'     differences of the log-transformed gene-count fold-change between the
#'     sample and the reference.
#'   * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#'   * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#'   * "CLR": centered log-ratio normalization.
#' @param ... other arguments passed to the corresponding normalization methods.
#' @seealso [edgeR::calcNormFactors()],[DESeq2::estimateSizeFactorsForMatrix()],
#' [metagenomeSeq::cumNorm()]
#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,otu_table-method normalize
#' @rdname normalize-methods
setMethod("normalize", "otu_table",
  function(object,
    method = c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR"),
    ...) {
    method <- match.arg(
      method,
      c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR")
    )
    switch (method,
      none = object,
      rafefy = norm_rarefy(object, ...),
      TSS = norm_tss(object),
      TMM = norm_tmm(object, ...),
      RLE = norm_rle(object, ...),
      CSS = norm_css(object, ...),
      CLR = norm_clr(object)
    )
  })

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,phyloseq-method
#' @rdname normalize-methods
setMethod("normalize", "phyloseq",
  function(object,
    method = c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR"),
    ...) {
    otu <- otu_table(object)

    otu_table(object) <- otu_table(
      normalize(otu),
      taxa_are_rows = taxa_are_rows(object)
    )

    object
  })


#' rarefying
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class] object
#' @param size,rng_seed,replace,trim_otus,verbose for details see
#'   [phyloseq::rarefy_even_depth()]
#' @keywords internal
#' @importFrom phyloseq rarefy_even_depth sample_sums
norm_rarefy <- function(object,
  size = min(sample_sums(object)),
  rng_seed = FALSE,
  replace = TRUE,
  trim_otus = TRUE,
  verbose = TRUE) {

  object_normed <- rarefy_even_depth(
    object,
    sample.size = size,
    rngseed = rng_seed,
    replace = replace,
    trimOTUs = trim_otus,
    verbose = verbose
  )

  object_normed
}

#' TSS normalization (relative abundance)
#' @param object object a [phyloseq::phyloseq-class] or
#'   [phyloseq::otu_table-class] object
#' @keywords internal
#' @importFrom phyloseq otu_table<-
norm_tss <- function(object) {
  otu <- otu_table(object)
  size <- colSums(otu)
  otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
  otu_table(object) <- otu_table(
    otu_normed,
    taxa_are_rows = taxa_are_rows(object)
  )

  object
}

#' css normalization (metagenomeSeq)
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param p The pth quantile, numeric or a function, if `p` is a function, it
#'   will be used to calculate the pth quantile, see[metagenomeSeq::cumNorm()].
#' @param log logical, whether or not to log2 transfrom scale,
#'   see [metagenomeSeq::MRcounts()].
#' @param sl The value to scale, see [metagenomeSeq::MRcounts()]
#' @keywords internal
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importFrom metagenomeSeq newMRexperiment cumNorm cumNormStatFast MRcounts
#' @seealso [metagenomeSeq::cumNorm()], [metagenomeSeq::MRcounts()]
norm_css <- function(object,
  p = metagenomeSeq::cumNormStatFast,
  log = FALSE,
  sl = 1000) {
  if (inherits(object, "phyloseq")) {
    object_mgs <- phyloseq_to_metagenomeSeq(object)
  } else if (inherits(object, "otu_table")) {
    # keep in accordance with the phyloseq::phyloseq_to_metagenomeSeq
    count <- round(as(object, "matrix"), digits = 0)
    object_mgs <- newMRexperiment(counts = count)
  }

  p <- ifelse(is.function(p), p(object_mgs), p)
  mgs_normed <- cumNorm(object_mgs, p = p)
  feature_normed <- MRcounts(mgs_normed, norm = TRUE, log = log, sl = sl)

  otu_table(object) <- otu_table(
    feature_normed,
    taxa_are_rows = taxa_are_rows(object)
  )

  object
}

#' Relative log expression (RLE) normalization
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param logcfunc a function to compute a location for a sample. By default,
#'   the median is used.
#' @param geo_means default `NULL`, which means the geometric means of the
#'   counts are used. A vector of geometric means from another count matrix can
#'   be provided for a "frozen" size factor calculation.
#' @param control_genes default `NULL`, which means all genes are used for size
#'   factor estimation, numeric or logical index vector specifying the genes
#'   used for size factor estimation (e.g. housekeeping or spike-in genes).
#' @seealso [DESeq2::estimateSizeFactorsForMatrix()]
#' @keywords internal
norm_rle <- function(object,
  locfunc = stats::median,
  geo_means = NULL,
  control_genes = NULL) {
  otu <- as(otu_table(object), "matrix")

  geo_means <- ifelse(is.null(geo_means), substitute(), geo_means)
  control_genes <- ifelse(is.null(control_genes), substitute(), control_genes)

  sf <- DESeq2::estimateSizeFactorsForMatrix(
    otu,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes
  )
  otu <- sweep(otu, 2, sf, FUN = "/")
  otu_table(object) <- otu_table(otu, taxa_are_rows = taxa_are_rows(object))

  object
}

#' TMM (trimmed mean of m-values) normalization
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param ref_column column to use as reference
#' @param logratio_trim amount of trim to use on log-ratios
#' @param sum_trim amount of trim to use on the combined absolute levels
#'   ("A" values)
#' @param do_weighting whether to compute the weights or not
#' @param Acutoff cutoff on "A" values to use before trimming
#' @seealso [edgeR::calcNormFactors()]
#' @keywords internal
#' @references https://github.com/biobakery/Maaslin2/blob/master/R/utility_scripts.R
norm_tmm <- function(object,
  ref_column = NULL,
  logratio_trim = 0.3,
  sum_trim = 0.05,
  do_weighting = TRUE,
  Acutoff = -1e10) {
  otu <- as(otu_table(object), "matrix")
  nf <- edgeR::calcNormFactors(
    otu,
    method = "TMM",
    refcolumn = ref_column,
    logratioTrim = logratio_trim,
    sumTrim = sum_trim,
    doWeighting = do_weighting,
    Acutoff = Acutoff
  )

  ef_nf <- colSums(otu) * nf

  # use the mean of the effective library size as a reference library size
  ref_nf <- mean(ef_nf)

  otu_norm <- sweep(otu, MARGIN = 2, ef_nf, "/") * ref_nf

  otu_table(object) <- otu_table(otu_norm, taxa_are_rows = taxa_are_rows(object))

  object
}

#' CLR (centered log-ratio) normalization
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @keywords internal
norm_clr <- function(object) {
  otu <- as(otu_table(object), "matrix")

  if (any(otu == 0)) {
    otu <- otu + 1
  }

  otu_norm <- apply(otu, 2, function(x) {log(x) - mean(log(x))})

  otu_table(object) <- otu_table(otu_norm, taxa_are_rows = taxa_are_rows(object))

  object

}
