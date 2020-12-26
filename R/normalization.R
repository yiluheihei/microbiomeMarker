#' Normalize the microbial abundance data
#'
#' @param object a matrix, data.frame, [phyloseq::phyloseq-class] or
#'   [phyloseq::otu_table-class] object
#' @param method the methods used to normalize the microbial abundance data.
#'   Options includes:
#'   * a integer, e.g. 1e6, indicating pre-sample normalization of the sum of
#'     the values to 1e6.
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
#' @importFrom phyloseq sample_data<-
#' @aliases normalize,phyloseq-method
#' @rdname normalize-methods
setMethod("normalize", "phyloseq",
  function(object,
           method = "TSS",
           ...) {
    otu <- otu_table(object)
    otu_normed <- normalize(otu, method = method, ...)

    # extract norm_factor attributes and prepend to the sample_data
    attrs <- attributes(otu_normed)
    nf_idx <- grepl("norm_factor", names(attrs))
    if (any(nf_idx)) {
      nf_name <- names(attrs)[nf_idx]
      nf <- attrs[[nf_name]]
      sample_data(object) <- cbind(
        sample_data(object),
        nf_name = attr(otu_normed, nf_name, exact = TRUE)
      )
    }

    otu_table(object) <- otu_table(
      otu_normed,
      taxa_are_rows = taxa_are_rows(object)
    )

    object
  }
)

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,otu_table-method normalize
#' @rdname normalize-methods
setMethod("normalize", "otu_table",
  function(object,
           method = "TSS",
           ...) {
    if (method %in% c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR")) {
      object_normed <- switch (method,
        none = object,
        rarefy = norm_rarefy(object, ...),
        TSS = norm_tss(object),
        TMM = norm_tmm(object, ...),
        RLE = norm_rle(object, ...),
        CSS = norm_css(object, ...),
        CLR = norm_clr(object)
      )
    } else if (is.numeric(method)) {
      object_normed <- norm_value(object, normalization = method)
    } else {
      stop("`method` must be one of none, rarefy, TSS, TMM, RLE, CSS, CLR, or an integer")
    }

    object_normed
  })

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,data.frame-method normalize
#' @rdname normalize-methods
setMethod("normalize", "data.frame",
  function(object,
           method = "TSS",
           ...) {
    otu <- otu_table(object, taxa_are_rows = TRUE)
    otu_norm <- normalize(otu, method, ...)

    as.data.frame(otu_norm)
  }
)

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,matrix-method normalize
#' @rdname normalize-methods
setMethod("normalize", "matrix",
  function(object,
           method = "TSS",
           ...) {
    otu <- as.data.frame(object)
    otu_norm <- normalize(otu, method, ...)

    as.matrix(otu_norm)
  }
)

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
#' @param log logical, whether or not to log2 transfrom scale,
#'   see [metagenomeSeq::MRcounts()].
#' @param sl The value to scale, see [metagenomeSeq::MRcounts()]
#' @keywords internal
#' @importFrom phyloseq sample_data<-
#' @importFrom metagenomeSeq newMRexperiment cumNorm cumNormStatFast MRcounts
#' @seealso [metagenomeSeq::cumNorm()], [metagenomeSeq::MRcounts()]
norm_css <- function(object,
                     log = FALSE,
                     sl = 1000) {
  if (inherits(object, "phyloseq")) {
    object_mgs <- phyloseq2metagenomeSeq(object)
  } else if (inherits(object, "otu_table")) {
    # keep in accordance with the phyloseq::phyloseq_to_metagenomeSeq
    # count <- round(as(object, "matrix"), digits = 0)
    object_mgs <- otu_table2metagenomeSeq(object)
  }

  # cumNormStatFast requires counts of all samples at least have two
  # non zero features
  count <- as(otu_table(object), "matrix")
  if (sum(colSums(count > 0) > 1) < ncol(count)) {
    p <- suppressMessages(metagenomeSeq::cumNormStat(object_mgs))
  }
  else {
    p <- suppressMessages(metagenomeSeq::cumNormStatFast(object_mgs))
  }
  object_mgs <- metagenomeSeq::cumNorm(object_mgs, p = p)

  count_normed <- MRcounts(object_mgs, norm = TRUE, log = log, sl = sl)
  otu_table(object) <- otu_table(
    count_normed,
    taxa_are_rows = taxa_are_rows(object)
  )

  # append normFactor to sample_data for model fitting if object is phyloseq
  # set the metagenomeSeq_norm_factor attributes if object is otu_table
  nf <- metagenomeSeq::normFactors(object_mgs)
  if (inherits(object, "phyloseq")) {
    sample_data(object) <- cbind(
      sample_data(object),
     metagenomeSeq_norm_factor = nf
    )
  } else {
    attr(object, "metagenomeSeq_norm_factor") <- nf
  }

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

#' Normalize the sum of values of each sample to a given value
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#' @param normalization  set the normalization value, return the feature itself
#'   if not supported.
#' @keywords internal
#' @importFrom phyloseq transform_sample_counts
norm_value <- function(object, normalization) {
  if (missing(normalization)) {
    return(object)
  }

  otu <- as(otu_table(object), "matrix") %>%
    as.data.frame()

  # whether the object is summarized
  hie <- check_tax_summarize(object)
  if (hie) {
    features <- row.names(otu)
    features_split <- strsplit(features, "|", fixed = TRUE)
    single_indx <- which(lengths(features_split) < 2)

    ## keep the counts of a sample identical with `normalization`
    ## if we norm the counts in two steps:
    ## 1. calculate scale size: norm_coef = normalization/lib_size
    ## 2. multiple the scale size value * norm_coef
    ## the counts of a sample colSums(otu) may not equal to the argument normalization
    ## e.g. normalization = 1e6, colSums(otu) = 999999
    ## Finally, the kruskal test may be inaccurate,
    ## e.g. https://github.com/yiluheihei/microbiomeMarker/issues/13
    ps_normed <- transform_sample_counts(
      object,
      function(x) x * normalization/ sum(x[single_indx])
    )

    # lib_size <- purrr::map_dbl(otu, ~ sum(.x[single_indx]))

    # `sum(abd)` must be greaer than 0 since the missing level is added
    # if (sum(abd) == 0) {
    #   abd <- purrr::map_dbl(feature, sum)
    # }
  } else {
    ps_normed <- transform_sample_counts(
      object,
      function(x) x * normalization / sum(x)
    )
  }

  otu_normed <- data.frame(otu_table(ps_normed))
  otu_normed <- purrr::map_df(otu_normed,
    function(x) {
      if (mean(x) && stats::sd(x)/mean(x) < 1e-10) {
        return(round(x*1e6)/1e6)
      } else {
        return(x)
      }
    }
  )

  # normed_coef <- normalization/lib_size

  # otu_normed <- purrr::map2_df(
  #   otu, normed_coef,
  #   function(x, y) {
  #     res <- x * y
  #     if (mean(res) && stats::sd(res)/mean(res) < 1e-10) {
  #       res <- round(res * 1e6)/1e6
  #     }
  #     res
  #   }
  # )

  otu_normed <- as.data.frame(otu_normed)
  row.names(otu_normed) <- row.names(otu)
  colnames(otu_normed) <- colnames(otu)
  otu_table(object) <- otu_table(otu_normed, taxa_are_rows = TRUE)

  object

}


#' normalize the summarized feature
#' @param feature otu table or data.frame
#' @param normalization  set the normalization value
#' @noRd
normalize_feature <- function(feature, normalization) {
  if (inherits(feature, "otu_table")) {
    if (!taxa_are_rows(feature)) {
      feature <- t(feature)
    }
    feature <- feature@.Data %>% data.frame()
  }
  if (is.null(normalization)) {
    return(feature)
  }

  feature_split <- strsplit(row.names(feature), "\\|")
  hie <- ifelse(any(lengths(feature_split) > 1), TRUE, FALSE)
  if (hie) {
    single_indx <- which(lengths(feature_split) < 2)
    abd <- purrr::map_dbl(feature, ~ sum(.x[single_indx]))

    # `sum(abd)` must be greaer than 0 since the missing level is added
    # if (sum(abd) == 0) {
    #   abd <- purrr::map_dbl(feature, sum)
    # }
  } else {
    abd <- purrr::map_dbl(feature, sum)
  }
  normed_coef <- normalization/abd
  normed_feature <- purrr::map2_df(
    feature, normed_coef,
    function(x, y) {
      res <- x * y
      if (mean(res) && stats::sd(res)/mean(res) < 1e-10) {
        res <- round(res * 1e6)/1e6
      }
      res
    }
  )

  # for row names setting, phyloseq requires otu_table and tax_table has the
  # same taxa
  normed_feature <- as.data.frame(normed_feature)
  row.names(normed_feature) <- row.names(feature)

  otu_table(normed_feature, taxa_are_rows = TRUE)
}
