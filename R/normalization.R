#' Normalize the microbial abundance data
#'
#' It is critical to normalize the feature table to eliminate any bias due to
#' differences in the sampling sequencing depth.This function implements six
#' widely-used normalization methods for microbial compositional data.
#'
#' @param object a matrix, data.frame, [phyloseq::phyloseq-class] or
#'   [phyloseq::otu_table-class] object.
#' @param method the methods used to normalize the microbial abundance data.
#'   Options includes:
#'   * "none": do not normalize.
#'   * "rarefy": random subsampling counts to the smallest library size in the
#'     data set.
#'   * "TSS": total sum scaling, also referred to as "relative abundance", the
#'     abundances were normalized by dividing the corresponding sample library
#'     size.
#'   * "TMM": trimmed mean of m-values. First, a sample is chosen as reference.
#'     The scaling factor is then derived using a weighted trimmed mean over
#'     the differences of the log-transformed gene-count fold-change between
#'     the sample and the reference.
#'   * "RLE", relative log expression, RLE uses a pseudo-reference calculated
#'     using the geometric mean of the gene-specific abundances over all
#'     samples. The scaling factors are then calculated as the median of the
#'     gene counts ratios between the samples and the reference.
#'   * "CSS": cumulative sum scaling, calculates scaling factors as the
#'     cumulative sum of gene abundances up to a data-derived threshold.
#'   * "CLR": centered log-ratio normalization.
#'   * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' @param ... other arguments passed to the corresponding normalization
#'   methods.
#' @seealso [edgeR::calcNormFactors()],[DESeq2::estimateSizeFactorsForMatrix()],
#' [metagenomeSeq::cumNorm()]
#' @importMethodsFrom BiocGenerics normalize
#' @importFrom phyloseq sample_data<-
#' @exportMethod normalize
#' @aliases normalize,phyloseq-method
#' @rdname normalize-methods
#' @return the same class with `object`.
#' @examples
#' data(caporaso)
#' normalize(caporaso, "TSS")
setMethod(
    "normalize", "phyloseq",
    function(object,
    method = "TSS",
    ...) {
        otu <- otu_table(object)
        otu_normed <- normalize(otu, method = method, ...)

        # extract norm_factor attributes and prepend to the sample_data
        nf <- attr(otu_normed, "norm_factor")
        if (!is.null(nf)) {
            sample_data(object) <- cbind(
                sample_data(object),
                norm_factor = nf
            )
        }

        # please note otu_table<- function will drop the norm_factor attribute
        otu_table(object) <- otu_normed

        object
    }
)

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,otu_table-method normalize
#' @rdname normalize-methods
setMethod(
    "normalize", "otu_table",
    function(object,
    method = "TSS",
    ...) {
        methods <- c("none", "rarefy", "TSS", "TMM", "RLE", "CSS", "CLR", "CPM")
        if (method %in% methods) {
            object_normed <- switch(method,
                none = object,
                rarefy = norm_rarefy(object, ...),
                TSS = norm_tss(object),
                TMM = norm_tmm(object, ...),
                RLE = norm_rle(object, ...),
                CSS = norm_css(object, ...),
                CLR = norm_clr(object),
                CPM = norm_cpm(object)
            )
        } else {
            stop(
                "`method` must be one of none, rarefy, TSS,",
                " TMM, RLE, CSS, CLR, or CPM",
                call. = FALSE
            )
        }

        object_normed
    }
)

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,data.frame-method normalize
#' @rdname normalize-methods
setMethod(
    "normalize", "data.frame",
    function(object,
    method = "TSS",
    ...) {
        otu <- otu_table(object, taxa_are_rows = TRUE)
        otu_norm <- normalize(otu, method, ...)
        nf <- attr(otu_norm, "norm_factor")

        res <- as.data.frame(otu_norm)
        if (!is.null(nf)) {
            attr(res, "norm_factor") <- nf
        }

        res
    }
)

#' @importMethodsFrom BiocGenerics normalize
#' @aliases normalize,matrix-method normalize
#' @rdname normalize-methods
setMethod(
    "normalize", "matrix",
    function(object,
    method = "TSS",
    ...) {
        otu <- as.data.frame(object)
        otu_norm <- normalize(otu, method, ...)
        nf <- attr(otu_norm, "norm_factor")

        res <- as.matrix(otu_norm)
        if (!is.null(nf)) {
            attr(res, "norm_factor") <- nf
        }

        res
    }
)

## Four normalization methods do not save the norm factor: value, rarefy,
## clr and tss; where three methods save the norm factor: css, rle, tmm.

#' Normalize feature table by rafefying such that all samples have the same
#' number of total counts (library size).
#'
#' For rarefying, reads in the different samples are randomly removed until
#' the same predefined number has been reached, to assure all samples have the
#' same library size. Rarefying normalization method is the standard in
#' microbial ecology. Please note that the authors of phyloseq do not advocate
#' using this rarefying a normalization procedure, despite its recent
#' popularity
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#' object.
#' @param size,rng_seed,replace,trim_otus,verbose extra arguments passed to
#'   [`phyloseq::rarefy_even_depth()`].
#' @export
#' @rdname normalize-methods
#' @aliases norm_rarefy
#' @importFrom phyloseq rarefy_even_depth sample_sums
#' @seealso [`phyloseq::rarefy_even_depth()`]
norm_rarefy <- function(object,
    size = min(sample_sums(object)),
    rng_seed = FALSE,
    replace = TRUE,
    trim_otus = TRUE,
    verbose = TRUE) {
    object_rarefied <- rarefy_even_depth(
        object,
        sample.size = size,
        rngseed = rng_seed,
        replace = replace,
        trimOTUs = trim_otus,
        verbose = verbose
    )

    # do not save the norm_factor
    # the norm factors can be calculated in the subsequently differential
    # analysis method, e.g. edgeR, DESeq
    object_rarefied
}

#' Total-Sum Scaling (TSS) method
#'
#' TSS simply transforms the feature table into relative abundance by dividing
#' the number of total reads of each sample.
#'
#' @param object object a [phyloseq::phyloseq-class] or
#'   [phyloseq::otu_table-class] object
#' @export
#' @rdname normalize-methods
#' @aliases norm_tss
#' @importFrom phyloseq otu_table<-
norm_tss <- function(object) {
    otu <- otu_table(object)
    size <- colSums(otu)
    otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
    otu_table(object) <- otu_table(
        otu_normed,
        taxa_are_rows = taxa_are_rows(object)
    )

    # do not save the norm_factor, the norm factors are calculated based on the
    # subsequently differential analysis method, e.g. edgeR, DESeq
    object
}

#' Cumulative-Sum Scaling (CSS) method
#'
#' CSS is based on the assumption that the count distributions in each sample
#' are equivalent for low abundant genes up to a certain threshold.  Only the
#' segment of each sample’s count distribution that is relatively invariant
#' across samples is scaled by CSS
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object.
#' @param sl The value to scale.
#' @importFrom phyloseq sample_data<-
#' @importFrom metagenomeSeq newMRexperiment cumNorm cumNormStatFast MRcounts
#' @seealso [metagenomeSeq::calcNormFactors()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_css
norm_css <- function(object, sl = 1000) {
    if (inherits(object, "phyloseq")) {
        object_mgs <- phyloseq2metagenomeSeq(object)
    } else if (inherits(object, "otu_table")) {
        object_mgs <- otu_table2metagenomeSeq(object)
    }

    # cumNormStatFast requires counts of all samples at least have two
    # non zero features. Thus, if there are samples with only one non-zer
    # features, cumNormStat is taken to compute the pth quantile.
    count <- as(otu_table(object), "matrix")
    fun_p <- select_quantile_func(count)
    nf <- metagenomeSeq::calcNormFactors(object_mgs, p = fun_p(object_mgs))
    nf <- unlist(nf) / sl
    object_nf <- set_nf(object, nf)

    object_nf
}

#' Relative log expression (RLE) normalization
#'
#' RLE assumes most features are not differential and uses the relative
#' abundances to calculate the normalization factor.
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param locfunc a function to compute a location for a sample. By default,
#'   the median is used.
#' @param type method for estimation: either "ratio"or "poscounts" (recommend).
#' @param geo_means default `NULL`, which means the geometric means of the
#'   counts are used. A vector of geometric means from another count matrix can
#'   be provided for a "frozen" size factor calculation.
#' @param control_genes default `NULL`, which means all taxa are used for size
#'   factor estimation, numeric or logical index vector specifying the taxa
#'   used for size factor estimation (e.g. core taxa).
#' @seealso [DESeq2::estimateSizeFactorsForMatrix()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_rle
norm_rle <- function(object,
    locfunc = stats::median,
    type = c("poscounts", "ratio"),
    geo_means = NULL,
    control_genes = NULL) {
    stopifnot(class(object) %in% c("phyloseq", "otu_table"))
    type <- match.arg(type, c("poscounts", "ratio"))

    # use substitute() to create missing argument
    geo_means <- ifelse(is.null(geo_means), substitute(), geo_means)
    control_genes <- ifelse(is.null(control_genes), substitute(), control_genes)

    otu <- as(otu_table(object), "matrix")
    nf <- estimateSizeFactorsForMatrix(
        otu,
        locfunc = locfunc,
        geoMeans = geo_means,
        controlGenes = control_genes,
        type = type
    )
    object_nf <- set_nf(object, nf)

    object_nf
}

# https://github.com/biobakery/Maaslin2/blob/master/R/utility_scripts.R
#
#' TMM (trimmed mean of m-values) normalization
#'
#' TMM calculates the normalization factor using a robust statistics based on
#' the assumption that most features are not differential and should, in
#' average, be equal between the samples. The TMM scaling factor is  calculated
#' as the weighted mean of log-ratios between each pair of samples, after
#' excluding the highest count OTUs and OTUs with the largest log-fold change.
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @param ref_column column to use as reference
#' @param logratio_trim amount of trim to use on log-ratios
#' @param sum_trim amount of trim to use on the combined absolute levels
#'   ("A" values)
#' @param do_weighting whether to compute the weights or not
#' @param Acutoff cutoff on "A" values to use before trimming
#' @seealso [edgeR::calcNormFactors()]
#' @export
#' @rdname normalize-methods
#' @aliases norm_tmm
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
    object_nf <- set_nf(object, nf)

    object_nf
}

#' CLR (centered log-ratio) normalization
#'
#' In CLR, the log-ratios are computed relative to the geometric mean of all
#' features.
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#'   object
#' @export
#' @rdname normalize-methods
#' @aliases norm_clr
norm_clr <- function(object) {
    otu <- as(otu_table(object), "matrix")
    otu_norm <- apply(otu, 2, trans_clr)

    otu_table(object) <- otu_table(
        otu_norm,
        taxa_are_rows = taxa_are_rows(object)
    )

    # do not save the norm_factor, the norm factors are calculated based on the
    # subsequently differential analysis method, e.g. edgeR, DESeq
    object
}

# from joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
gm_mean <- function(x, na.rm = TRUE) {
    # The geometric mean, with some error-protection bits.
    exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
}

trans_clr <- function(x, base = exp(1)) {
    x <- log((x / gm_mean(x)), base)
    x[!is.finite(x) | is.na(x)] <- 0.0
    return(x)
}

#' Normalize the sum of values of each sample to million (counts per million)
#'
#' `norm_cpm`: This normalization method is from the original LEfSe algorithm,
#' recommended when very low values are present (as shown in the LEfSe galaxy).
#'
#' @param object a [phyloseq::phyloseq-class] or [phyloseq::otu_table-class]
#' @export
#' @rdname normalize-methods
#' @aliases norm_cpm
#' @importFrom phyloseq transform_sample_counts
norm_cpm <- function(object) {
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
        ## 1. calculate scale size: norm_coef = normalization/lib_size;
        ## 2. multiple the scale size value * norm_coef
        ## the counts of a sample colSums(otu) may not equal to the argument
        ## normalization.
        ## e.g. normalization = 1e6, colSums(otu) = 999999
        ## Finally, the kruskal test may be inaccurate,
        ## e.g. https://github.com/yiluheihei/microbiomeMarker/issues/13
        ps_normed <- transform_sample_counts(
            object,
            function(x) x * 1e+06 / sum(x[single_indx])
        )
    } else {
        ps_normed <- transform_sample_counts(
            object,
            function(x) x * 1e+06 / sum(x)
        )
    }

    otu_normed <- data.frame(otu_table(ps_normed))
    otu_normed <- purrr::map_df(
        otu_normed,
        function(x) {
            if (mean(x) && stats::sd(x) / mean(x) < 1e-10) {
                return(round(x * 1e6) / 1e6)
            } else {
                return(x)
            }
        }
    )

    otu_normed <- as.data.frame(otu_normed)
    row.names(otu_normed) <- row.names(otu)
    colnames(otu_normed) <- colnames(otu)
    otu_table(object) <- otu_table(otu_normed, taxa_are_rows = TRUE)

    # do not save the norm_factor, the norm factors are calculated based on the
    # subsequently differential analysis method, e.g. edgeR, DESeq
    object
}

# set the norm factors of the object, if object is in phyloseq-class, add a
# var `norm_factor` in the sample_data to save the norm factors of each sample;
# if object in otu_table-class, add a attributes `norm_factor` to the object
# to save the norm factors.
#' @importFrom phyloseq sample_data<-
#' @noRd
set_nf <- function(object, nf) {
    # ensure norm factors from sample_data and attributes of otu_table are
    # identical
    names(nf) <- NULL

    if (inherits(object, "phyloseq")) {
        sample_data(object) <- cbind(sample_data(object), norm_factor = nf)
        # to keep accordance with otu_table,
        # we also add the attributes `norm_factor` to otu_table
        ot <- otu_table(object)
        attr(ot, "norm_factor") <- nf
        otu_table(object) <- ot
    } else if (inherits(object, "otu_table")) {
        attr(object, "norm_factor") <- nf
    } else {
        stop("object must be a `phloseq` or `otu_table` object")
    }

    object
}


#' Extract the normalization factors
#'
#' This function will be used to extract the normalization factors. After
#' dividing the observed feature table by normalization factors (eliminate
#' sequencing biases), we will obtain the normalized feature table.
#'
#' @param object a [`phyloseq::phyloseq-class`], [phyloseq::otu_table-class]
#'   object.
#' @return a numeric vector with the length equal to the number of samples, or
#'  `NULL` if the `object` has not been normalized.
#' @noRd
get_norm_factors <- function(object) {
    if (inherits(object, "phyloseq")) {
        nf <- sample_data(object)$norm_factor
    } else {
        nf <- attr(object, "norm_factor")
    }

    nf
}


# Deprecated functions ----------------------------------------------------

# This function is deprecated
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
    } else {
        abd <- purrr::map_dbl(feature, sum)
    }
    normed_coef <- normalization / abd
    normed_feature <- purrr::map2_df(
        feature, normed_coef,
        function(x, y) {
            res <- x * y
            if (mean(res) && stats::sd(res) / mean(res) < 1e-10) {
                res <- round(res * 1e6) / 1e6
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
