# In the vignette of DESeq2:
# The values in the matrix should be un-normalized counts or estimated counts
# of sequencing reads (for single-end RNA-seq) or fragments (for paired-end
# RNA-seq). The RNA-seq workflow describes multiple techniques for preparing
# such count matrices. It is important to provide count matrices as input for
# DESeq2’s statistical model (Love, Huber, and Anders 2014) to hold, as only
# the count values allow assessing the measurement precision correctly. The
# DESeq2 model internally corrects for library size, so transformed or
# normalized values such as counts scaled by library size should not be used
# as input.
#
# DESeq2 contrast: https://github.com/tavareshugo/tutorial_DESeq2_contrasts
#
# reference source code:
# from biocore/qiime/blob/master/qiime/support_files/R/DESeq2_nbinom.r
# https://github.com/hbctraining/DGE_workshop/blob/master/schedule/1.5-day.md
#
#
## p value and logFC from LRT
# From https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# https://support.bioconductor.org/p/133804/#133856
# By default the Wald test is used to generate the results table, but DESeq2
# also offers the LRT which is used to identify any genes that show change in
# expression across the different levels. The LRT is comparing the full model
# to the reduced model to identify significant genes. The p-values are
# determined solely by the difference in deviance between the ‘full’ and
# "reduced" model formula (not log2 fold changes).
#
# The log2 fold change LRT results  is calculated using Wald (two groups
# comparison).
#
#
#
#' Perform DESeq differential analysis
#'
#' Differential expression analysis based on the Negative Binomial distribution
#' using **DESeq2**.
#'
#' @param ps  a [`phyloseq::phyloseq-class`] object.
#' @param group  character, the variable to set the group, must be one of
#'   the var of the sample metadata.
#' @param confounders character vector, the confounding variables to be adjusted.
#'   default `character(0)`, indicating no confounding variable.
#' @param contrast this parameter only used for two groups comparison while
#'   there are multiple groups. For more please see the following details.
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of
#'   `phyloseq::rank_names(phyloseq)`, or "all" means to summarize the taxa by
#'   the top taxa ranks (`summarize_taxa(ps, level = rank_names(ps)[1])`), or
#'   "none" means perform differential analysis on the original taxa
#'   (`taxa_names(phyloseq)`, e.g., OTU or ASV).
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation
#'     (default).
#'   * "log10", the transformation is `log10(object)`, and if the data contains
#'     zeros the transformation is `log10(1 + object)`.
#'   * "log10p", the transformation is `log10(1 + object)`.
#' @param norm the methods used to normalize the microbial abundance data. See
#'   [`normalize()`] for more details.
#'   Options include:
#'   * "none": do not normalize.
#'   * "rarefy": random subsampling counts to the smallest library size in the
#'     data set.
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
#' @param norm_para arguments passed to specific normalization methods. Most
#'   users will not need to pass any additional arguments here.
#' @param fitType,sfType,betaPrior,modelMatrixType,useT,minmu these seven
#'   parameters are inherited form [`DESeq2::DESeq()`].
#'   - `fitType`, either "parametric", "local", "mean", or "glmGamPoi" for the
#'     type of fitting of dispersions to the mean intensity.
#'   - `sfType`, either "ratio", "poscounts", or "iterate" for the type of size
#'     factor estimation. We recommend to use "poscounts".
#'   - `betaPrior`, whether or not to put a zero-mean normal prior on the
#'     non-intercept coefficients.
#'   - `modelMatrixType`, either "standard" or "expanded", which describe how
#'     the model matrix,
#'   - `useT`, logical, where Wald statistics are assumed to follow a standard
#'     Normal.
#'   - `minmu`, lower bound on the estimated count for fitting gene-wise
#'     dispersion.
#'
#'   For more details, see [`DESeq2::DESeq()`].  Most users will not need to
#'   set this arguments (just use the defaults).
#'
#' @param p_adjust method for multiple test correction, default `none`, for
#'   more details see [stats::p.adjust].
#' @param pvalue_cutoff pvalue_cutoff numeric, p value cutoff, default 0.05.
#' @param ... extra parameters passed to [`DESeq2::DESeq()`].
#'
#' @details
#' **Note**: DESeq2 requires the input is raw counts (un-normalized counts), as
#' only the counts values allow assessing the measurement precision correctly.
#' For more details see the vignette of DESeq2 (`vignette("DESeq2")`).
#'
#' Thus, this function only supports "none", "rarefy", "RLE", "CSS", and
#' "TMM" normalization methods. We strongly recommend using the "RLE" method
#' (default normalization method in the DESeq2 package). The other
#' normalization methods are used for expert users and comparisons among
#' different normalization methods.
#'
#' For two groups comparison, this function utilizes the Wald test (defined by
#' [`DESeq2::nbinomWaldTest()`]) for hypothesis testing. A Wald test statistic
#' is computed along with a probability (p-value) that a test statistic at least
#' as extreme as the observed value were selected at random. `contrasts` are
#' used to specify which two groups to compare. The order of the names
#' determines the direction of fold change that is reported.
#'
#' Likelihood ratio test (LRT) is used to identify the genes that significantly
#' changed across all the different levels for multiple groups comparisons. The
#' LRT identified the significant features by comparing the full model to the
#' reduced model. It is testing whether a feature removed in the reduced
#' model explains a significant variation in the data.
#'
#' `contrast` must be a two length character or `NULL` (default). It is only
#' required to set manually for two groups comparison when there are multiple
#' groups. The order determines the direction of comparison, the first element
#' is used to specify the reference group (control). This means that, the first
#' element is the denominator for the fold change, and the second element is
#' used as baseline (numerator for fold change). Otherwise, users do required
#' to concern this parameter (set as default `NULL`), and if there are
#' two groups, the first level of groups will set as the reference group; if
#' there are multiple groups, it will perform an ANOVA-like testing to find
#' markers which difference in any of the groups.
#'
#' @export
#' @return a [`microbiomeMarker-class`] object.
#' @seealso [`DESeq2::results()`],[`DESeq2::DESeq()`]
#' @importFrom stats formula coef
#' @importFrom DESeq2 dispersions<-
#' @importMethodsFrom S4Vectors mcols
#' @importMethodsFrom BiocGenerics sizeFactors<- counts
#' @references
#' Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation
#' of fold change and dispersion for RNA-seq data with DESeq2." Genome
#' biology 15.12 (2014): 1-21.
#' @examples
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")) %>%
#'     phyloseq::subset_taxa(Phylum %in% c("Firmicutes"))
#' run_deseq2(ps, group = "Enterotype")
run_deseq2 <- function(ps,
    group,
    confounders = character(0),
    contrast = NULL,
    taxa_rank = "all",
    norm = "RLE",
    norm_para = list(),
    transform = c("identity", "log10", "log10p"),
    # test = c("Wald", "LRT"),
    fitType = c("parametric", "local", "mean", "glmGamPoi"),
    sfType = "poscounts",
    betaPrior = FALSE,
    modelMatrixType,
    useT = FALSE,
    minmu = ifelse(fitType == "glmGamPoi", 1e-06, 0.5),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    ...) {
    ps <- check_rank_names(ps) %>% 
        check_taxa_rank( taxa_rank)
    
    norm_methods <- c("none", "rarefy", "RLE", "CSS", "TMM")
    if (!norm %in% norm_methods) {
        stop(
            "`norm` must be one of 'none', 'rarefy', 'RLE', 'CSS', or 'TMM'",
            call. = FALSE
        )
    }
    
    if (length(confounders)) {
        confounders <- check_confounder(ps, group, confounders)
    }

    # groups
    meta <- sample_data(ps)
    groups <- meta[[group]]
    groups <- make.names(groups)
    if (!is.null(contrast)) {
        contrast <- make.names(contrast)
    }
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    groups <- set_lvl(groups, contrast)
    sample_data(ps)[[group]] <- groups
    lvl <- levels(groups)
    n_lvl <- length(lvl)

    if (n_lvl < 2) {
        stop("Differential analysis requires at least two groups.")
    }

    # contrast, test method, name of effect size
    if (n_lvl == 2) { # two groups
        if (!is.null(contrast)) {
            warning(
                "`contrast` is ignored, you do not need to set it",
                call. = FALSE
            )
        }
        contrast_new <- c(group, lvl[2], lvl[1])
    } else {
        if (!is.null(contrast)) {
            if (!is.character(contrast) || length(contrast) != 2) {
                stop("`contrast` must be a two length character", call. = FALSE)
            }
            idx <- match(contrast, lvl, nomatch = 0L)
            if (!all(idx)) {
                stop(
                    "all elements of `contrast` must be contained in `groups`",
                    call. = FALSE
                )
            }
            contrast_new <- c(group, contrast[2], contrast[1])
        }
    }
    test <- ifelse(n_lvl > 2 && is.null(contrast), "LRT", "Wald")
    ef_name <- ifelse(test == "Wald", "logFC", "F")

    fitType <- match.arg(fitType, c("parametric", "local", "mean", "glmGamPoi"))
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )

    if (!sfType %in% c("ratio", "poscounts", "iterate")) {
        stop("`sfType` muste be one of poscounts, ratio, or iterate")
    }

    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # prenormalize the data
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    
    if (!length(confounders)) {
        dsg <- formula(paste("~", group))
    } else {
        dsg <- formula(paste(
            "~", 
            paste(c(confounders, group), collapse = " + ")
        ))
    }
    
    dds_summarized <- phyloseq2DESeq2(
        ps_summarized,
        design = dsg
    )
    nf <- get_norm_factors(ps_normed)
    if (!is.null(nf)) {
        sizeFactors(dds_summarized) <- nf
    }


    # error: all gene-wise dispersion estimates are within 2 orders of magnitude
    # from the minimum value, which indicates that the count are not
    # overdispersed
    #
    # If dispersion values are less than 1e-6  (minimal value is 1e-8),
    # it would be problematic to fit a dispersion trend in DESeq2.
    # The reason for a minimal value, is that for a given row of the count 
    # matrix, the maximum likelihood estimate can tend to 0 (and so we have a 
    # rule to stop after 1e-8)
    # https://support.bioconductor.org/p/63845/
    # https://support.bioconductor.org/p/122757/
    # from biocore/qiime/blob/master/qiime/support_files/R/DESeq2_nbinom.r

    # LRT is used to analyze all levels of a factor at once, and the
    # The p values are determined solely by the difference in deviance between
    # the "full" and "reduced" model formula (not log2 fold changes). Only Wast
    # method was used for pair-wise comparison. Thus, for pair-wise comparison,
    # we use Wald test. Moreover, you can set the argument `test` in `results()`
    # when extract the results from LRT, and it is equivalent to Wast test.
    #
    # However, even though there are fold changes  present in the results of
    # LRT, they are not directly associated with the actual hypothesis test (
    # actually determined by the arguments name or contrast).
    if (test == "Wald") { # two groups comparison
        res_deseq <- try(
            DESeq2::DESeq(
                dds_summarized,
                test = test,
                fitType = fitType,
                sfType = sfType,
                quiet = TRUE,
                betaPrior = betaPrior,
                modelMatrixType = modelMatrixType,
                useT = useT,
                minmu = minmu,
                ...
            ),
            silent = TRUE
        )

        if (inherits(res_deseq, "try-error") && fitType != "local") {
            warning("data is not overdispered, try `fitType = 'local'`")
            res_deseq <- try(
                DESeq2::DESeq(
                    dds_summarized,
                    test = test,
                    fitType = "local",
                    sfType = sfType,
                    quiet = TRUE,
                    betaPrior = betaPrior,
                    modelMatrixType = modelMatrixType,
                    useT = useT,
                    minmu = minmu,
                    ...
                ),
                silent = TRUE
            )
        }

        if (inherits(res_deseq, "try-error") && fitType != "mean") {
            warning("data is not overdispered, try `fitType = 'mean'`")
            res_deseq <- try(
                DESeq2::DESeq(
                    dds_summarized,
                    test = test,
                    fitType = "mean",
                    sfType = sfType,
                    quiet = TRUE,
                    betaPrior = betaPrior,
                    modelMatrixType = modelMatrixType,
                    useT = useT,
                    minmu = minmu,
                    ...
                ),
                silent = TRUE
            )
        }
        if (inherits(res_deseq, "try-error")) {
            warning(
                "data is not overdispered, use gene-wise estimates ",
                "as final estimates"
            )
            dds_summarized <- DESeq2::estimateDispersionsGeneEst(dds_summarized)
            DESeq2::dispersions(dds_summarized) <- 
                mcols(dds_summarized)$dispGeneEst

            dds_summarized <- DESeq2::nbinomWaldTest(
                dds_summarized,
                betaPrior = betaPrior,
                quiet = TRUE,
                modelMatrixType = modelMatrixType,
                useT = useT,
                minmu = minmu
            )
        } else {
            dds_summarized <- res_deseq
        }

        res <- DESeq2::results(
            object = dds_summarized,
            contrast = contrast_new,
            pAdjustMethod = p_adjust
        )
        # rename log2FoldChange to logFC, use base R rather than dplyr::rename
        names(res)[names(res) == "log2FoldChange"] <- "logFC"
    } else {
        dds_summarized <- DESeq2::DESeq(
            dds_summarized,
            test = test,
            fitType = fitType,
            sfType = sfType,
            quiet = TRUE,
            minmu = minmu,
            reduced = ~1,
            ...
        )
        res <- DESeq2::results(
            object = dds_summarized,
            pAdjustMethod = p_adjust
        )
    }

    # Why p value is NA?
    # By default, independent filtering is performed to select a set of genes
    # for multiple test correction which maximizes the number of adjusted p
    # values less than a given critical value alpha (by default 0.1).
    # The adjusted p-values for the genes which do not pass the filter threshold
    # are set to NA.
    # By default, results assigns a p-value of NA to genes containing count
    # outliers, as identified using Cook's distance.

    # normalized counts
    counts_normalized <- DESeq2::counts(dds_summarized, normalized = TRUE)

    # one way anova f statistic for LRT
    if (test == "LRT") {
        temp_count <- data.frame(t(counts_normalized))
        f_stat <- vapply(
            temp_count,
            calc_ef_md_f,
            FUN.VALUE = 0.0,
            group = groups
        )
        res[["F"]] <- f_stat
    }

    res <- data.frame(res)
    # enrich group
    if (test == "Wald") {
        enrich_group <- ifelse(res$logFC > 0, contrast_new[2], contrast_new[3])
    } else {
        cf <- coef(dds_summarized)
        
        # extract coef of interested var
        target_idx <- grepl(group, colnames(cf))
        cf <- cf[, target_idx]
        # the first coef is intercept, bind the coef of the reference group as 0
        # (the first column)
        cf <- cbind(0, cf)
        enrich_idx <- apply(
            cf, 1,
            function(x) ifelse(any(is.na(x)), NA, which.max(x))
        )
        enrich_group <- lvl[enrich_idx]
        enrich_group <- enrich_group[match(row.names(res), row.names(cf))]
    }
    res$enrich_group <- enrich_group
    # order according to padj
    res_ordered <- res[order(res$padj), ]
    # filter sig feature
    padj <- res_ordered$padj
    res_ordered <- cbind(feature = row.names(res_ordered), res_ordered)
    # rownames in the form of marker*
    row.names(res_ordered) <- paste0("marker", seq_len(nrow(res_ordered)))

    # reorder columns: feature, enrich_group, other columns
    other_col <- setdiff(names(res_ordered), c("feature", "enrich_group"))
    res_ordered <- res_ordered[, c("feature", "enrich_group", other_col)]
    row.names(res_ordered) <- paste0("marker", seq_len(nrow(res_ordered)))

    # only keep five variables: feature, enrich_group, effect_size (logFC),
    # pvalue, and padj
    keep_var <- c("feature", "enrich_group", ef_name, "pvalue", "padj")
    res_ordered <- res_ordered[keep_var]
    names(res_ordered)[3] <- paste0("ef_", ef_name)
    sig_res <-  res_ordered[!is.na(padj) & padj < pvalue_cutoff, ]
    marker <- return_marker(sig_res, res_ordered)

    marker <- microbiomeMarker(
        marker_table = marker,
        # if no pre-calculated size factors, DESeq2 will calculate the 
        # size factors internally, so norm method shoule be RLE
        norm_method = ifelse(is.null(nf), "RLE", get_norm_method(norm)),
        diff_method = paste0("DESeq2: ", test),
        sam_data = sample_data(ps_normed),
        tax_table = tax_table(ps_summarized),
        otu_table = otu_table(counts_normalized, taxa_are_rows = TRUE)
    )

    marker
}

#' Convert `phyloseq-class` object to `DESeqDataSet-class` object
#'
#' This function convert [phyloseq::phyloseq-class`] to
#' [`DESeq2::DESeqDataSet-class`], which can then be tested using
#' [`DESeq2::DESeq()`].
#'
#' @param ps the [phyloseq::phyloseq-class`] object to convert, which must have
#'   a [`phyloseq::sample_data()`] component.
#' @param design a `formula` or `matrix`, the formula expresses how the counts
#'   for each gene depend on the variables in colData. Many R formula are valid,
#'   including designs with multiple variables, e.g., `~ group + condition`.
#'   This argument is passed to [`DESeq2::DESeqDataSetFromMatrix()`].
#' @param ... additional arguments passed to
#'   [`DESeq2::DESeqDataSetFromMatrix()`], Most users will not need to pass any
#'   additional arguments here.
#' @export
#' @return a [`DESeq2::DESeqDataSet-class`] object.
#' @seealso [`DESeq2::DESeqDataSetFromMatrix()`],[`DESeq2::DESeq()`]
#' @examples
#' data(caporaso)
#' phyloseq2DESeq2(caporaso, ~SampleType)
phyloseq2DESeq2 <- function(ps, design, ...) {
    stopifnot(inherits(ps, "phyloseq"))
    ps <- keep_taxa_in_rows(ps)

    # sample data
    samp <- sample_data(ps, errorIfNULL = FALSE)
    if (is.null(samp)) {
        stop(
            "`sample_data` of `ps` is required,",
            " for specifying experimental design.",
            call. = FALSE
        )
    }
    # count data
    ct <- as(otu_table(ps), "matrix")

    # deseq2 requires raw counts, means the counts must be integer
    dds <- tryCatch(
        DESeq2::DESeqDataSetFromMatrix(
            countData = ct,
            colData = data.frame(samp),
            design = design,
            ...
        ),
        error = function(e) e
    )
    if (inherits(dds, "error") &&
        conditionMessage(dds) == "some values in assay are not integers") {
        warning(
            "Some counts are non-integers, they are rounded to integers.\n",
            "Raw count is recommended for reliable results for deseq2 method.",
            call. = FALSE
        )
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = round(ct),
            colData = data.frame(samp),
            design = design,
            ...
        )
    }


    dds
}

# Modified from `DESeq2::estimateFactorsForMatrix()` directly
# for `estimateSizeFactors`:
# `sizeFactors(estimateSizeFactors(dds, type = "poscounts"))` is identical to
# `sizeFactors(estimateSizeFactors(dds, geoMeans = geoMeans))`
#
# The original function of `DESeq2::estimateFactorsForMatrix()` does not
# stabilize size factors to have geometric mean of 1 while `type = "poscounts"`.
# This modified function is to make
# `estimateSizeFactorsForMatrix(counts(diagdds2),geoMeans = geoMeans)` is equal
# to `estimateSizeFactorsForMatrix(counts(diagdds2), type = "poscounts")` by
# stabilize size factors if `type = "poscounts"`.
estimateSizeFactorsForMatrix <- function(counts,
    locfunc = stats::median,
    geoMeans,
    controlGenes,
    type = c("ratio", "poscounts")) {
    type <- match.arg(type, c("ratio", "poscounts"))
    if (missing(geoMeans)) {
        incomingGeoMeans <- FALSE
        if (type == "ratio") {
            loggeomeans <- rowMeans(log(counts))
        } else if (type == "poscounts") {
            lc <- log(counts)
            lc[!is.finite(lc)] <- 0
            loggeomeans <- rowMeans(lc)
            allZero <- rowSums(counts) == 0
            loggeomeans[allZero] <- -Inf
        }
    } else {
        incomingGeoMeans <- TRUE
        if (length(geoMeans) != nrow(counts)) {
            stop("geoMeans should be as long as the number of rows of counts")
        }
        loggeomeans <- log(geoMeans)
    }
    if (all(is.infinite(loggeomeans))) {
        stop(
            "every gene contains at least one zero ",
            "cannot compute log geometric means",
            call. = FALSE
        )
    }
    sf <- if (missing(controlGenes)) {
        apply(counts, 2, function(cnts) {
            exp(locfunc((log(cnts) - loggeomeans)[
                is.finite(loggeomeans) & cnts > 0]))
        })
    } else {
        if (!(is.numeric(controlGenes) | is.logical(controlGenes))) {
            stop("controlGenes should be either a numeric or logical vector")
        }
        loggeomeansSub <- loggeomeans[controlGenes]
        apply(
            counts[controlGenes, , drop = FALSE], 2,
            function(cnts) {
                idx <- is.finite(loggeomeansSub) & cnts > 0
                exp(locfunc((log(cnts) - loggeomeansSub)[idx]))
            }
        )
    }
    if (incomingGeoMeans | type == "poscounts") {
        # stabilize size factors to have geometric mean of 1
        sf <- sf / exp(mean(log(sf)))
    }
    sf
}
