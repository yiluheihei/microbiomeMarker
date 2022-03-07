#' Identify biomarkers using supervised leaning (SL) methods
#'
#' Identify biomarkers using logistic regression, random forest, or support
#' vector machine.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param group character, the variable to set the group.
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
#'   * "CPM": pre-sample normalization of the sum of the values to 1e+06.
#' @param norm_para  named `list`. other arguments passed to specific
#'   normalization methods.  Most users will not need to pass any additional
#'   arguments here.
#' @param nfolds the number of splits in CV.
#' @param nrepeats the number of complete sets of folds to compute.
#' @param sampling a single character value describing the type of additional
#'   sampling that is conducted after resampling (usually to resolve class
#'   imbalances). Values are "none", "down", "up", "smote", or "rose". For
#'   more details see [`caret::trainControl()`].
#' @param tune_length an integer denoting the amount of granularity in the
#'   tuning parameter grid. For more details see [`caret::train()`].
#' @param top_n an integer denoting the top `n` features as the biomarker
#'   according the importance score.
#' @param method supervised learning method, options are "LR" (logistic
#'   regression), "RF" (rando forest), or "SVM" (support vector machine).
#' @param ... extra arguments passed to the classification. e.g., `importance`
#'   for `randomForest::randomForest`.
#'
#' @details Only support two groups comparison in the current version. And the
#'  marker was selected based on its importance score. Moreover, The
#'  hyper-parameters are selected automatically by a grid-search based method
#'  in the N-time K-fold cross-validation. Thus, the identified biomarker based
#'  can be biased due to model overfitting for small datasets (e.g., with less
#'  than 100 samples).
#'
#'  The argument `top_n` is used to denote the number of markers based on the
#'  importance score. There is no rule or principle on how to select `top_n`,
#'  however, usually it is very useful to try a different `top_n` and compare
#'  the performance of the marker predictions for the testing data.
#'
#' @return a [microbiomeMarker-class] object.
#'
#' @seealso [`caret::train()`],[`caret::trainControl()`]
#'
#' @author Yang Cao
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' # small example phyloseq object for test
#' ps_small <- phyloseq::subset_taxa(
#'     enterotypes_arumugam,
#'     Phylum %in% c("Firmicutes", "Bacteroidetes")
#' )
#'
#' set.seed(2021)
#' mm <- run_sl(
#'     ps_small,
#'     group = "Gender",
#'     taxa_rank = "Genus",
#'     nfolds = 2,
#'     nrepeats = 1,
#'     top_n = 15,
#'     norm = "TSS",
#'     method = "LR",
#' )
#' mm
run_sl <- function(ps,
    group,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "none",
    norm_para = list(),
    nfolds = 3,
    nrepeats = 3,
    sampling = NULL,
    tune_length = 5,
    top_n = 10,
    method = c("LR", "RF", "SVM"),
    ...) {
    sample_meta <- sample_data(ps)
    meta_nms <- names(sample_meta)
    if (!group %in% meta_nms) {
        stop(
            group, " are not contained in the `sample_data` of `ps`",
            call. = FALSE
        )
    }
    ps <- check_rank_names(ps)
    ps <- check_taxa_rank(ps, taxa_rank)

    # In current version, sl just for two groups comparisons
    groups <- sample_meta[[group]]
    group_n <- length(unique(groups))
    if (group_n != 2) {
        stop(
            "Supervised method only support for two groups comparisons",
            call. = FALSE
        )
    }

    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    method <- match.arg(method, choices = c("LR", "RF", "SVM"))
    full_method <- switch(method,
        LR = "logistic regression",
        RF = "random forest",
        SVM = "support vector machine"
    )
    train_method <- switch(method,
        LR = "glmnet",
        RF = "ranger",
        SVM = "svmLinear"
    )


    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # normalization, write a function here
    # fitZig fitFeatureModel
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)

    # summarize data
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    counts_tab <- abundances(ps_summarized, norm = TRUE)
    tax_tab <- as.data.frame(tax_table(ps_summarized))
    # in the animalcules, the counts were transferred as cpm, counts per million

    # number of markers must smaller than features
    feature_n <- nrow(tax_tab)
    if (feature_n < top_n) {
        stop(
            "There are ", feature_n, " features, ",
            "`top_n` must be smaller than number of features",
            call. = FALSE
        )
    }

    # transpose for modeling train
    counts_tab <- transpose_and_2df(counts_tab)
    colnames(counts_tab) <- tax_tab[, 1]

    # filter zero or near zero-variance predictors
    # https://topepo.github.io/caret/pre-processing.html#nzv
    # {stackoverflow}questions/47060233/stacking-models-from-different-packages

    # add target variable
    # to compute the varImp of svm model y must be a factor,
    # or there is a error: Error in y - mean(y, rm.na = TRUE): non-numeric
    # argument to binary operator
    counts_tab$y <- factor(groups)

    # set up classification model parameters
    fit_control <- caret::trainControl(
        method = "repeatedcv", # k fold cross validation
        number = nfolds,
        repeats = nrepeats,
        classProbs = TRUE,
        summaryFunction = caret::twoClassSummary, # is only for two classes
        sampling = sampling,
        savePredictions = TRUE
    )

    model_fit <- caret::train(
        y ~ .,
        data = counts_tab,
        method = train_method,
        tuneLength = tune_length,
        trControl = fit_control,
        metric = "ROC",
        ...
    )
    imp_df <- caret::varImp(model_fit)$importance
    # remove backtick
    feature <- gsub("`", "", row.names(imp_df))
    imp_df <- data.frame(
        feature = feature,
        imp = imp_df[, 1]
    )
    imp_df <- imp_df[order(imp_df$imp, decreasing = TRUE), ]
    marker <- imp_df[seq_len(top_n), ]

    # remove target variable
    counts_tab$y <- NULL

    # get the enrich_group
    marker_idx <- match(marker$feature, tax_tab[, 1])
    counts_tab_marker <- counts_tab[, marker_idx]
    enrich_group <- get_sl_enrich_group(
        counts_tab_marker,
        groups,
        sample_in_cols = FALSE
    )
    marker$enrich_group <- enrich_group
    marker <- marker[c("feature", "enrich_group", "imp")]
    names(marker) <- c("feature", "enrich_group", "ef_imp")
    ot <- otu_table(t(counts_tab), taxa_are_rows = TRUE)
    tt <- tax_table(ps_summarized)
    row.names(ot) <- row.names(tt)

    mm <- microbiomeMarker(
        marker_table = marker_table(marker),
        norm_method = get_norm_method(norm),
        diff_method = full_method,
        otu_table = ot,
        sam_data = sample_data(ps_summarized),
        tax_table = tt
    )

    mm
}

get_sl_enrich_group <- function(counts_tab, groups, sample_in_cols = TRUE) {
    if (sample_in_cols) {
        counts_tab <- t(counts_tab)
    }

    counts_mean <- by(counts_tab, groups, colMeans)
    counts_mean <- do.call(cbind, counts_mean)
    idx_enrich <- apply(counts_mean, 1, which.max)
    group_enrich <- colnames(counts_mean)[idx_enrich]

    group_enrich
}
