#' ROC curve of microbiome marker from supervised learning methods
#'
#' Show the ROC curve of the microbiome marker calculated by `run_sl`.
#'
#' @param mm a [microbiomeMarker-class] object.
#' @param group,nfolds,nrepeats,tune_length,... same with the `run_sl()`.
#'
#' @importFrom plotROC geom_roc style_roc calc_auc
#' @importFrom ggplot2 coord_equal annotate
#' @return  a [`ggplot2::ggplot`] object.
#' @seealso [`run_sl()`]
#' @export
#' @examples
#' data(enterotypes_arumugam)
#' # small example phyloseq object for test
#' ps_s <- phyloseq::subset_taxa(
#'     enterotypes_arumugam,
#'     Phylum %in% c("Firmicutes", "Bacteroidetes")
#' )
#'
#' set.seed(2021)
#' mm <- run_sl(
#'     ps_s,
#'     group = "Gender",
#'     taxa_rank = "Genus",
#'     nfolds = 2,
#'     nrepeats = 1,
#'     top_n = 15,
#'     norm = "TSS",
#'     method = "LR",
#' )
#' plot_sl_roc(mm, group = "Gender")
plot_sl_roc <- function(mm,
    group,
    nfolds = 3,
    nrepeats = 3,
    tune_length = 5,
    ...) {

    # sl method
    diff_method <- mm@diff_method
    diff_method <- switch(diff_method,
        "logistic regression" = "LR",
        "random forest" =  "RF",
        "support vector machine" = "SVM"
    )
    train_method <- switch(diff_method,
        LR = "glmnet",
        RF = "ranger",
        SVM = "svmLinear"
    )

    count_tab <- as(otu_table(mm), "matrix")
    marker <- marker_table(mm)
    count_marker <- count_tab[rownames(count_tab) %in% marker$feature, ]

    fitControl <- caret::trainControl(
        method = "repeatedcv",
        number = nfolds,
        repeats = nrepeats,
        classProbs = TRUE,
        summaryFunction = caret::twoClassSummary,
        savePredictions = TRUE
    )


    dat <- transpose_and_2df(count_marker)
    dat$y <- factor(sample_data(mm)[[group]])
    fit <- caret::train(
        y ~ .,
        data = dat,
        method = train_method,
        trControl = fitControl,
        tuneLength = tune_length,
        metric = "ROC",
        ...
    )

    prob_pred <- as.numeric(fit$pred$obs)
    prob_pred[prob_pred == 1] <- 0
    prob_pred[prob_pred == 2] <- 1
    df_roc <- data.frame(
        m = fit$pred[, colnames(fit$pred) == levels(fit$pred$obs)[2]],
        d = prob_pred,
        stringsAsFactors = FALSE
    )

    p <- ggplot(df_roc, aes(m = .data$m, d = .data$d)) +
        geom_roc(n.cuts = 0) +
        coord_equal() +
        style_roc()
    p <- p +
        annotate(
            "text",
            x = 0.75, y = 0.25,
            label = paste("AUC =", round((calc_auc(p))$AUC, 2))
        )

    p
}
