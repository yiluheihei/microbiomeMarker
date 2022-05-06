# The module of comparing differential analysis is inspired from DAtest
# https://github.com/Russel88/DAtest
# If you use this function please cite the original paper:
# Russel, Jakob, et al. "DAtest: a framework for choosing differential abundance
# or expression method." BioRxiv (2018): 241802.

#' Comparing the results of differential analysis methods by Empirical power
#' and False Discovery Rate
#'
#' Calculating power, false discovery rates, false positive rates and auc (
#' area under the receiver operating characteristic (ROC) curve)
#' for various DA methods.
#'
#' @param ps,group,taxa_rank main arguments of all differential analysis
#'   methods. `ps`: a [`phyloseq::phyloseq-class`] object; `group`, character,
#'   the variable to set the group, must be one of the var of the sample
#'   metadata; `taxa_rank`: character, taxonomic rank, please not that **since
#'   the abundance table is spiked in the lowest level, only
#'   `taxa_rank = "none"` is allowed**.
#' @param methods character vector, differential analysis methods to be
#'   compared, available methods are "aldex", "ancom", "ancombc", "deseq2",
#'   "edger", "lefse", "limma_voom", "metagenomeseq", "simple_stat".
#' @param args named list, which used to set the extra arguments of the
#'   differential analysis methods, so the names must be contained in `methods`.
#'   For more see details below.
#' @param n_rep integer, number of times to run the differential analyses.
#' @param effect_size numeric, the effect size for the spike-ins. Default 5.
#' @param k numeric vector of length 3, number of features to spike in each
#'   tertile (lower, mid, upper), e.g. `k=c(5,10,15)` means 5 features spiked
#'   in low abundance tertile, 10 features spiked in mid abundance tertile and
#'   15 features spiked in high abundance tertile. Default `NULL`, which will
#'   spike 2 percent of the total amount of features in each tertile (a total
#'   of 6 percent), but minimum c(5,5,5).
#' @param relative logical, whether rescale the total number of individuals
#'   observed for each sample to the original level after spike-in. Default
#'   `TRUE`.
#' @param BPPARAM [`BiocParallel::BiocParallelParam`] instance defining the
#'   parallel back-end.
#'
#' @details
#' To make this function support for different arguments for a certain DA method
#' `args` allows list of list of list e.g. `args = list(lefse = list(list(norm = "CPM"), list(norm = "TSS")))`, which specify to compare the different norm
#' arguments for lefse analysis.
#'
#' For `taxa_rank`, only `taxa_rank = "none"` is supported, if this argument is
#' not "none", it will be forced to "none" internally.
#'
#'
#' @return an `compareDA` object, which contains a two-length list of:
#'   - `metrics`: `data.frame`, FPR, AUC and spike detection rate for each run.
#'   - `mm`: differential analysis results.
#'
#' @importFrom phyloseq `otu_table<-`
#' @importFrom stats median
#' @export
compare_DA <- function(ps,
                       group,
                       taxa_rank = "none",
                       methods,
                       args = list(),
                       n_rep = 20,
                       effect_size = 5,
                       k = NULL,
                       relative = TRUE,
                       BPPARAM = BiocParallel::SnowParam(progressbar = TRUE)) {
    stopifnot(inherits(ps, "phyloseq"))
    # check methods
    avlb_methods <- c("aldex", "ancom", "ancombc", "deseq2", "edger", "lefse",
                      "limma_voom", "metagenomeseq", "simple_stat")
    out_methods <- setdiff(methods, avlb_methods)
    if (length(out_methods)) {
        stop("methods ", paste(out_methods, collapse = ", "),
             " not available. \n",
             "Please check your `methods`.\n",
             paste(strwrap(paste("Available methods:",
                           paste(avlb_methods, collapse = ", ")),
                           width = 0.9 * getOption("width")),
                   collapse = paste("\n", space(nchar("Available methods:")))),
             ".\n",
             call. = FALSE)
    }

    ps_var_name <- deparse(substitute(ps))

    # support for different arguments for a DA method, list of list of list
    # e.g. args = list(lefse = list(list(norm = "CPM"), list(norm = "TSS")))
    # different norm arguments for lefse analysis
    new_args <- generate_compare_args(methods, args)
    methods <- new_args$methods
    args <- new_args$args

    meta <- sample_data(ps)
    groups <- meta[[group]] |> factor()
    n_lvl <- nlevels(groups)

    if (n_lvl == 2) {
        if ("test_multiple_groups" %in% methods) {
            warning("There are two categories of interested variable ", group,
                    ", method `test_multiple_groups` are dropped.")
            methods <- setdiff(methods, "test_multiple_groups")
        }
    } else if (n_lvl >= 3) {
        if ("test_two_groups" %in% methods) {
            warning("There are more than two categories of interested variable ",
                    group,
                    ", method `test_two_groups` are dropped.")
            methods <- setdiff(methods, "test_two_groups")
        }
    } else {
        stop("Only one category  of interested variable ", group, ".")
    }

    # taxa_rank must be "none"
    if (taxa_rank != "none") {
        warning("since the abundance table is spiked in the lowest level, ",
             "`taxa_rank` was forced set as 'none'",
             call. = FALSE)
        taxa_rank <- "none"
    }

    # taxa_ranks <- vapply(args, `[[`, "taxa_rank", FUN.VALUE = character(1))
    # wrong_taxa_rank <- taxa_ranks != "none"
    # if (any(wrong_taxa_rank)) {
    #     warning("Set `taxa_rank` of all methods to 'none'")
    #     for (i in which(wrong_taxa_rank)) {
    #         args[[i]]$taxa_rank <- "none"
    #     }
    # }

    count_tab <- otu_table(ps)
    features <- rownames(count_tab)

    # spike in, differential features
    n_feature <- nrow(count_tab)

    if (is.null(k)) {
        k <- rep(round(n_feature * 0.02), 3)
        if (sum(k) < 15) {
            k <- rep(5, 3)
        }
    }
    if (sum(k) == n_feature) {
        stop("Set to spike all features, can't calculate FDR or AUC",
             call. = FALSE)
    }
    if (sum(k) > n_feature) {
        stop("Set to spike more features than are present in the data",
             call. = FALSE)
    }
    if (sum(k) < 15 && sum(k) >= 10 && n_rep <= 10) {
        warning("Few features are spiked, increase `k` or set `n_rep` to ",
                "more than 10 to ensure proper estimation of AUC and FPR",,
                call. = FALSE)
    }
    if (sum(k) < 10 && sum(k) >= 5 && n_rep <= 20) {
        warning("Few features are spiked, increase `k` or set `n_rep` to ",
                "more than 20 to ensure proper estimation of AUC and FPR",
                call. = FALSE)
    }
    if (sum(k) < 5 && n_rep <= 50) {
        warning("Very few features are spiked, increase `k` set `n_rep` to ",
                "more than 50 to ensure proper estimation of AUC and FPR",
                call. = FALSE)
    }
    if (sum(k) > n_feature/2) {
        warning("Set to spike more than half of the features, ",
                "which might give unreliable estimates")
   }

    # if(verbose) cat("Spikeing...\n")
    # shuffle predictor
    # predictor <- sample_data(ps)[[group[[1]]]]
    rands <- lapply(seq_len(n_rep), \(x) sample(groups))
    # spikeins
    spikeds <- lapply(rands,
                      \(x) spikein(count_tab, x, effect_size, k, relative))
    count_tabs <- lapply(spikeds, `[[`, 1)
    spiked_features <- lapply(spikeds, `[[`, 2)
    spiked_features <- rep(spiked_features, each = length(methods))

    # spiked phyloseq objects
    generate_spiked_ps <- function(spiked_count, rand, group, ps = ps) {
        otu_table(ps) <- otu_table(spiked_count, taxa_are_rows = TRUE)
        meta <- sample_data(ps)
        meta[[group]] <- rand
        sample_data(ps) <- meta

        ps
    }
    pss <- mapply(generate_spiked_ps, 
                  count_tabs, rands,
                  MoreArgs = list(group = group, ps = ps))
    pss <- rep(pss, each = length(methods))
    # rep methods
    rep_methods <- rep(methods, n_rep)
    all_methods <- methods

    # fun for performance metrics of DA methods
    for (item in names(args)) {
        args[[item]]$group <- group
        args[[item]]$taxa_rank <- taxa_rank
    }

    ## for ancombc, use bpmapply will raise an error:
    # `argument "formula" is missing, with no default`. This error could be due
    # to ANCOMBC package. We first run the ancombc method sequentially, and then
    # run other methods parallelly, finnally bind the results
    ancombc_md_idx <- grepl("ancombc", rep_methods, fixed = TRUE)
    if (any(ancombc_md_idx)) {
        ancombc_mds <- rep_methods[ancombc_md_idx]
        ancombc_args <- args[ancombc_md_idx]
        ancombc_pss <- pss[ancombc_md_idx]
        ancombc_spiked_features <- spiked_features[ancombc_md_idx]

        methods <- rep_methods[!ancombc_md_idx]
        args <- args[!ancombc_md_idx]
        pss <- pss[!ancombc_md_idx]
        spiked_features <- spiked_features[!ancombc_md_idx]
    }

    calc_da_metrics <- function(ps, method, args, features,
                                spiked_features, ps_var_name,
                                effect_size) {
        args <- args[[method]]
        args$ps <- ps
        fun <- paste0("run_", method)
        # remove number suffix for different args for a certern method
        fun <- gsub("(.*)_\\d+$", "\\1", fun)
        tm <- system.time(mm <- do.call(fun, args))

        marker <- data.frame(marker_table(mm))

        # pseduo pvalue of ancom
        if (method == "ancom" && !is.null(marker)) {
            w <- marker$W
            cf <- 0.05 * min(w)
        }

        if (!is.null(marker)) {
            spiked <- rep("no", nrow(marker))
            spiked[marker$feature %in% spiked_features] <- "yes"
            marker$spiked <- spiked
        }

        # confusion matrix
        n_pos <- ifelse(is.null(marker), 0, nrow(marker))
        n_feature <- length(features)
        neg_feature <- setdiff(features, spiked_features)
        n_neg <- n_feature - n_pos
        # for effect_size = 1
        true_neg <- n_neg
        true_pos <- 0
        false_neg <- 0
        # for effect size != 1
        if (effect_size != 1) {
            true_pos <- ifelse(is.null(marker), 0, sum(marker$spiked == "yes"))
            false_neg <- sum(neg_feature %in% spiked_features)
        }
        false_pos <- n_pos - true_pos
        true_neg <- n_neg - false_neg

        # fpr
        fpr <- ifelse((false_pos + true_neg) != 0,
                      false_pos / (false_pos + true_neg),
                      0)

        # sdr
        # sum(k) == length(spiked_features)
        sdr <- true_pos / length(spiked_features)

        # auc
        if(effect_size != 1) {
            test_roc <- NULL
            tryCatch(test_roc <- pROC::roc(
                as.numeric(marker$feature %in% spiked_features) ~ marker$pvalue,
                auc = TRUE,
                direction = ">",
                quiet = TRUE), error = function(e) NULL)
            auc <- ifelse(is.null(test_roc), 0.5, as.numeric(test_roc$auc))
        } else {
            auc <- 0.5
        }

        # fdr
        fdr <- ifelse(n_pos != 0, false_pos / n_pos, 0)

        # create call
        cmd_args <- args
        cmd_args$ps <- ps_var_name
        cmd_args$fun <- fun
        # reorder
        args_nms <- names(cmd_args)
        head_nms <- c("fun", "ps", "group", "taxa_rank")
        new_nms <- c(head_nms, setdiff(args_nms, head_nms))
        cmd_args <- cmd_args[new_nms]
        cmd_chr <- deparse1(as.call(cmd_args))
        cmd_chr <- gsub(paste0("\"(", fun, ")\""), "\\1", cmd_chr)

        metrics <- data.frame(auc = auc,
                              fpr = fpr,
                              fdr = fdr,
                              power = sdr,
                              # method = gsub("(.*)_\\d+$", "\\1", method),
                              method = method,
                              call = cmd_chr,
                              time_min = round(tm[3] / 60, 4))
        rownames(metrics) <- NULL

        list(metrics = metrics, mm = mm)
    }

    if (any(ancombc_md_idx)) {
        ancombc_out <- mapply(calc_da_metrics,
                              ps = ancombc_pss,
                              method = ancombc_mds,
                              spiked_features = ancombc_spiked_features,
                              MoreArgs = list(args = ancombc_args,
                                              features = features,
                                              ps_var_name = ps_var_name,
                                              effect_size = effect_size),
                              SIMPLIFY = FALSE)
    }
    if (!all(ancombc_md_idx)) {
        bp_out <- BiocParallel::bpmapply(calc_da_metrics,
                               ps = pss,
                               method = methods,
                               spiked_features = spiked_features,
                               MoreArgs = list(args = args,
                                               features = features,
                                               ps_var_name = ps_var_name,
                                               effect_size = effect_size),
                               BPPARAM = BPPARAM,
                               SIMPLIFY = FALSE)
    }

    if (all(ancombc_md_idx)) {
        da_out <- ancombc_out
    } else {
        if (any(ancombc_md_idx)) {
            da_out <- c(ancombc_out, bp_out)
        } else {
            da_out <- bp_out
        }
    }


    # order the out, to the original order
    idx <- order(c(which(ancombc_md_idx), which(!ancombc_md_idx)))
    da_out <- da_out[idx]

    da_metrics <- do.call(rbind, lapply(da_out, `[[`, "metrics"))
    da_metrics$run <- rep(seq_len(n_rep), each = length(all_methods))
    mms <- lapply(da_out, `[[`, "mm")

    # detail <- list(n_feature = n_feature,
    #                n_sample = ncol(count_tab),
    #                effect_size = effect_size,
    #                spike = paste0(c("Low:","Mid:","High:"), k,
    #                               collapse = ", "))

    out <- list(metrics = da_metrics, mm = mms)
    class(out) <- "compareDA"

    out
}


#' Summary differential analysis methods comparison results
#'
#' @param object an `compareDA` object, output from [`compare_DA()`].
#' @param sort character string specifying sort method. Possibilities are
#'   "score" which is calculated as \eqn{(auc - 0.5) * power - fdr}, "auc" for
#'   area under the ROC curve, "fpr" for false positive rate, "power" for
#'   empirical power.
#' @param boot logical, whether use bootstrap for confidence limites of the
#'   score, default `TRUE`. Recommended to be `TRUE` unless `n_rep` is larger
#'   then 100 in [`compare_DA()`].
#' @param boot_n integer, number of bootstraps, default 1000L.
#' @param prob two length numeric vector, confidence limits for score, default
#'   `c(0.05, 0.95)`.
#' @param ... extra arguments affecting the summary produced.
#' @return a `data.frame` containing measurements for differential analysis
#'    methods:
#'    - `call`: differential analysis commands.
#'    - `auc`: area under curve of ROC.
#'    - `fpr`: false positive rate
#'    - `power`: empirical power.
#'    - `fdr`: false discover7y rate.
#'    - `score`: score whch is calculated as \eqn{(auc - 0.5) * power - fdr}.
#'    - `score_*`: confidence limits of score.
#' @export
summary.compareDA <- function(object,
                              sort = c("score", "auc", "fpr", "power"),
                              boot = TRUE,
                              boot_n = 1000L,
                              prob = c(0.05, 0.95),
                              ...) {
    stopifnot(inherits(object, "compareDA"))
    sort <- match.arg(sort, c("score", "auc", "fpr", "power"))

    # medians
    metrics <- object$metrics
    calls <- metrics$call
    new_metrics <- metrics[c('auc', 'fpr', 'power', 'fdr')]
    metric_med <- stats::aggregate(new_metrics,
                                   by = list(call = calls),
                                   FUN = median)
    # score
    metric_med$score <- (metric_med$auc - 0.5) * metric_med$power -
        metric_med$fdr
    # interval
    new_metrics$score <- (new_metrics$auc - 0.5) * new_metrics$power -
        new_metrics$fdr
    metrics$score <- (metrics$auc - 0.5) * metrics$power - metrics$fdr

    if (boot) {
        boots <- dplyr::group_by(metrics, call) |>
            dplyr::group_modify(~ .x[sample(rownames(.x),
                                            boot_n,
                                            replace = TRUE),
                                    ])
        score_cl <- stats::aggregate(score ~ call,
                                  data = boots,
                                  FUN = function(x)
                                      stats::quantile(x, probs = prob))
    } else {
        score_cl <- stats::aggregate(score ~ call,
                                     data = metrics,
                                     FUN = \(x) stats::quantile(x, probs = prob))
    }
    score_cl <- data.frame(call = score_cl$call,
                           score_cl$score[, 1],
                           score_cl$score[, 2])
    names(score_cl) <- c("call", paste0("score_", prob))
    out <- merge(metric_med, score_cl, by = "call")
    # reorder score descreasing
    out <- out[order(out$score,
                     out[[paste0("score_", prob[1])]],
                     out[[paste0("score_", prob[2])]],
                     decreasing = TRUE),
               ]
    if (out$score[1] <= 0) {
        warning("Best score is <= 0.\n",
                "You might require to preprocessing your data or ",
                "re-run with a higher effect size.",
                call. = FALSE)
    }

    # mat <- vector("logical", nrow(out))
    # for (i in seq_len(nrow(out))) {
    #     mat[i] <- out$score[i] >= out[[paste0("score_", prob[1])]][1]
    # }
    # out$` ` <- " "
    # out[mat,]$` ` <- "*"


    if (sort == "auc") {
        out <- out[order(out$auc, decreasing = TRUE), ]
    }
    if (sort == "fpr") {
        out <- out[order(out$fpr, decreasing = FALSE), ]
    }
    if (sort == "power") {
        out <- out[order(out$power, decreasing = TRUE), ]
    }

    out
}

# compare_DA <- function(...,
#                        n_rep = 20,
#                        effect_size = 5,
#                        k = NULL,
#                        n_core = parallel::detectCores() -1,
#                        check_core = TRUE,
#                        relative = TRUE,
#                        verbose = TRUE) {
#     if (check_core) {
#         if (check_core > 20) {
#             ANSWER <- readline(paste("You are about to run compareDA using",
#                                      n_core,
#                                      "cores. Enter y to proceed "))
#             if (ANSWER != "y") {
#                 stop("Process aborted")
#             }
#         }
#     }
#
#     exp_chrs <- list(...)
#     t_start <- proc.time()
#     calls <- lapply(exp_chrs, \(x) standardise_call(str2lang(x)))
#
#     # extract the ps object and target variable
#     pss <- lapply(calls, `[[`, "ps")
#     group <- lapply(calls, `[[`, "group")
#     # all the ps and target variable must be the same
#     if (sum(duplicated(pss)) != (length(pss) - 1)) {
#         stop("`ps` objects in DA analysis must be the same")
#     }
#     if (sum(duplicated(group)) != (length(group) - 1)) {
#         stop("`group` var in all DA analysis must be the same")
#     }
#
#     # DA methods comparison only support for taxa_rank = "none", since
#     #  the abundance table is spiked in the lowest level
#     # full_paras <- lapply(calls, \(x) formals(match.fun(x[[1]])))
#     taxa_ranks <- lapply(calls, `[[`, "taxa_rank")
#     if (sum(duplicated(taxa_ranks)) != (length(group) - 1)) {
#         stop("`taxa_rank` objects in all DA analysis must be the same")
#     }
#     if (!is.null(taxa_ranks[[1]]) && taxa_ranks[[1]] != "none") {
#         stop("since the abundance table is spiked in the lowest level, ",
#              "`taxa_rank` must be 'none'",
#              call. = FALSE)
#     }
#
#      if (verbose) {
#         message("Comparing differential methods may take a long time")
#         message("Running on ", n_core, "cores")
#     }
#
#     # differential analysis functions
#     funs <- lapply(calls, `[[`, 1)
#     funs_chr <- vapply(funs, as.character, FUN.VALUE = character(1))
#
#     ps <- eval(pss[[1]], envir = parent.frame())
#     count_tab <- otu_table(ps)
#     features <- rownames(count_tab)
#
#     # spike in differential features
#     if (is.null(k)) {
#         k <- rep(round(nrow(count_tab) * 0.02), 3)
#         if (sum(k) < 15) {
#             k <- rep(5, 3)
#         }
#     }
#     n_feature <- nrow(count_tab)
#     if (sum(k) == n_feature) {
#         stop("Set to spike all features, can't calculate FDR or AUC",
#              call. = FALSE)
#     }
#     if (sum(k) > n_feature) {
#         stop("Set to spike more features than are present in the data",
#              call. = FALSE)
#     }
#     if (sum(k) < 15 && sum(k) >= 10 && n_rep <= 10) {
#         warning("Few features are spiked, increase `k` or set `n_rep` to ",
#                 "more than 10 to ensure proper estimation of AUC and FPR",,
#                 call. = FALSE)
#     }
#     if (sum(k) < 10 && sum(k) >= 5 && n_rep <= 20) {
#         warning("Few features are spiked, increase `k` or set `n_rep` to ",
#                 "more than 20 to ensure proper estimation of AUC and FPR",
#                 call. = FALSE)
#     }
#     if (sum(k) < 5 && n_rep <= 50) {
#         warning("Very few features are spiked, increase `k` set `n_rep` to ",
#                 "more than 50 to ensure proper estimation of AUC and FPR",
#                 call. = FALSE)
#     }
#     if (sum(k) > n_feature/2) {
#         warning("Set to spike more than half of the features, ",
#                 "which might give unreliable estimates")
#     }
#
#     if(verbose) cat("Spikeing...\n")
#     # shuffle predictor
#     predictor <- sample_data(ps)[[group[[1]]]]
#     rands <- lapply(seq_len(n_rep), \(x) sample(predictor))
#
#     # spikeins
#     spikeds <- lapply(rands,
#                       \(x) spikein(count_tab, x, effect_size, k, relative))
#     count_tabs <- lapply(spikeds, `[[`, 1)
#
#     if (verbose) {
#         cat(paste("Testing", length(exp_chrs),
#                   "methods", n_rep, "times each ...\n"))
#     }
#
#     # progress bar
#     # da_par <- paste(rep(seq_len(n_rep), each = 2),
#     #                funs_chr, sep = "-")
#     cmds <- rep(exp_chrs, n_rep)
#     run_no <- rep(seq_len(n_rep), each = length(exp_chrs))
#     pb <- utils::txtProgressBar(max = length(cmds), style = 3)
#     progress <- function(n) setTxtProgressBar(pb, n)
#     opts <- list(progress = progress)
#
#     # config parallel
#     if (n_core == 1) {
#         foreach::registerDoSEQ()
#     } else {
#         cl <- parallel::makeCluster(n_core)
#         doSNOW::registerDoSNOW(cl)
#         on.exit(parallel::stopCluster(cl))
#     }
#
#     # run the DA analysis in parallel
#     res <- foreach::foreach(exp_chr = cmds, i = run_no,
#                             .export = c("otu_table", "otu_table<-", funs_chr),
#                             .options.snow = opts) %dopar% {
#         t1_sub <- proc.time()
#         # construct new ps with spiked feature abundance table
#         new_count_tab <- count_tabs[[i]]
#         otu_table(ps) <- otu_table(new_count_tab, taxa_are_rows = TRUE)
#         res_sub <- eval(str2expression(exp_chr),
#                         list(ps = ps),
#                         enclos = parent.frame())
#
#         run_time_sub <- (proc.time() - t1_sub)[3]
#         return(list(res_sub, run_time_sub))
#     }
#     run_times <- lapply(res, `[[`, "run_time_sub")
#     da_res <- lapply(res, `[[`, "res_sub")
#
#     n_da <- length(exp_chrs)
#     r <- NULL
#     final_res <- foreach::foreach(r = seq_len(n_rep)) %do% {
#         da_sub <- da_res[(1 + (r - 1) * n_da):(r * n_da)]
#         curr_cmds <- cmds[(1 + (r - 1) * n_da):(r * n_da)]
#         curr_spiked_features <- spikeds[[r]][[2]]
#         # insert spiked column
#         rsp <- NULL
#         da_sub <- foreach::foreach(rsp = seq_along(da_sub)) %do% {
#             tmp <- da_sub[[rsp]]
#             tmp_marker <- data.frame(marker_table(tmp))
#
#             # psedudo pvalue of ancom
#             if (grepl("ancom(", curr_cmds[rsp], fixed = TRUE)) {
#                 w <- tmp_marker$W
#                 cf <- 0.05 * min(w)
#                 tmp_marker$pvalue <- (1 / w) * cf
#             }
#
#             tmp_spiked <- rep("no", nrow(tmp_marker))
#             tmp_spiked[tmp_marker$feature %in% curr_spiked_features] <- "yes"
#             tmp_marker$spiked <- tmp_spiked
#             return(tmp_marker)
#         }
#
#         # confusion matrix
#         n_pos <- vapply(da_sub, nrow, FUN.VALUE = integer(1))
#         neg_feature <- lapply(da_sub, \(x) setdiff(features, x$feature))
#         n_neg <- n_feature - n_pos # vapply(neg_feature, length, integer(1))
#         # for effect_size = 1
#         true_neg <- n_neg
#         true_pos <- 0
#         false_neg <- 0
#         # for effect_size != 1
#         if (effect_size != 1) {
#             true_pos <- vapply(da_sub,
#                                \(x) sum(da_sub$spiked == "yes"),
#                                FUN.VALUE = integer(1))
#             false_neg <- vapply(neg_feature,
#                                \(x) sum(x %in% curr_spiked_features),
#                                FUN.VALUE = integer(1))
#         }
#         false_pos <- n_pos - true_pos
#         true_neg <- n_neg - false_neg
#
#         # FPR: false positive rate
#         fprs <- vapply(seq_along(da_sub),
#                        \(x) ifelse((false_pos[x] + true_neg[x]) != 0,
#                                    false_pos[x] / (false_pos[x] + true_neg[x]),
#                                    0),
#                        FUN.VALUE = numeric(1))
#
#         # sdr: spike detection rate
#         sdrs <- vapply(true_pos, \(x) x / sum(k), FUN.VALUE = numeric(1))
#
#         # auc
#         aucs <- vapply(da_sub, \(x) {
#             if (effect_size != 1) {
#                 test_roc <- NULL
#                 spiked_idx <- as.numeric(x$feature %in% curr_spiked_features)
#                 tryCatch(
#                     test_roc <- pROC::roc(spiked_idx ~ x$pvalue,
#                                           auc = TRUE,
#                                           direction = ">",
#                                           quiet = TRUE),
#                     error = function(e) NULL)
#                     res <- ifelse(is.null(test_roc),
#                                   0.5,
#                                   as.numeric(test_roc$auc))
#             } else {
#                 res <- 0.5
#             }
#
#             res
#         }, FUN.VALUE = numeric(1))
#
#         # fdrs
#         fdrs <- vapply(seq_along(da_sub),
#                        \(x) ifelse(n_pos != 0, false_pos / n_pos, 0))
#
#         df_combine <- data.frame(call = curr_cmds,
#                                  AUC = aucs,
#                                  FPR = fprs,
#                                  FDR = fdrs,
#                                  Power = sdrs,
#                                  run = r)
#         rownames(df_combine) <- NULL
#
#         return(df_combine, da_sub)
#     }
#
#     out_res <- do.call(rbind, lapply(final_res, `[[`, 1))
#     out_res_marker <- lapply(final_res, `[[`, 2)
#
#     # running time
#     run_secs <- (proc.time() - t_start)[3]
#     if ((run_secs)/60/60 > 1) {
#         run_time <- paste(round((run_secs)/60/60, 2), "Hours")
#     } else {
#         run_time <- paste(round((run_secs)/60,2),"Minutes")
#     }
#
#     out_detail <- data.frame(n_feature = nrow(count_tab),
#                              n_sample = ncol(count_tab),
#                              run_time = run_time,
#                              effect_size = effect_size,
#                              spiked = paste0(c("Low:","Mid:","High:"), k,
#                                              collapse = ", "))
#     out_detail <- as.data.frame(t(out_detail))
#     names(out_detail) <- NULL
#
#     # run times
#     run_times <- data.frame(DA = cmds,
#                             minites = round(unlist(run_times) / 60, 4))
#
#     out <- list(res = out_res,
#                 marker = out_res_marker,
#                 detail = out_detail,
#                 run_time = run_times)
#
#     out
# }


# spike in features
spikein <- function(count_tab,
                    predictor,
                    effect_size = 2,
                    k,
                    relative = TRUE) {
    if (effect_size < 0) {
        stop("Effect size should be positive")
    }

    spike_method <- ifelse(effect_size == 1, "none", "mult")
    if (is.null(rownames(count_tab))) {
        rownames(count_tab) <- seq_len(nrow(count_tab))
    }
    count_tab <- as.data.frame(count_tab)
    predictor <- as.numeric(as.factor(predictor)) - 1

    # Choose Features to spike
    propcount <- sweep(count_tab, 2, colSums(count_tab), "/")
    # propcount <- apply(count_tab, 2, function(x) x/sum(x))
    count_abundances <- sort(rowSums(propcount)/ncol(propcount))

    # Only spike Features present in cases (except if predictor is numeric)
    case_count_tab <- count_tab[
        rowSums(count_tab[, predictor == 1]) > 0, predictor == 1]
    approved_count_abundances <- count_abundances[
        names(count_abundances) %in% row.names(case_count_tab)]

    # Which to spike in each tertile
    lower_tert <- names(approved_count_abundances[
        approved_count_abundances < quantile(approved_count_abundances,1/3)])
    mid_tert <- names(approved_count_abundances[
      approved_count_abundances >= quantile(approved_count_abundances,1/3) &
          approved_count_abundances < quantile(approved_count_abundances,2/3)])
    upper_tert <- names(approved_count_abundances[
        approved_count_abundances >= quantile(approved_count_abundances,2/3)])

    spike_features <- c(sample(lower_tert, k[1]),
                        sample(mid_tert, k[2]),
                        sample(upper_tert,k[3]))
    spike_feature_index <- which(row.names(count_tab) %in% spike_features)

    # Spike Features by multiplication
    old_sums <- colSums(count_tab)

    if (spike_method == "mult"){
        count_tab[spike_feature_index, predictor==1] <-
            count_tab[spike_feature_index, predictor==1] * effect_size
    }

    # Rescale to original sample sums
    new_sums <- colSums(count_tab)
    if (relative) {
        count_tab <- round(sweep(count_tab, 2, old_sums/new_sums, "*"))
    }

    list(count_tab, spike_features)
}


# from pryr: Standardise a function call
standardise_call <- function(call, env = parent.frame()) {
    stopifnot(is.call(call))
    f <- eval(call[[1]], env)
    if (is.primitive(f)) {
        return(call)
    }

    return(match.call(f, call))
}

# To make compare_DA() support for different arguments for a certain DA method,
# args allows list of list of list
# e.g. args = list(lefse = list(list(norm = "CPM"), list(norm = "TSS"))),
# represents compare the different norm arguments for lefse analysis. So we
# need to flattern the args and extend methods for DA analysis:
# methods = c("lefse", "lefse"),
# args = list(list(norm = "CPM"), list(norm = "TSS"))
#
# For method with no args provided, set it to list(), e.g. list(ancom = list()).
generate_compare_args <- function(methods, args) {
    # check args
    args_nms <- names(args)
    if (length(args)) {
        out_args <- setdiff(args_nms, methods)
        if (length(out_args)) {
            stop("names of `args` must be contained in `methods`.\n",
                paste(args[out_args], collapse = ", "), " in names of `args` ",
                "does not match DA methods",
                call. = FALSE)
        }
    }

    # create args list for each method
    method_no_args <- setdiff(methods, args_nms)
    for (i in seq_along(method_no_args)) {
        args[[method_no_args[i]]] <- list()
    }

    new_args <- list()
    n_arg <- vector("integer", length(args))
    for (i in seq_along(args)) {
        curr_arg <- args[i]
        if (purrr::vec_depth(curr_arg) > 4) {
            stop("`args` could be 'list of list', ",
                 "'list of list of list' to support for different arguments ",
                 "for a certain DA method")
        }
        if (purrr::vec_depth(curr_arg) == 4) {
            curr_arg <- unlist(curr_arg, recursive = FALSE)
            names(curr_arg) <- paste(names(args)[i],
                                     seq_along(curr_arg),
                                     sep = "_")

        }
        new_args <- c(new_args, curr_arg)
        n_arg[i] <- length(curr_arg)
    }
    methods <- rep(methods, times = n_arg)
    methods_suffix <- lapply(n_arg, \(x) {
        if (x > 1) {
            as.character(paste0("_", seq_len(x)))
        } else {
            ""
        }
    }) |> unlist()
    methods <- paste(methods, methods_suffix, sep = "")

    return(list(methods = methods, args = new_args))
}
