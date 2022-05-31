# enrich group of the feature ---------------------------------------------

#' get the mean abundances of each class for each feature
#' @noRd
get_feature_enrich_group <- function(class, feature) {
    feature$class <- class
    feature_mean <- group_by(feature, class) %>%
        group_modify(~ purrr::map_df(.x, mean)) %>%
        ungroup()

    feature_enrich_index <- select(feature_mean, -class) %>%
        purrr::map_dbl(which.max)
    feature_enrich_group <- feature_mean$class[feature_enrich_index]
    names(feature_enrich_group) <- names(feature)[names(feature) != "class"]
    feature_max_mean <- purrr::map2_dbl(
        select(feature_mean, -class),
        feature_enrich_index,
        ~ .x[.y]
    )
    feature_max_mean[feature_max_mean < 1] <- 1

    return(list(
        group = feature_enrich_group,
        log_max_mean = log10(feature_max_mean)
    ))
}


# bootstrap iteration LDA -------------------------------------------------

#' check bootstrap sample, not contast within class and too fewer in each
#' class
#'
#' @param feature_abudance data.frame, significant feature abundance,
#'   where columns represents samples and rows represents features
#' @param sample_indx numeric vector, sample index be bootstraped
#' @param sample_min integer, min samples in each class
#' @param class character vector, class of all samples
#'
#' @noRd
check_bootstrap_sample <- function(feature_abundance,
    sample_indx,
    sample_min,
    class) {
    if ("class" %in% names(feature_abundance)) {
        feature_abundance$class <- NULL
    }
    feature_abundance <- feature_abundance[sample_indx, ]
    class_n <- length(unique(class))
    class <- class[sample_indx]

    if (length(unique(class)) < class_n) {
        return(FALSE)
    }

    for (cls in unique(class)) {
        if (sum(class == cls) < sample_min) {
            return(FALSE)
        }

        # sig feature smaller than min sample count
        cls_abundance <- feature_abundance[class == cls, ]

        for (i in seq_along(ncol(cls_abundance))) {
            unique_abd <- length(unique(cls_abundance[[i]]))

            if ((unique_abd <= sample_min && sample_min > 1) ||
                (unique_abd <= 1 && sample_min == 1)) {
                return(FALSE)
            }
        }
    }

    return(TRUE)
}

#' bootstrap iteration of samples for lda analysis
#' @noRd
bootstap_lda <- function(feature_abundance,
    boot_n,
    class,
    sample_fract,
    seed = 2020) {
    # Bioconductor not allows set.seed
    ldas <- purrr::rerun(
        boot_n,
        bootstap_lda_one(
            feature_abundance,
            class,
            sample_fract
        )
    ) %>%
        purrr::transpose() %>%
        purrr::map(~ do.call(bind_rows, .x)) %>%
        bind_rows()

    mean_lds <- colMeans(ldas)
    mean_lds <- sign(mean_lds) * log10(1 + abs(mean_lds))

    mean_lds
}

bootstap_lda_one <- function(feature_abundance,
    class,
    sample_fract) {
    sample_groups <- unique(class)
    class_count <- table(class)
    feature_abundance$class <- class
    feature_abundance <- preprocess_feature_all(feature_abundance, class)

    sample_n <- nrow(feature_abundance)
    random_n <- floor(sample_n * sample_fract)
    class_n <- length(sample_groups)
    sample_min <- floor(
        min(class_count) * sample_fract * sample_fract * 0.5) %>%
        max(1)

    # class vs class
    pairs <- utils::combn(sample_groups, 2, simplify = FALSE) %>%
        purrr::map(sort, decreasing = TRUE)

    for (i in seq_len(1000)) {
        # random select samples using bootstrap method
        sample_indx <- sample(sample_n, random_n, replace = TRUE)

        is_checked <- check_bootstrap_sample(
            feature_abundance,
            sample_indx,
            sample_min,
            class
        )
        if (is_checked) {
            break
        }
    }

    if (!is_checked) {
        stop(
            "Too small samples in each class",
            " or the variance of feature abundances within a",
            " class too small (zero or near zero)",
            call. = FALSE
        )
    }

    lda <- purrr::map(
        pairs,
        ~ cal_pair_lda(feature_abundance, sample_indx, .x)
    )
    names(lda) <- purrr::map(pairs, paste, collapse = " -VS- ")

    lda
}

#' calculate lda score of single pair groups
#' @noRd
cal_pair_lda <- function(feature_abundance,
    sample_indx,
    pair) {
    sample_feature_abundance <- feature_abundance[sample_indx, ]

    # reference lefse.py in lefse
    lda_res <- suppressWarnings(
        MASS::lda(
            class ~ .,
            data = sample_feature_abundance,
            tol = 1.0e-10
        )
    )
    w <- lda_res$scaling[, 1]
    w_unit <- w / sqrt(sum(w^2))
    feature_remove_class <- sample_feature_abundance[-1]

    # not support subclass and subject argument in lefse

    ld <- as.matrix(feature_remove_class) %*% w_unit
    group1_indx <- sample_feature_abundance$class == pair[1]
    effect_size <- abs(mean(ld[group1_indx]) - mean(ld[-group1_indx]))
    wfinal <- w_unit * effect_size
    lda_means <- lda_res$means
    lda_row_nms <- row.names(lda_means)
    feature_n <- ncol(lda_means)
    coeff <- ifelse(is.nan(wfinal), 0, abs(wfinal))

    res <- purrr::map(
        pair,
        function(x) {
            if (x %in% lda_row_nms) {
                # fixes #7, Since `pair` is a level, and `lda_means[pair[i], ]` 
                # corced pair[i]` to numeric rather than use the corresponding 
                # level of pair[i]
                ind <- match(x, lda_row_nms)
                lda_means[ind, ]
            } else {
                rep(0, feature_n)
            }
        }
    )

    names(res) <- pair

    feature <- names(feature_remove_class)
    lda_score <- purrr::map_dbl(
        seq_along(feature),
        function(i) {
            gm <- abs(res[[1]][i] - res[[2]][i])
            return(gm + coeff[i] * 0.5)
        }
    )
    names(lda_score) <- feature

    lda_score
}


#' feature abundance preprocess
#' @noRd
preprocess_feature_all <- function(x, class) {
    res <- group_by(x, class) %>%
        group_modify(~ purrr::map_df(.x, preprocess_feature)) %>%
        ungroup()

    res
}

preprocess_feature <- function(x) {
    if (length(unique(x)) <= max(length(x) * 0.5, 4)) {
        x <- purrr::map_dbl(x, ~ abs(.x + rnorm(1, 0, max(.x * 0.05, 0.01))))
    }

    x
}

# wilcoxon test within the same class-------------------------------------------

#' wilcoxon test in each subclass
#'
#' @param subcls character vector, length equal to the number of samples
#' @param cls_hie list, length equal to the number of classes, class hierarchy
#' @param feats_abd numeric vector, abundance profile of a given feature
#' @param feats_name character, feature names
#' @param wilcoxon_cutoff the cutoff for the wilcoxon test, default 0.05
#' @param multicls_strat logical, for multiple class tasks, whether the test is
#'   performed in a one-against one (more strict) or in a one-against all
#'   setting, default `FALSE`.
#' @param strict multiple testing options, 0 for no correction (default), 1 for
#'   independent comparisons, 2 for independent comparison
#' @param sample_min integer, minimum number of samples per subclass for
#'   performing wilcoxon test, default 10
#' @param only_same_subcls logical, whether perform the wilcoxon test only
#'   among the subclasses with the same name, default `FALSE`
#' @param curv logical, whether perform the wilcoxon test using the
#'   Curtis's approach, defalt `FALSE`
#'
#' @noRd
test_rep_wilcoxon <- function(subcls,
    cls_hie,
    feats_abd,
    feats_name,
    strict = 0,
    wilcoxon_cutoff = 0.05,
    multicls_strat = FALSE,
    sample_min = 10,
    only_same_subcls = FALSE,
    curv = FALSE) {
    if (!strict %in% c(0, 1, 2)) {
        stop("`strict` must be 0, 1 or 2")
    }

    cls_nms <- names(cls_hie)
    pairs <- utils::combn(cls_nms, 2, simplify = FALSE)
    tot_ok <- 0

    all_diff <- list()

    for (pair in pairs) {
        dir_cmp <- "not_set"

        subcls1 <- cls_hie[[pair[1]]]
        subcls1_n <- length(subcls1)
        subcls2 <- cls_hie[[pair[2]]]
        subcls2_n <- length(subcls2)

        # multiple tests
        if (strict != 0) {
            wilcoxon_cutoff <- ifelse(
                strict == 2,
                wilcoxon_cutoff * subcls1_n * subcls2_n,
                1 - (1 - wilcoxon_cutoff)^(subcls1_n * subcls2_n)
            )
        }

        ok <- 0
        curv_sign <- 0
        first <- TRUE

        for (i in seq_along(subcls1)) {
            br <- FALSE
            for (j in seq_along(subcls2)) {
                if (only_same_subcls &&
                    gsub(pair[1], "", subcls1[i]) != 
                        gsub(pair[2], "", subcls2[j])) {
                    ok <- ok + 1
                    next
                }

                cls1_abd <- feats_abd[subcls == subcls1[i]]
                cls2_abd <- feats_abd[subcls == subcls2[j]]
                med_comp <- FALSE

                if (length(cls1_abd) < sample_min || 
                        length(cls2_abd) < sample_min) {
                    med_comp <- TRUE
                }

                sx <- stats::median(cls1_abd)
                sy <- stats::median(cls2_abd)

                if (cls1_abd[1] == cls2_abd[1] &&
                    length(unique(cls1_abd)) == 1 &&
                    length(unique(cls2_abd)) == 1
                ) {
                    tres <- FALSE
                    first <- FALSE
                } else if (!med_comp) {
                    x <- c(cls1_abd, cls2_abd)
                    y <- factor(
                        c(
                            rep(1, length(cls1_abd)),
                            rep(2, length(cls2_abd))
                        )
                    )
                    pv <- coin::wilcox_test(
                        x ~ y,
                        data = data.frame(x, y)
                    )
                    pv <- coin::pvalue(pv)
                    tres <- pv < wilcoxon_cutoff * 2
                }

                if (first) {
                    first <- FALSE

                    if (!curv && (med_comp || tres)) {
                        dir_cmp <- sx < sy
                    } else if (curv) {
                        dir_cmp <- NULL
                        if (med_comp || tres) {
                            curv_sign <- curv_sign + 1
                            dir_cmp <- sx < sy
                        }
                    } else {
                        br <- TRUE
                    }
                } else if (!curv && med_comp) {
                    if ((sx < sy) != dir_cmp || sx == sy) {
                        br <- TRUE
                    }
                } else if (curv) {
                    if (tres && is.null(dir_cmp)) {
                        curv_sign <- curv_sign + 1
                        dir_cmp <- sx < sy
                    }

                    if (tres && dir_cmp != (sx < sy)) {
                        br <- TRUE
                        curv_sign <- curv_sign - 1
                    }
                } else if (!tres || (sx < sy) != dir_cmp || sx == sy) {
                    br <- TRUE
                }

                if (br) {
                    break
                }
                ok <- ok + 1
            }
            if (br) {
                break
            }
        }

        diff <- ifelse(
            curv,
            curv_sign > 0,
            ok == subcls1_n * subcls2_n
        )

        if (diff) tot_ok <- tot_ok + 1
        if (!diff && multicls_strat) {
            return(FALSE)
        }
        if (diff && !multicls_strat) all_diff <- c(all_diff, pair)
    }

    if (!multicls_strat) {
        tot_k <- length(cls_hie)

        for (k in names(cls_hie)) {
            nk <- 0
            for (a in all_diff) {
                if (k %in% a) nk <- nk + 1
            }

            if (nk == tot_k - 1) {
                return(TRUE)
            }
        }

        return(FALSE)
    }

    return(TRUE)
}

# format input ------------------------------------------------------------

#' format lefse input
#'
#' @param sample_meta a data.frame like object, sample metadata
#' @param cls character variable of class
#' @param subcls character variable of subclass, default `NULL`, no subclass
#'
#' @noRd
#'
#' @return  a list, contains class, subclass, and class hierarchy
lefse_format_grp <- function(sample_meta, group, subgroup = NULL) {
    groups <- sample_meta[[group]]
    group_nms <- unique(groups)

    if (is.null(subgroup)) {
        subgroups <- paste0(groups, "_subgrp")
    } else {
        subgroups <- paste(groups, sample_meta[[subgroup]], sep = "_")
    }

    group_hie <- split(subgroups, groups) %>%
        purrr::map(unique)

    return(list(group = groups, subgroup = subgroups, group_hie = group_hie))
}

#' add missing levels, used for summarized taxa
#' @param feature feature data, a [phyloseq::otu_table-class]
#' @noRd
#' @description this function require the row names of `feature` is the
#' summarized taxa
#' @return a data frame, where taxa in rows
add_missing_levels <- function(feature) {
    if (!taxa_are_rows(feature)) {
        feature <- t(feature)
    }
    feature_nms <- row.names(feature)
    feature <- feature@.Data %>% data.frame()

    # the missing feature names
    feature_nms2 <-
        strsplit(feature_nms, "|", fixed = TRUE) %>%
        purrr::map(
            ~ Reduce(
                function(x, y) paste(x, y, sep = "|"),
                .x,
                accumulate = TRUE
            )
        )

    unq_nms <- unlist(feature_nms2) %>% unique()
    missing_nms <- setdiff(unq_nms, feature_nms)

    if (length(missing_nms) == 0) {
        return(feature)
    }

    missing_nms_split <- strsplit(missing_nms, split = "|", fixed = TRUE)
    missing_mns_level <- lengths(missing_nms_split)
    missing_level_range <- range(missing_mns_level)

    # only sum the next level of tax, so we need first add the missing tax at
    # the most depth level
    for (i in missing_level_range[2]:missing_level_range[1]) {
        missing_nms_i <- missing_nms[missing_mns_level == i]
        taxs <- row.names(feature)

        indx <- purrr::map(
            missing_nms_i,
            ~ purrr::map_lgl(taxs, function(x) grepl(.x, x, fixed = TRUE))
        )

        # only sum the next level of tax
        feature_nms_level <- strsplit(taxs, split = "|", fixed = TRUE) %>%
            lengths()
        indx <- purrr::map(indx, ~ .x & feature_nms_level == (i + 1))
        abd <- purrr::map_df(
            feature,
            ~ purrr::map_dbl(indx, function(x) sum(.x[x]))
        )

        feature <- rbind(feature, abd)
        row.names(feature) <- c(taxs, missing_nms_i)
    }

    otu_table(feature, taxa_are_rows = TRUE)
}

# check whether tax have level prefix, such as `p__`
check_tax_prefix <- function(taxa_nms) {
    prefix <- paste0(c("k", "p", "c", "o", "f", "g", "s"), "__")
    has_prefix <- purrr::map_lgl(prefix, 
        ~ any(grepl(.x, taxa_nms, fixed = TRUE))
    )

    any(has_prefix)
}
