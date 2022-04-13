#' check whether tax abundance table is summarized or not
#' @noRd
check_tax_summarize <- function(ps) {
    taxa <- row.names(otu_table(ps))
    # whether taxa is separated by `|`,
    # may be required to add extra separate strings in the future
    has_separate <- any(grepl("[|]", taxa))

    has_separate
}


# whether picrust functional profile
is_picrust2 <- function(ps) {
    ps_ranks <- rank_names(ps)
    if ("Picrust_trait" %in% ps_ranks) TRUE else FALSE
}

#' check whether all names of taxonomic ranks include in available_ranks
#' @noRd
check_rank_names <- function(ps) {
    ps_ranks <- rank_names(ps)
    if (is_picrust2(ps)) {
        picrust_rank <- c("Picrust_trait", "Picrust_description")
        diff_rank <- setdiff(ps_ranks, picrust_rank)
        if (length(diff_rank)) {
            stop("ranks of picrust2 functional profile must be one of ",
                 paste(picrust_rank, collapse = ", "))
        }
    } else {
        if (!all(ps_ranks %in% available_ranks)) {
            stop(
                "ranks of taxonimic profile must be one of ",
                paste(available_ranks, collapse = ", ")
            )
        }
    }

    invisible(ps)
}

#' only first letter in lower case
#' @noRd
upper_firstletter <- function(x) {
    paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
}

#' add prefix of taxonomic ranks for summarized data construct from original
#' lefse (galaxy lefse or python app) input, p__, k__
#' @param ps a [`phyloseq::phyloseq-class`] object
#' @param ranks character vector, prefix of ranks to add, e.g. "p", "c"
#' @importFrom phyloseq taxa_names<-
#' @noRd
add_prefix_summarized <- function(ps, ranks_prefix, sep = "|") {
    tax <- tax_table(ps)@.Data[, 1]
    tax_split <- strsplit(tax, split = sep, fixed = TRUE)

    if (max(lengths(tax_split)) != length(ranks_prefix)) {
        stop(
            "The length of `ranks_prefix` muste be",
            " equal to number of taxonomic ranks.",
            call. = FALSE
        )
    }

    # ensure the ranks_prefix is contained in available_ranks
    # and in descending order
    available_prefix <- get_available_prefix(available_ranks)
    if (!all(ranks_prefix %in% available_prefix)) {
        stop("all elements of ranks_prefix must be contained ", "
            in available_ranks")
    }
    ranks_prefix <- keep_prefix_desc(ranks_prefix, type = "ranks_prefix")

    tax_prefix <- purrr::map(
        tax_split,
        ~ paste0(ranks_prefix[seq_along(.x)], "__", .x) %>%
            paste0(collapse = sep)
    )
    tax_prefix <- do.call(rbind, tax_prefix)
    colnames(tax_prefix) <- paste0(ranks_prefix, collapse = sep)
    tax_table(ps) <- tax_table(tax_prefix)
    taxa_names(ps) <- tax_prefix[, 1]

    ps
}

# extract the first letter of taxonomic ranks as the prefixes of the ranks
get_available_prefix <- function(ranks) {
    substr(ranks, 1, 1) %>%
        tolower()
}

# keep prefix in descending order: "k" "p" "c" "o" "f" "g" "s"
keep_prefix_desc <- function(ranks_prefix, type = c("ranks", "ranks_prefix")) {
    type <- match.arg(type, choices = c("ranks", "ranks_prefix"))
    available_prefix <- get_available_prefix(available_ranks)
    idx_desc <- sort(match(ranks_prefix, available_prefix))

    if (type == "ranks") {
        return(available_ranks[idx_desc])
    } else {
        return(available_prefix[idx_desc])
    }
}

# check whether var in sample meta data, or raise an error
check_var_in_meta <- function(var, sample_meta) {
    stopifnot(inherits(sample_meta, "sample_data"))
    meta_nms <- names(sample_meta)
    if (!var %in% meta_nms) {
        stop(var, " must be one of variable of `sample_meta`", call. = FALSE)
    }
}

################################################################################
## preprocessing ps object
################################################################################

# preprocess of phyloseq object, including keep taxa in rows,
# filter taxa whose abundance is zero, fix duplicated tax
# filter samples whose library size is zero
#' @importFrom phyloseq prune_samples
preprocess_ps <- function(ps) {
    zero_sample <- check_samples(ps)
    if (!is.null(zero_sample)) {
        warning(
            "The library size of sample(s): ",
            paste(zero_sample, collapse = ", "),
            " is/are zero, and will be removed in the subsequent analysis."
        )

        keep <- setdiff(sample_names(ps), zero_sample)
        ps <- prune_samples(keep, ps)
    }

    # keep taxa in rows
    ps <- keep_taxa_in_rows(ps)
    # filter the taxa whose abundance is zero
    ps <- phyloseq_qc(ps)
    # fix duplicated tax
    ps <- fix_duplicate_tax(ps)

    ps
}

#' phyloseq quality control, remove otu/asv of which abundance is zero
#' @noRd
phyloseq_qc <- function(ps) {
    prune_taxa(taxa_sums(ps) > 0, ps)
}

#' Transpose the phyloseq object to ensure taxa are in rows
#' @param ps a [phyloseq::phyloseq-class] object
#' @importMethodsFrom phyloseq t
#' @keywords internal
#' @noRd
keep_taxa_in_rows <- function(ps) {
    if (!taxa_are_rows(ps)) {
        ps <- t(ps)
    }

    ps
}

# https://github.com/lch14forever/microbiomeVizb
# /94cbfe452a735aadf88733b27b8221a03f450a55/R/utils.R#L68-L86
#
#' Duplicated taxa: e.g. maybe multiple species (s__uncultured)
#' belong to different genera. append the upper level taxa to the taxa to
#' distinguish this duplicated taxa
#' @param ps [phyloseq::phyloseq-class] object or
#'   [phyloseq::taxonomyTable-class] object
#' @importFrom phyloseq tax_table<-
#' @keywords internal
#' @noRd
fix_duplicate_tax <- function(ps) {
    # convert na to Unknown first
    ps <- fix_na_tax(ps)

    tax <- tax_table(ps)
    if (ncol(tax) == 1) {
        return(ps)
    }

    for (i in 2:ncol(tax)) {
        tax_uniq <- unique(tax[, i])
        for (j in seq_along(tax_uniq)) {
            if (is.na(tax_uniq[j])) next
            ind <- which(tax[, i] == as.character(tax_uniq[j]))
            if (length(unique(tax[ind, i - 1])) > 1) {
                tax[ind, i] <- paste(tax[ind, i - 1], tax[ind, i], sep = "_")
            }
        }
    }

    tax_table(ps) <- tax

    ps
}

#' set NA (missing) tax to its prefix, e.g. s__ (or s__unknown?)
#' @keywords internal
#' @noRd
fix_na_tax <- function(ps) {
    tax <- as.data.frame(tax_table(ps))

    tax_fixed <- purrr::imap_dfc(
        tax,
        ~ ifelse(is.na(.x), get_prefix(.y), .x)
    ) %>%
        as.matrix()
    row.names(tax_fixed) <- taxa_names(ps)
    tax_table(ps) <- tax_fixed

    ps
}

# extract the prefix of taxonomic ranks
get_prefix <- function(ranks) {
    prefix <- substr(ranks, 1, 1) %>%
        tolower() %>%
        paste("__", sep = "")

    prefix
}


# `metagenomeSeq::cumNormStatFast()` requires counts of all samples at least
# have two non zero features. Thus, if there are samples with only one non-zer
# features, `cumNormStat()` is taken to compute the pth quantile.
# This function was used to select the function to calculate the quantile used
# for CSS norm factors estimation in metagenomeSeq.
select_quantile_func <- function(counts) {
    if (sum(colSums(counts > 0) > 1) < ncol(counts)) {
        fun_p <- metagenomeSeq::cumNormStat
    } else {
        fun_p <- metagenomeSeq::cumNormStatFast
    }

    fun_p
}


get_norm_method <- function(norm) {
    new_norm <- ifelse(
        is.numeric(norm),
        paste("per-sample normalized (sum of all taxa) to", norm),
        norm
    )

    new_norm
}


# check samples in ps, make sure at least one non zero features in a sample
check_samples <- function(ps) {
    if (!taxa_are_rows(ps)) {
        ps <- t(ps)
    }
    lib_size <- colSums(otu_table(ps))
    zero_ind <- which(lib_size == 0)

    if (length(zero_ind) == 0) {
        return(NULL)
    }

    return(sample_names(ps)[zero_ind])
}


# remove samples with missing values for any variable specified in the group
remove_na_samples <- function(ps, group) {
    groups <- sample_data(ps)[[group]]
    na_idx <- is.na(groups)
    if (all(!na_idx)) {
        return(ps)
    }

    sample_nms <- sample_names(ps)
    warning(
        "remove sample(s): ", paste(sample_nms[na_idx], collapse = ","),
        ", which with missing value in `", group, "`"
    )
    ps <- phyloseq::prune_samples(sample_nms[!na_idx], ps)

    ps
}


## calculate coef for edgeR, metagenomeSeq
# if contrast is a two length character, set the first element as the first level
# (reference group), the second element as the second level, return a single
# integer
#
# if contrast is null, return a integer vector (number of levels - 1)
check_contrast <- function(contrast) {
    if (!is.null(contrast)) {
        if (!is.character(contrast) || length(contrast) != 2) {
            stop("`contrast` must be a two length character", call. = FALSE)
        }
    }

    contrast
}

set_lvl <- function(groups, contrast) {
    if (!is.factor(groups)) {
        stop("`groups` must be a factor", call. = FALSE)
    }
    # this will will change the elements simultaneously
    # levels(groups) <- c(contrast, setdiff(levels(groups), contrast))
    groups <- factor(
        groups,
        levels =  c(contrast, setdiff(levels(groups), contrast))
    )

    groups
}

create_design <- function(groups, meta, confounders = character(0)) {
    if (inherits(meta, "sample_data")) {
        meta <- data.frame(meta)
    }
    model_data <- data.frame(group = groups)
    if (!length(confounders)) {
        design <- stats::model.matrix(~ group, data = model_data)
    } else {
        model_data[confounders] <- meta[confounders]
        design <- stats::model.matrix(
            formula(paste(
                "~ + ",
                paste(c(confounders, "group"), collapse = " + "))),
            data = model_data
        )
    }

    design
}

calc_coef <- function(groups, design, contrast = NULL) {
    contrast <- check_contrast(contrast)
    groups <- set_lvl(groups, contrast)
    lvl <- levels(groups)
    n_lvl <- length(lvl)
    n_design <- ncol(design)

    if (n_lvl < 2) {
        stop("Differential analysis requires at least two groups.")
    }

    if (n_lvl == 2) { # two groups
        if (!is.null(contrast)) {
            warning(
                "`contrast` is ignored, you do not need to set it",
                call. = FALSE
            )
        }
        coef <- n_design
    } else { # multiple groups
        if (!is.null(contrast)) {
            coef <- n_design - n_lvl + 2L
        } else {
            coef <- (n_design - n_lvl + 2L):n_design
        }
    }

    coef
}

# create_contrast <- function(groups, contrast = NULL) {
#     if (!is.factor(groups)) {
#         groups <- factor(groups)
#     }
#     lvl <- levels(groups)
#     n_lvl <- length(lvl)
#     if (n_lvl < 2) {
#         stop("Differential analysis requires at least two groups.")
#     }
#
#     if (n_lvl == 2) { # two groups
#         if (!is.null(contrast)) {
#             warning(
#                 "`contrast` is ignored, you do not need to set it",
#                 call. = FALSE
#             )
#         }
#         design <- rep(0, n_lvl)
#         design[1] <- -1
#         design[2] <- 1
#     } else { # multiple groups
#         if (!is.null(contrast)) {
#             if (!is.character(contrast) || length(contrast) != 2) {
#                 stop("`contrast` must be a two length character", call. = FALSE)
#             }
#
#             idx <- match(contrast, lvl, nomatch = 0L)
#             if (!all(idx)) {
#                 stop(
#                     "all elements of `contrast` must be contained in `groups`",
#                     call. = FALSE
#                 )
#             }
#             design <- rep(0, n_lvl)
#             design[idx[1]] <- -1
#             design[idx[2]] <- 1
#             design <- matrix(design)
#             row.names(design) <- lvl
#             colnames(design) <- paste0(contrast[2], "-", contrast[1])
#         } else {
#             design <- create_pairwise_contrast(lvl)
#         }
#     }
#
#     design
# }
#
# # create all pair-wise comparisons (contrasts) for anova-like test
# create_pairwise_contrast <- function(groups) {
#     groups <- factor(groups)
#     lvl <- levels(groups)
#     n <- length(lvl)
#
#     design <- matrix(0, n, choose(n, 2))
#     rownames(design) <- lvl
#     colnames(design) <- seq_len(choose(n, 2))
#     k <- 0
#     for (i in seq_len(n - 1)) {
#         for (j in (i + 1):n) {
#             k <- k + 1
#             design[j, k] <- 1
#             design[i, k] <- -1
#             colnames(design)[k] <- paste(lvl[j], "-", lvl[i], sep = "")
#         }
#     }
#     design
# }


# extract the specify rank of phyloseq object, return a phyloseq object
# with only one rank
extract_rank <- function(ps, taxa_rank) {
    ranks <- rank_names(ps)
    if (!taxa_rank %in% c("none", ranks)) {
        stop(
            "`taxa_rank` must be one of options: none, ",
            paste(rank_names(ps), collapse = ", "),
            call. = FALSE
        )
    }
    if (taxa_rank != "none") {
        ps <- aggregate_taxa(ps, taxa_rank)
        new_tax_table <- tax_table(ps)[, taxa_rank]
    } else {
        taxon <- taxa_names(ps)
        new_tax_table <- tax_table(matrix(taxon))
        colnames(new_tax_table) <- "otu"
        rownames(new_tax_table) <- taxon
    }

    tax_table(ps) <- new_tax_table

    # set the taxa names as the corresponding names
    if (taxa_rank != "none") {
        taxa_names(ps) <- new_tax_table[, 1]
    }

    ps
}

# only used for check the argument taxa_rank which is used to specify
# taxonomic rank to perform differential analysis on
check_taxa_rank <- function(ps, taxa_rank) {
    ranks <- rank_names(ps)
    all_taxa_rank <- c("all", "none", ranks)
    if (!taxa_rank %in% all_taxa_rank) {
        stop(
            "`taxa_rank` must be one of ",
            paste(all_taxa_rank, collapse = ", "),
            call. = FALSE
        )
    }

    invisible(ps)
}

# preprocess the ps according to para taxa_rank
pre_ps_taxa_rank <- function(ps, taxa_rank) {
    if (is_picrust2(ps)) {
        warning("para `taxa_rank` is not worked for picrust2 function profile",
                " and it will be ignored")
        return(ps)
    }

    ps <- check_taxa_rank(ps, taxa_rank)
    if (taxa_rank == "all") {
        ps_orig_summarized <- summarize_taxa(ps)
    } else if (taxa_rank == "none") {
        ps_orig_summarized <- extract_rank(ps, taxa_rank)
    } else {
        ps_orig_summarized <- aggregate_taxa(ps, taxa_rank) %>%
            extract_rank(taxa_rank)
    }

    return(ps_orig_summarized)
}

# return the marker_table, if no significant marker return NULL
return_marker <- function(sig_feature, all_feature) {
    if (nrow(sig_feature)) {
        row.names(sig_feature) <- paste0("marker", seq_len(nrow(sig_feature)))
        marker <- marker_table(sig_feature)

    } else {
        warning("No marker was identified", call. = FALSE)
        marker <- NULL
    }

    marker
}


# For multiple groups comparison of LRT test of DESeq2.
# Only fold changes of pair-wise comparisons are supported in DESse2.
# https://support.bioconductor.org/p/131272/#131450
# https://github.com/qiime2/q2-composition/blob/HEAD/q2_composition/_ancom.py
#'
#' Calculate effect size, mean differences for two groups, and f statistic of
#' one-way anova for multiple groups.
#' @noRd
#' @importFrom stats lm aov
calc_ef_md_f <- function(feature_abd, group) {
    group_n <- length(unique(group))
    if (group_n < 2) {
        stop("The number of group must be greater than 2", call. = FALSE)
    }

    if (group_n == 2) {
        ef <- abs(lm(feature_abd ~ group)$coefficients[2])
    } else if (group_n > 2) {
        # f statistic from aov
        ef <- summary(aov(feature_abd ~ group))[[1]]$`F value`[1]
    }

    ef
}

# create phyloseq from microbiomeMarker object,
# and keep only nodes correlated with significant features
create_ps_from_mm <- function(mm, only_marker = TRUE) {
    ot <- otu_table(mm)
    tt <- tax_table(mm)
    st <- sample_data(mm)
    mt <- marker_table(mm)
    sig_features <- mt$feature

    # extract all nodes correlated with the significant features
    # First, all parent nodes of marker
    down_nodes <- strsplit(sig_features, "|", fixed = TRUE) %>%
        purrr::map(~ purrr::map_chr(
            seq_along(.x), function(y) paste(.x[1:y], collapse = "|")))
    down_nodes <- unique(unlist(down_nodes))
    # Two, all children nodes of marker
    all_features <- tt@.Data[, 1]
    up_nodes <- purrr::map(sig_features,
               ~ all_features[grepl(.x, all_features, fixed = TRUE)])
    up_nodes <- unique(unlist(up_nodes))
    idx <- match(unique(c(down_nodes, up_nodes)), all_features)

    if (only_marker) {
        ot <- ot[idx, ]
        tt <- tt[idx, ]
    }
    ps <- phyloseq(ot, tt, st)

    ps
}

# check confounder
check_confounder <- function(ps, target_var, confounders = NULL) {
    meta <- sample_data(ps)
    vars <- names(meta)

    if (! target_var %in% vars) {
        stop(
            "the interested var `target_var` must be contained in the meta data",
            call. = FALSE
        )
    }

    other_vars <- setdiff(vars, target_var)

     if (is.null(confounders)) {
        confounders <- other_vars
        if (! length(confounders)) {
            stop("No confounding var in sample meta data")
        }
    } else {
        out_confounder <- setdiff(confounders, other_vars)
        if (length(out_confounder)) {
            stop("var(s) `", paste(out_confounder, collapse = "`, ` "),
                 "` not be contained in the sample meta data")
        }
    }

    confounders
}

# generate n spaces character
space <- function(n) {
    paste(rep(" ", each = n), collapse = "")
}
