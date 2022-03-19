#' Perform differential analysis using ANCOM
#'
#' Perform significant test by comparing the pairwise log ratios between all
#' features.
#'
#' @param ps a \code{\link[phyloseq]{phyloseq-class}} object.
#' @param group character, the variable to set the group.
#' @param confounders character vector, the confounding variables to be adjusted.
#'   default `character(0)`, indicating no confounding variable.
#' @param taxa_rank character to specify taxonomic rank to perform
#'   differential analysis on. Should be one of
#'   `phyloseq::rank_names(phyloseq)`, or "all" means to summarize the taxa by
#'   the top taxa ranks (`summarize_taxa(ps, level = rank_names(ps)[1])`), or
#'   "none" means perform differential analysis on the original taxa
#'   (`taxa_names(phyloseq)`, e.g., OTU or ASV).
#' @param transform character, the methods used to transform the microbial
#'   abundance. See [`transform_abundances()`] for more details. The
#'   options include:
#'   * "identity", return the original data without any transformation.
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
#' @param norm_para  named `list`. other arguments passed to specific
#'   normalization methods.  Most users will not need to pass any additional
#'   arguments here.
#' @param p_adjust method for multiple test correction, default `none`,
#' for more details see [stats::p.adjust].
#' @param pvalue_cutoff significance level for each of the statistical tests,
#'   default 0.05.
#' @param W_cutoff lower bound for the proportion for the W-statistic, default
#'   0.7.
#'
#' @details
#' In an experiment with only two treatments, this tests the following
#' hypothesis for feature \eqn{i}:
#'
#' \deqn{H_{0i}: E(log(\mu_i^1)) =  E(log(\mu_i^2))}
#'
#' where \eqn{\mu_i^1} and \eqn{\mu_i^2} are the mean abundances for feature
#' \eqn{i} in the two groups.
#'
#' The developers of this method recommend the following significance tests
#' if there are 2 groups, use non-parametric Wilcoxon rank sum test
#' [`stats::wilcox.test()`]. If there are more than 2 groups, use nonparametric
#' [`stats::kruskal.test()`] or one-way ANOVA [`stats::aov()`].
#'
#' @return a [microbiomeMarker-class] object, in which the `slot` of
#' `marker_table` contains four variables:
#' * `feature`, significantly different features.
#' * `enrich_group`, the class of the differential features enriched.
#' *  `effect_size`, differential means for two groups, or F statistic for more
#'   than two groups.
#' * `W`, the W-statistic, number of features that a single feature is tested
#'   to be significantly different against.
#'
#' @references Mandal et al. "Analysis of composition of microbiomes: a novel
#' method for studying microbial composition", Microbial Ecology in Health
#' & Disease, (2015), 26.
#' @author Huang Lin, Yang Cao
#' @export
#' @examples
#' \donttest{
#' data(enterotypes_arumugam)
#' ps <- phyloseq::subset_samples(
#'     enterotypes_arumugam,
#'     Enterotype %in% c("Enterotype 3", "Enterotype 2")
#' )
#' run_ancom(ps, group = "Enterotype")
#' }
run_ancom <- function(ps,
    group,
    confounders = character(0),
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "TSS",
    norm_para = list(),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    pvalue_cutoff = 0.05,
    W_cutoff = 0.75) {
    stopifnot(inherits(ps, "phyloseq"))
    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )
    
    ps <- check_rank_names(ps) %>% 
        check_taxa_rank( taxa_rank)
    
    if (length(confounders)) {
        confounders <- check_confounder(ps, group, confounders)
    }

    # check whether group is valid, write a function
    meta <- sample_data(ps)
    meta_nms <- names(meta)
    groups <- meta[[group]]
    groups <- make.names(groups)
   
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    sample_data(ps)[[group]] <- groups
    lvl <- levels(groups)
    n_lvl <- length(lvl)
    
    if (!length(confounders)) {
       tfun <- ifelse(n_lvl > 2, stats::kruskal.test, stats::wilcox.test)
       fml <- paste("x ~ ", group)
    } else {
        tfun <- stats::aov
        fml <- paste("x ~ ", group, "+",  paste(confounders, collapse = " + "))
    }

    # preprocess phyloseq object
    ps <- preprocess_ps(ps)
    ps <- transform_abundances(ps, transform = transform)

    # normalize the data
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(normalize, norm_para)
    ps_summarized <- pre_ps_taxa_rank(ps_normed, taxa_rank)
    feature_table <- abundances(ps_summarized, norm = TRUE)

    # effect size: CLR mean_difference or aov f statistic
    feature_table_clr <- norm_clr(
        otu_table(feature_table, taxa_are_rows = TRUE)
    )
    feature_table_clr <- data.frame(t(feature_table_clr))
    ef <- vapply(
        feature_table_clr,
        calc_ef_md_f,
        FUN.VALUE = 0.0,
        group = groups
    )

    # enrich_group
    group_enriched <- vapply(
        feature_table_clr,
        get_ancom_enrich_group,
        FUN.VALUE = character(1),
        group = groups
    )

    # ANCOM requires log transformation
    feature_table <- log(as.matrix(feature_table) + 1)
    n_taxa <- nrow(feature_table)
    taxa_id <- row.names(feature_table)
    n_samp <- ncol(feature_table)

    # Calculate the p-value for each pairwise comparison of taxa.
    # para group is just for the main var in the formula
    test_var_dat <- data.frame(groups)
    names(test_var_dat) <- group
    if (length(confounders)) {
        test_var_dat[[confounders]] <- meta[[confounders]]
    }
    p <- calc_ancom_pmat(
        feature_table, 
        test_var_dat,
        tfun, 
        fml
    )

    # Multiple comparisons correction.
    p_adjusted <- vapply(
        data.frame(p),
        p.adjust,
        FUN.VALUE = numeric(n_taxa),
        method = p_adjust
    )

    # Calculate the W statistic of ANCOM.
    # For each taxon, count the number of q-values < pvalue_cutoff.
    W <- apply(p_adjusted, 2, function(x) sum(x < pvalue_cutoff))

    # Organize outputs
    out_comp <- data.frame(
        feature = taxa_id,
        enrich_group = group_enriched,
        ef = ef,
        W = W,
        row.names = NULL,
        check.names = FALSE
    )
    # Declare a taxon to be differentially abundant based on the quantile of W
    # statistic. We perform (n_taxa - 1) hypothesis testings on each taxon, so
    # the maximum number of rejections is (n_taxa - 1).
    sig_out <- out_comp[out_comp$W > W_cutoff * (n_taxa - 1), ]
    if (n_lvl == 2) {
        names(sig_out)[3] <- "ef_CLR_diff_mean"
    } else {
        names(sig_out)[3] <- "ef_CLR_F_statistic"
    }

    marker <- return_marker(sig_out, out_comp)
    tax <- matrix(taxa_id) %>%
        tax_table()
    row.names(tax) <- row.names(feature_table)

    mm <- microbiomeMarker(
        marker_table = marker,
        norm_method = get_norm_method(norm),
        diff_method = "ANCOM",
        otu_table = otu_table(feature_table, taxa_are_rows = TRUE),
        sam_data = sample_data(ps_normed),
        tax_table = tax
    )

    mm
}

#' Calculates pairwise pvalues between all features
#' @param feature_table matrix-like, logged feature table.
#' @param test_var_dat data.frame, variables data (sample meta data)
#' @param test  character, the test to determine the p value of log ratio,
#'   one of "aov", "wilcox.test",  "kruskal.test".
#' @param ... extra arguments passed to the test.
#' @references 
#' github/biocore/scikit-bio/blob/master/skbio/stats/composition.py#L811
#' @noRd
calc_ancom_pmat <- function(feature_table, test_var_dat, test, fml) {

    taxas <- row.names(feature_table)
    feature_table <- data.frame(t(feature_table))
    taxa_n <- ncol(feature_table)
    p <- matrix(NA, nrow = taxa_n, ncol = taxa_n)
    row.names(p) <- taxas
    colnames(p) <- taxas


    for (i in seq_len(taxa_n - 1)) {
        new_table <- -(feature_table[(i + 1):taxa_n] - feature_table[[i]])
        p[-(seq_len(i)), i] <- vapply(
            new_table,
            calc_ancom_p,
            FUN.VALUE = numeric(1),
            test_var_dat = test_var_dat, test = test, fml = fml
        )
    }

    # Complete the p-value matrix.
    # What we got from above iterations is a lower triangle matrix of p-values.
    p[upper.tri(p)] <- t(p)[upper.tri(p)]
    diag(p) <- 1 # let p-values on diagonal equal to 1
    p[is.na(p)] <- 1 # let p-values of NA equal to 1

    p
}

#' calculate the p value of a pair-wise log ratio
#' @param log_ratio  a numeric vector, a pair-wise log ratio.
#' @param classes character vector, the same length with `log_ratio`.
#' @param test  character, the test to dtermine the p value of log ratio,
#'   one of "aov", "wilcox.test",  "kruskal.test".
#' @param ... extra arguments passed to the test.
#' @noRd
calc_ancom_p <- function(log_ratio, test_var_dat, test, fml) {
    # fist var is the target var (main var)
    group <- names(test_var_dat)[1]
    test_dat <- cbind(x = log_ratio, test_var_dat)
    fml <- stats::formula(fml)
    if (identical(test, stats::aov)) {
        fit = test(fml, 
                   data = test_dat, 
                   na.action = na.omit)
        p = summary(fit)[[1]][group, "Pr(>F)"]
    } else {
        suppressWarnings(p <- test(fml, data = test_dat)$p.value)
    }
   
    p
}


#' Identify structural zeros
#' from "FrederickHuangLin/ANCOMBC/R/get_struc_zero.R"
#'
#' @author Huang Lin, Yang Cao
#' @noRd
get_struc_zero <- function(ps, group, neg_lb) {
    stopifnot(inherits(ps, "phyloseq"))
    stopifnot(is.logical(neg_lb))
    stopifnot(length(group) == 1 & is.character(group))

    meta_tab <- sample_data(ps)
    check_var_in_meta(group, meta_tab)
    groups <- factor(meta_tab[[group]])

    feature_tab <- as(otu_table(ps), "matrix")
    present_tab <- feature_tab
    present_tab[is.na(present_tab)] <- 0
    present_tab[present_tab != 0] <- 1
    n_taxa <- nrow(feature_tab)
    n_group <- nlevels(groups)

    p_hat <- matrix(NA, nrow = n_taxa, ncol = n_group)
    rownames(p_hat) <- rownames(feature_tab)
    colnames(p_hat) <- levels(groups)
    samp_size <- p_hat

    for (i in seq_len(n_taxa)) {
        p_hat[i, ] <- tapply(
            present_tab[i, ],
            groups,
            function(x) mean(x, na.rm = TRUE)
        )
        samp_size[i, ] <- tapply(
            feature_tab[i, ],
            groups,
            function(x) length(x[!is.na(x)])
        )
    }

    p_hat_lo <- p_hat - 1.96 * sqrt(p_hat * (1 - p_hat) / samp_size)

    zero_ind <- p_hat == 0

    if (neg_lb) {
        zero_ind[p_hat_lo <= 0] <- TRUE
    }
    colnames(zero_ind) <- paste0(
        "structural_zero (", group, " = ", colnames(zero_ind), ")"
    )

    data.frame(zero_ind)
}

#' enrich group for ancom, rewrite this function in the later
#' split get_feature_enrich_group into two funcitons: enrich_group and
#' log max mean
#' @noRd
get_ancom_enrich_group <- function(feature_abd, group) {
    abd_split <- split(feature_abd, group)
    abd_mean_group <- vapply(abd_split, mean, FUN.VALUE = 0.0)
    enrich_group <- names(abd_split)[which.max(abd_mean_group)]

    enrich_group
}

#' preprocess feature data using methods of ANCOM-II
#' @noRd
#' @importFrom stats dnorm lm na.omit quantile residuals sd
preprocess_ancom <- function(feature_table,
    meta_data,
    sample_var,
    lib_cut,
    neg_lb,
    group = NULL,
    out_cut = 0.05,
    zero_cut = 0.90) {
    
    feature_table <- data.frame(feature_table, check.names = FALSE)
    meta_data <- data.frame(meta_data, check.names = FALSE)
    # Drop unused levels
    meta_data[] <- lapply(
        meta_data,
        function(x) if (is.factor(x)) factor(x) else x
    )
    # Match sample IDs between metadata and feature table
    sample_ID <- intersect(meta_data[, sample_var], colnames(feature_table))
    feature_table <- feature_table[, sample_ID]
    meta_data <- meta_data[match(sample_ID, meta_data[, sample_var]), ]

    # 1. Identify outliers within each taxon
    if (!is.null(group)) {
        groups <- meta_data[, group]
        z <- feature_table + 1 # Add pseudo-count (1)
        f <- log(z)
        f[f == 0] <- NA
        f <- colMeans(f, na.rm = TRUE)
        f_fit <- lm(f ~ groups)
        e <- rep(0, length(f))
        e[!is.na(groups)] <- residuals(f_fit)
        y <- t(t(z) - e)

        outlier_check <- function(x) {
            # Fitting the mixture model using the algorithm of Peddada, S. Das,
            # and JT Gene Hwang (2002)
            mu1 <- quantile(x, 0.25, na.rm = TRUE)
            mu2 <- quantile(x, 0.75, na.rm = TRUE)
            sigma1 <- quantile(x, 0.75, na.rm = TRUE) -
                quantile(x, 0.25, na.rm = TRUE)
            sigma2 <- sigma1
            pi <- 0.75
            n <- length(x)
            epsilon <- 100
            tol <- 1e-5
            score <- pi * dnorm(x, mean = mu1, sd = sigma1) /
                ((1 - pi) * dnorm(x, mean = mu2, sd = sigma2))
            while (epsilon > tol) {
                grp1_ind <- (score >= 1)
                mu1_new <- mean(x[grp1_ind])
                mu2_new <- mean(x[!grp1_ind])
                sigma1_new <- sd(x[grp1_ind])
                if (is.na(sigma1_new)) sigma1_new <- 0
                sigma2_new <- sd(x[!grp1_ind])
                if (is.na(sigma2_new)) sigma2_new <- 0
                pi_new <- sum(grp1_ind) / n

                para <- c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
                if (any(is.na(para))) break

                score <- pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new) /
                    ((1 - pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

                epsilon <- sqrt(
                    (mu1 - mu1_new)^2 +
                        (mu2 - mu2_new)^2 +
                        (sigma1 - sigma1_new)^2 +
                        (sigma2 - sigma2_new)^2 +
                        (pi - pi_new)^2
                )
                mu1 <- mu1_new
                mu2 <- mu2_new
                sigma1 <- sigma1_new
                sigma2 <- sigma2_new
                pi <- pi_new
            }

            if (mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2) {
                if (pi < out_cut) {
                    out_ind <- grp1_ind
                } else if (pi > 1 - out_cut) {
                    out_ind <- (!grp1_ind)
                } else {
                    out_ind <- rep(FALSE, n)
                }
            } else {
                out_ind <- rep(FALSE, n)
            }
            return(out_ind)
        }
        out_ind <- matrix(
            FALSE,
            nrow = nrow(feature_table),
            ncol = ncol(feature_table)
        )
        out_ind[, !is.na(groups)] <- t(apply(
            y, 1,
            function(i) {
                unlist(tapply(i, groups, function(j) outlier_check(j)))
            }
        ))

        feature_table[out_ind] <- NA
    }

    # 2. Discard taxa with zeros  >=  zero_cut
    zero_prop <- apply(
        feature_table, 1,
        function(x) sum(x == 0, na.rm = TRUE) / length(x[!is.na(x)])
    )
    taxa_del <- which(zero_prop >= zero_cut)
    if (length(taxa_del) > 0) {
        feature_table <- feature_table[-taxa_del, ]
    }

    # 3. Discard samples with library size < lib_cut
    lib_size <- colSums(feature_table, na.rm = TRUE)
    if (any(lib_size < lib_cut)) {
        subj_del <- which(lib_size < lib_cut)
        feature_table <- feature_table[, -subj_del]
        meta_data <- meta_data[-subj_del, ]
    }

    # 4. Identify taxa with structure zeros
    if (!is.null(group)) {
        groups <- factor(meta_data[, group])
        present_table <- as.matrix(feature_table)
        present_table[is.na(present_table)] <- 0
        present_table[present_table != 0] <- 1

        p_hat <- t(apply(
            present_table, 1,
            function(x) {
                unlist(tapply(x, groups, function(y) mean(y, na.rm = TRUE)))
            }
        ))
        samp_size <- t(apply(
            feature_table, 1,
            function(x) {
                unlist(tapply(x, groups, function(y) length(y[!is.na(y)])))
            }
        ))
        p_hat_lo <- p_hat - 1.96 * sqrt(p_hat * (1 - p_hat) / samp_size)

        struc_zero <- (p_hat == 0) * 1
        # Whether we need to classify a taxon into structural zero by its
        # negative lower bound?
        if (neg_lb) struc_zero[p_hat_lo <= 0] <- 1

        # Entries considered to be structural zeros are set to be 0s
        struc_ind <- struc_zero[, groups]
        feature_table <- feature_table * (1 - struc_ind)

        colnames(struc_zero) <- paste0(
            "structural_zero (",
            colnames(struc_zero),
            ")"
        )
    } else {
        struc_zero <- NULL
    }

    # 5. Return results
    res <- list(
        feature_table = feature_table,
        meta_data = meta_data,
        structure_zeros = struc_zero
    )

    res
}
