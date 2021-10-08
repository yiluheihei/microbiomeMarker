context("bar plot of effect size")

test_that("feature label in bar plot", {
    feature <- "Bacteria|Bacteroidetes|Bacteroidia|Bacteroidales|Bacteroidaceae"
    short_feature <- "Bacteria|Bacteroidetes"

    expect_equal(get_feature_label(feature, 1), "Bacteroidaceae")
    expect_equal(get_feature_label(feature, 2), "Bacteroidales|Bacteroidaceae")
    expect_equal(
        get_feature_label(feature, 7),
        "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
    )
    expect_equal(
        get_feature_label(feature, 0),
        "Bacteria|Bacteroidetes|Bacte..Bacteroidales|Bacteroidaceae"
    )
    expect_equal(get_feature_label(short_feature, 0), "Bacteria|Bacteroidetes")

    # replace "Unknown" in the species level as "sp."
    expect_equal(
        replace_unknown_species("s__Leucobacter_Unknown"),
        "s__Leucobacter_sp."
    )
    expect_equal(
        replace_unknown_species("g__abcd|s__Leucobacter_Unknown"),
        "g__abcd|s__Leucobacter_sp."
    )
    expect_equal(
        replace_unknown_species("g__abcd"), "g__abcd"
    )
    expect_equal(
        replace_unknown_species("g__abcd|s__Leucobacter_sp."),
        "g__abcd|s__Leucobacter_sp."
    )
})

mm <- microbiomeMarker(
    marker_table = marker_table(data.frame(
        feature = c("speciesA", "speciesB"),
        enrich_group = c("groupA", "groupB"),
        ef_logFC = c(-2, 2),
        pvalue = c(0.01, 0.01),
        padj = c(0.01, 0.01),
        row.names = c("marker1", "marker2")
    )),
    otu_table = otu_table(matrix(
        c(4, 1, 1, 4),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("speciesA", "speciesB"), c("sample1", "sample2"))
    ),
    taxa_are_rows = TRUE
    ),
    tax_table = tax_table(matrix(
        c("speciesA", "speciesB"),
        nrow = 2,
        dimnames = list(c("speciesA", "speciesB"), "Species")
    ))
)


test_that("label of x, and effect size in descending order", {
    # logFC, such as edgeR, DESeq2 for two groups comparison
    p_logfc <- plot_ef_bar(mm)
    expect_identical(p_logfc$labels$x, "log2 Fold Change")

    # lefse - lda
    mt <- marker_table(data.frame(
        feature = c("speciesA", "speciesB"),
        enrich_group = c("groupA", "groupA"),
        ef_lda = c(2, 3),
        pvalue = c(0.01, 0.01),
        padj = c(0.01, 0.01),
        row.names = c("marker1", "marker2")
    ))
    marker_table(mm) <- mt
    p_lda <- plot_ef_bar(mm)
    expect_identical(p_lda$labels$x, "LDA score (log10)")
    # descending order
    ef <- p_lda$data$effect_size
    expect_true(all(diff(ef) >= 0))

    # two groups test, diff_mean
    names(mt)[3] <- "ef_diff_mean"
    marker_table(mm) <- mt
    p_diff_mean <- plot_ef_bar(mm)
    expect_identical(p_diff_mean$labels$x, "Differential means")

    # multiple group
    names(mt)[3] <- "ef_eta_squared"
    marker_table(mm) <- mt
    p_eta_squared <- plot_ef_bar(mm)
    expect_identical(p_eta_squared$labels$x, "Eta squared")

    # CLR diff mean
    names(mt)[3] <- "ef_CLR_diff_mean"
    marker_table(mm) <- mt
    p_clr_diff <- plot_ef_bar(mm)
    expect_identical(p_clr_diff$labels$x, "CLR differential means")

    # CLR F statistic
    names(mt)[3] <- "ef_CLR_F_statistic"
    marker_table(mm) <- mt
    p_clr_f <- plot_ef_bar(mm)
    expect_identical(p_clr_f$labels$x, "CLR F statistic")

    # W statistic
    names(mt)[3] <- "ef_W"
    marker_table(mm) <- mt
    p_w <- plot_ef_bar(mm)
    expect_identical(p_w$labels$x, "W")

    # importance
    names(mt)[3] <- "ef_imp"
    marker_table(mm) <- mt
    p_imp <- plot_ef_bar(mm)
    expect_identical(p_imp$labels$x, "Importance score")

    # likelihood ratio statistic
    names(mt)[3] <- "ef_LR"
    marker_table(mm) <- mt
    p_lr <- plot_ef_bar(mm)
    expect_identical(p_lr$labels$x, "Likelihood ratio statistic")

    # F statistic
    names(mt)[3] <- "ef_F"
    marker_table(mm) <- mt
    p_f <- plot_ef_bar(mm)
    expect_identical(p_f$labels$x, "F statistic")
})
