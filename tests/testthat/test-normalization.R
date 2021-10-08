context("taxa abundance normalization")

ct <- as(otu_table(pediatric_ibd), "matrix")
gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans <- apply(ct, 1, gm_mean)

test_that("ensure the results are the same for object in different class ", {
    ot <- otu_table(enterotypes_arumugam)
    df <- as.data.frame(ot)
    mat <- as.matrix(df)


    # tss ---------------------------------------------------------------------
    # results are the same
    tss_ps <- norm_tss(enterotypes_arumugam)
    tss_ot <- norm_tss(ot)
    tss_df <- normalize(df, "TSS")
    tss_mat <- normalize(mat, "TSS")
    expect_equal(tss_ot, otu_table(tss_ps))
    expect_equal(as.data.frame(tss_ot), tss_df)
    expect_equal(as.matrix(tss_df), tss_mat)

    # no norm_factor
    expect_true(is.null(attr(tss_df, "norm_factor")))
    expect_true(is.null(attr(tss_mat, "norm_factor")))
    expect_true(is.null(attr(tss_ot, "norm_factor")))
    expect_true(is.null(attr(otu_table(tss_ps), "norm_factor")))
    expect_false("norm_factor" %in% names(sample_data(tss_ps)))


    # rarefy ------------------------------------------------------------------
    rarefy_ps <- norm_rarefy(enterotypes_arumugam, rng_seed = 2020)
    rarefy_ot <- norm_rarefy(ot, rng_seed = 2020)
    rarefy_df <- normalize(df, "rarefy", rng_seed = 2020)
    rarefy_mat <- normalize(mat, "rarefy", rng_seed = 2020)
    expect_equal(otu_table(rarefy_ps), rarefy_ot)
    expect_equal(as.data.frame(rarefy_ot), rarefy_df)
    expect_equal(as.matrix(rarefy_df), rarefy_mat)

    # no norm_factor
    expect_true(is.null(attr(rarefy_df, "norm_factor")))
    expect_true(is.null(attr(rarefy_mat, "norm_factor")))
    expect_true(is.null(attr(rarefy_ot, "norm_factor")))
    expect_true(is.null(attr(otu_table(rarefy_ps), "norm_factor")))
    expect_false("norm_factor" %in% names(sample_data(rarefy_ps)))


    # css ---------------------------------------------------------------------
    css_ps <- norm_css(enterotypes_arumugam)
    css_ot <- norm_css(ot)
    css_df <- normalize(df, "CSS")
    css_mat <- normalize(mat, "CSS")
    expect_identical(otu_table(css_ps), css_ot)
    expect_equivalent(as.data.frame(css_ot), css_df)
    expect_equivalent(as.matrix(css_df), css_mat)

    # norm factor
    css_nf_ot <- attr(css_ot, "norm_factor")
    css_nf_ps <- sample_data(css_ps)$norm_factor
    css_nf_df <- attr(css_df, "norm_factor")
    css_nf_mat <- attr(css_mat, "norm_factor")
    expect_identical(css_nf_ps, css_nf_ot)
    expect_identical(css_nf_ot, css_nf_df)
    expect_identical(css_nf_mat, css_nf_df)


    # rle ---------------------------------------------------------------------
    rle_ps <- norm_rle(enterotypes_arumugam)
    rle_ot <- norm_rle(ot)
    rle_df <- normalize(df, "RLE")
    rle_mat <- normalize(mat, "RLE")
    expect_identical(otu_table(rle_ps), rle_ot)
    expect_equivalent(as.data.frame(rle_ot), rle_df)
    expect_equivalent(as.matrix(rle_df), rle_mat)

    # norm factor
    rle_nf_ot <- attr(rle_ot, "norm_factor")
    rle_nf_ps <- sample_data(rle_ps)$norm_factor
    rle_nf_df <- attr(rle_df, "norm_factor")
    rle_nf_mat <- attr(rle_mat, "norm_factor")
    expect_identical(rle_nf_ps, rle_nf_ot)
    expect_identical(rle_nf_ot, rle_nf_df)
    expect_identical(rle_nf_mat, rle_nf_df)


    # TMM ---------------------------------------------------------------------
    tmm_ps <- norm_tmm(enterotypes_arumugam)
    tmm_ot <- norm_tmm(ot)
    tmm_df <- normalize(df, "TMM")
    tmm_mat <- normalize(mat, "TMM")
    expect_identical(otu_table(tmm_ps), tmm_ot)
    expect_equivalent(as.data.frame(tmm_ot), tmm_df)
    expect_equivalent(as.matrix(tmm_df), tmm_mat)

    # norm factor
    tmm_nf_ot <- attr(tmm_ot, "norm_factor")
    tmm_nf_ps <- sample_data(tmm_ps)$norm_factor
    tmm_nf_df <- attr(tmm_df, "norm_factor")
    tmm_nf_mat <- attr(tmm_mat, "norm_factor")
    expect_identical(tmm_nf_ps, tmm_nf_ot)
    expect_identical(tmm_nf_ot, tmm_nf_df)
    expect_identical(tmm_nf_mat, tmm_nf_df)


    # clr ---------------------------------------------------------------------
    clr_ps <- norm_clr(enterotypes_arumugam)
    clr_ot <- norm_clr(ot)
    clr_df <- normalize(df, "CLR")
    clr_mat <- normalize(mat, "CLR")
    expect_identical(otu_table(clr_ps), clr_ot)
    expect_identical(as.data.frame(clr_ot), clr_df)
    expect_identical(as.matrix(clr_df), clr_mat)

    # no norm_factor
    expect_true(is.null(attr(clr_df, "norm_factor")))
    expect_true(is.null(attr(clr_mat, "norm_factor")))
    expect_true(is.null(attr(clr_ot, "norm_factor")))
    expect_true(is.null(attr(otu_table(clr_ps), "norm_factor")))
    expect_false("norm_factor" %in% names(sample_data(clr_ps)))


    # cpm -------------------------------------------------------------------
    value_ps <- norm_cpm(enterotypes_arumugam)
    value_ot <- norm_cpm(ot)
    value_df <- normalize(df, "CPM")
    value_mat <- normalize(mat, "CPM")
    expect_identical(otu_table(value_ps), value_ot)
    expect_identical(as.data.frame(value_ot), value_df)
    expect_identical(as.matrix(value_df), value_mat)

    # no norm_factor
    expect_true(is.null(attr(value_df, "norm_factor")))
    expect_true(is.null(attr(value_mat, "norm_factor")))
    expect_true(is.null(attr(value_ot, "norm_factor")))
    expect_true(is.null(attr(otu_table(value_ps), "norm_factor")))
    expect_false("norm_factor" %in% names(sample_data(value_ps)))
})


test_that(paste0(
    "`geoMeans` and `type = 'poscounts'` in ",
    "`estimateSizeFactorsForMatrix` are equal"
), {
    expect_equal(
        estimateSizeFactorsForMatrix(ct, geoMeans = geoMeans),
        estimateSizeFactorsForMatrix(ct, type = "poscounts")
    )
})

test_that(paste0(
    "the size factors from `estimateSizeFactorForMatrix()` ",
    "and `DESeq2::estimateSizeFactors` are equal"
), {
    suppressWarnings(dds <- phyloseq2DESeq2(pediatric_ibd, ~Class))
    sf1 <- DESeq2::estimateSizeFactors(dds, type = "poscounts") %>%
        DESeq2::sizeFactors()
    sf2 <- DESeq2::estimateSizeFactors(dds, geoMeans = geoMeans) %>%
        DESeq2::sizeFactors()
    sf3 <- estimateSizeFactorsForMatrix(ct, geoMeans = geoMeans)
    sf4 <- estimateSizeFactorsForMatrix(ct, type = "poscounts")

    expect_identical(sf1, sf2)
    expect_equal(sf3, sf4)
    expect_equal(sf1, sf3)
})


# clr normalization
test_that("gm_mean works well with NA values", {
    expect_equal(gm_mean(c(1, NA, 2, 3)), gm_mean(c(1, 0, 2, 3)))
    expect_equal(gm_mean(c(1, -1, 2, 3), na.rm = FALSE), gm_mean(c(1, 0, 2, 3)))
})

test_that("trans_clr works well with NA/infinite values", {
    expect_equal(trans_clr(c(1, NA, 2, 3)), trans_clr(c(1, 0, 2, 3)))
    # Actually，Inf will unlikely to appear
    expect_equal(trans_clr(c(1, Inf, 2, 3)), c(0, 0, 0, 0))
})
