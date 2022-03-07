data(enterotypes_arumugam)
enterotype <- phyloseq::subset_samples(
    enterotypes_arumugam,
    Enterotype %in% c("Enterotype 3", "Enterotype 2", "Enterotype 1")
)

tukey_res <- run_posthoc_test(enterotype, "Enterotype", method = "tukey")

test_that("etaseq effect size", {
    etasq <- calc_etasq(c(1, 2, 1.2, 3, 4, 1.4), c("a", "b", "c", "a", "b", "c"))
    expect_equal(signif(etasq, 3), 0.421)
})

test_that("test multiple group enterotype result", {

    # error group
    expect_error(
        run_test_multiple_groups(enterotype, "Entertype"),
        regexp = "`group` must in the field of sample meta data",
        fixed = TRUE
    )
})

test_that("test post hoc test result", {
    res_test <- tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]] %>% 
        data.frame
    expect_snapshot(print(res_test, digits = 5))
})

test_that("test visualization of post hoc test, p value significance level ", {
    expect_equal(
        pvalue2siglevel(c(0.05, 0.01, 0.0001, 0.06)),
        c("*", "**", "***", "NS.")
    )
})

test_that(
    "test visualization of posthoc test, data of signicance level annotation",
    {
        # single feature
        abd <- tukey_res@abundance
        group <- abd$group
        pht_df <- data.frame(tukey_res@result[["p__Bacteroidetes|g__Bacteroides"]])
        annotation_single <- get_sig_annotation_single(
            abd[["p__Bacteroidetes|g__Bacteroides"]],
            pht_df,
            group
        )
        annotation_single$y_position <- formatC(
            annotation_single$y_position,
            format = "g",
            digits = 5
        )
        expect_snapshot(annotation_single)

        # all features
        annotation_all <- get_sig_annotation(tukey_res)
        annotation_all$y_position <- formatC(
            annotation_all$y_position,
            format = "g",
            digits = 5
        )
        expect_snapshot(head(annotation_all))
    }
)
