test_that("convert mc instances", {
    instance <- list(
        sample1 = data.frame(inst1 = runif(10), inst2 = runif(10)),
        sample2 = data.frame(inst1 = runif(10), inst2 = runif(10)),
        sample3 = data.frame(inst1 = runif(10), inst2 = runif(10))
    )
    instance_converted <- convert_instance(instance, 2)

    expect_identical(length(instance_converted), 2L)
    expect_identical(names(instance_converted[[1]]), paste0("sample", 1:3))
    expect_equal(
        instance_converted[[1]],
        tibble::tibble(
            sample1 = instance[[1]][[1]],
            sample2 = instance[[2]][[1]],
            sample3 = instance[[3]][[1]]
        )
    )
})


# fast test
data(kostic_crc)
x <- as(phyloseq::otu_table(kostic_crc), "matrix")[1:20, 1:20]
groups <- phyloseq::sample_data(kostic_crc)[["DIAGNOSIS"]][1:20]
idx1 <- groups == "Tumor"
idx2 <- groups == "Healthy"
x_clr <- suppressWarnings(ALDEx2::aldex.clr(x, groups, mc.samples = 12))
instance <- convert_instance(x_clr@analysisData, 12)[[1]]
# keep the same number of samples for two groups
idx2 <- which(idx2)[1:4]
instance <- instance[c(which(idx1), idx2)]
new_groups <- groups[c(which(idx1), idx2)]

test_that("t_fast equal to t.test and ALDEx2:::t.fast", {
    # unpaired
    t_res_up <-  apply(
        instance, 1,
        function(x) t.test(x[1:4], x[5:8], alternative = "greater")$p.value
    )
    expect_equal(t_res_up, t_fast(instance, new_groups))
    expect_equal(
       t_res_up,
       ALDEx2:::t.fast(instance, new_groups, paired = FALSE)$p
    )

    # paired
    t_res_p <- apply(
        instance, 1,
        function(x) t.test(x[1:4], x[5:8],
                           paired = TRUE,
                           alternative = "greater")$p.value
    )
    expect_equal(
        t_res_p,
        t_fast(instance, new_groups, paired = TRUE)
    )
    expect_equal(
        t_res_p,
        ALDEx2:::t.fast(instance, new_groups, paired = TRUE)$p
    )

})

test_that("wilcox_fast equal to wilcox.test", {
    # unpaired
    wilcox_res_p <- apply(
        instance, 1,
        function(x) wilcox.test(x[1:4], x[5:8],
                                alternative = "greater")$p.value
    )
    expect_equal(
        wilcox_res_p,
        wilcox_fast(instance, new_groups)
    )
    expect_equal(
        wilcox_res_p,
        ALDEx2:::wilcox.fast(instance, new_groups, paired = FALSE))

    # paired
    wilcox_res_p <- apply(
        instance, 1,
        function(x) wilcox.test(x[1:4], x[5:8],
                                paired = TRUE,
                                alternative = "greater")$p.value
    )
    expect_equal(
        wilcox_res_p,
        wilcox_fast(instance, new_groups, paired = TRUE)
    )
    expect_equal(
        wilcox_res_p,
        ALDEx2:::wilcox.fast(instance, new_groups, paired = TRUE)
    )

    expect_equal(
        apply(
            instance, 1, function(x) {
                wilcox.test(x[1:4], x[5:8],
                            alternative = "greater",
                            paired = TRUE, exact = TRUE)$p.value
            }
        ),
        wilcox_fast(instance, new_groups, paired = TRUE)
    )
})

test_that("aldex_t", {
    aldex_test_out <- ALDEx2::aldex.ttest(x_clr)
    t_out <- aldex_t(x_clr, groups, method = "t.test", mc_samples = 12)
    expect_equal(t_out$pvalue, aldex_test_out$we.ep)

    wilcox_out <- aldex_t(x_clr, groups, method = "wilcox.test", mc_samples = 12)
    expect_equal(wilcox_out$pvalue, aldex_test_out$wi.ep)
})

test_that("aldex_kw", {
    aldex_kw_out <- ALDEx2::aldex.kw(x_clr)
    kw_out <- aldex_kw(x_clr, groups, method = "kruskal", mc_samples = 12)
    expect_equal(kw_out$pvalue, aldex_kw_out$kw.ep)

    glm_out <- aldex_kw(x_clr, groups, method = "glm_anova", mc_samples = 12)
    expect_equal(glm_out$pvalue, aldex_kw_out$glm.ep)
})
