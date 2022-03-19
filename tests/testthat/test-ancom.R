if (FALSE) {
    zero_neg_lb <- ANCOMBC:::get_struc_zero(
        data.frame(otu_table(ecam)),
        as(sample_data(ecam), "matrix"),
        "delivery",
        neg_lb = TRUE
    )
    write.csv(
        data.frame(zero_neg_lb),
        file = test_path("data", "ancom-zero_neg_lb.csv")
    )
    zero <- ANCOMBC:::get_struc_zero(
        data.frame(otu_table(ecam)),
        as(sample_data(ecam), "matrix"),
        "delivery",
        neg_lb = FALSE
    )
    write.csv(
        data.frame(zero),
        file = test_path("data", "ancom-zero.csv")
    )
}

data(ecam)
test_that("identify structural zeros", {
    expect_identical(
        get_struc_zero(ecam, "delivery", TRUE),
        read.csv(test_path("data/ancom-zero_neg_lb.csv"), row.names = 1)
    )

    expect_identical(
        get_struc_zero(ecam, "delivery", FALSE),
        read.csv(test_path("data/ancom-zero.csv"), row.names = 1)
    )
})

test_that("ancom result", {
    ancom_res <- run_ancom(
        ecam, "delivery",
        p_adjust = "BH",
        W_cutoff = 0,
        taxa_rank = "Class"
    )
    curr_marker <- marker_table(ancom_res)
    
    expect_snapshot(print(head(curr_marker), digits = 5))
})
