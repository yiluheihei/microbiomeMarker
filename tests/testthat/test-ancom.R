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
  skip_on_bioc()
  ancom_res <- run_ancom(
    ecam, "delivery",
    test = "wilcox.test",
    p_adjust = "BH",
    W_cutoff = 0
  )
  curr_marker <- round_DF(marker_table(ancom_res))
  expect_known_output(
    head(curr_marker),
    test_path("out/test-ancom_marker.txt"),
    print = TRUE
  )
})
