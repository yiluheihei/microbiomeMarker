test_that("comparing da methods", {
    data(ecam)
    expect_error(compare_DA(ecam, "delivery", methods = c("test", "ancombc")),
                 "methods test not available")
})

test_that("flatten args",{
    expect_error(generate_compare_args("lefse", list(a=1)),
                 "does not match DA methods")

    args <-  list(lefse = list(list(norm = "CPM"), list(norm = "TSS")))
    out <- generate_compare_args("lefse", args)
    expect_identical(out$methods, c("lefse_1", "lefse_2"))
    expect_identical(out$args,
                     list(lefse_1 = list(norm = "CPM"),
                          lefse_2 = list(norm = "TSS")))
    expect_identical(out$methods, names(out$args))

    out <- generate_compare_args(c("lefse", "ancom"), args)
    expect_identical(out$methods, c("lefse_1", "lefse_2", "ancom"))
    expect_identical(out$args,
                     list(lefse_1 = list(norm = "CPM"),
                          lefse_2 = list(norm = "TSS"),
                          ancom = list()))


    args <-  list(lefse = list(list(norm = "CPM", pvalue_cutoff = 0.05),
                               list(norm = "TSS", pvalue_cutoff = 0.05)),
                  ancom = list(norm = "CPM", W_cutoff = 0.75))
    out <- generate_compare_args(c("lefse", "ancom"), args)
    expect_identical(out$methods, c("lefse_1", "lefse_2", "ancom"))
    expect_identical(out$args,
                     list(lefse_1 = list(norm = "CPM", pvalue_cutoff = 0.05),
                          lefse_2 = list(norm = "TSS", pvalue_cutoff = 0.05),
                          ancom = list(norm = "CPM", W_cutoff = 0.75)))
    expect_identical(out$methods, names(out$args))


})
