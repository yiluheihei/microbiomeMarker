context("lefse input preprocessing")

test_that("add missing levels: keep abundance is lower than 1", {
    data(oxygen)
    oxygen_feature <- otu_table(oxygen)
    feature <- add_missing_levels(oxygen_feature)

    expect_true(max(feature) <= 1)
    expect_true(min(feature) >= 0)
})

test_that("check whether taxa has level prefix", {
    prefix <- paste0(c("k", "p", "c", "o", "f", "g", "s"), "__")
    tax_nms1 <- c("Bacteria|Verrucomicrobia|Verrucomicrobiae|")
    tax_nms2 <- paste0(prefix, tax_nms1)
    check1 <- check_tax_prefix(tax_nms1)
    check2 <- purrr::map_lgl(tax_nms2, check_tax_prefix)

    expect_false(check1)
    expect_true(all(check2))
})

test_that("format lefse input: subgroups, group_hie", {
    sample_meta <- data.frame(group = rep(c("cr", "exp"), 3),
                              subgroup = rep(c("a", "b"), each = 3))
    group_info <- lefse_format_grp(sample_meta, "group", "subgroup")
    expect_identical(group_info$subgroup, 
                     paste(sample_meta$group, sample_meta$subgroup, sep = "_"))
    expect_identical(group_info$group_hie, 
                     list(cr = c("cr_a", "cr_b"), exp = c("exp_a", "exp_b")))
})
