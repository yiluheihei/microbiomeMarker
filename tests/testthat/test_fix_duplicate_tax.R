context("fix duplicate tax")

test_that("fix duplicate tax", {
    ps <- readRDS("data/data_tax_duplicate.rds")
    ps_fixed <- fix_duplicate_tax(ps)

    expect_fixed_length_qual <- function(ps, ps_fixed, level) {
        tax <- tax_table(ps)@.Data
        l <- length(unique(apply(tax[, 1:level], 1, paste, collapse = "_")))

        tax_fixed <- tax_table(ps_fixed)@.Data
        l_fixed <- length(unique(tax_fixed[, level]))

        expect_equal(l, l_fixed)
    }

    expect_fixed_length_qual(ps, ps_fixed, 7)
    expect_fixed_length_qual(ps, ps_fixed, 6)
    expect_fixed_length_qual(ps, ps_fixed, 5)
    expect_fixed_length_qual(ps, ps_fixed, 4)
    expect_fixed_length_qual(ps, ps_fixed, 3)
    expect_fixed_length_qual(ps, ps_fixed, 2)
})
