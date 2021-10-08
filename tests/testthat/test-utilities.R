context("test utilities")

test_that("check upper first letter", {
    expect_equal(
        upper_firstletter(c("abc", "ABC", "Abc")),
        c("Abc", "Abc", "Abc")
    )
})

test_that("check whether taxonomic ranks include in available_ranks", {
    ps <- phyloseq::phyloseq(
        otu_table = otu_table(
            matrix(
                sample(100, 40),
                nrow = 4,
                dimnames = list(
                    paste0("otu", 1:4),
                    paste0("sample", 1:10)
                )
            ),
            taxa_are_rows = TRUE
        ),
        tax_table = tax_table(
            matrix(
                c(rep("g1", 4), rep(c("s1", "s2"), 2)),
                nrow = 4, byrow = FALSE,
                dimnames = list(paste0("otu", 1:4), c("Genus", "Species"))
            )
        )
    )
    ps2 <- phyloseq::phyloseq(
        otu_table = otu_table(
            matrix(
                sample(100, 40),
                nrow = 4,
                dimnames = list(
                    paste0("otu", 1:4),
                    paste0("sample", 1:10)
                )
            ),
            taxa_are_rows = TRUE
        ),
        tax_table = tax_table(
            matrix(
                c(rep("g1", 4), rep(c("s1", "s2"), 2)),
                nrow = 4, byrow = FALSE,
                dimnames = list(paste0("otu", 1:4), c("xxx", "Species"))
            )
        )
    )
    expect_true(check_rank_names(ps))

    expect_false(check_rank_names(ps2))
})

test_that(
    "taxa prefix",
    expect_identical(get_prefix(c("Phylum", "Genus")), c("p__", "g__"))
)

test_that("var in sample_data", {
    expect_error(
        check_var_in_meta("abc", sample_data(ecam)),
        "variable of `sample_meta`"
    )
    expect_silent(check_var_in_meta("delivery", sample_data(ecam)))
})


test_that("get_norm_method works well", {
    expect_identical(
        get_norm_method(100),
        "per-sample normalized (sum of all taxa) to 100"
    )
    expect_identical(get_norm_method("a"), "a")
})


test_that("check_samples, at least one non zero features in a sample", {
    test_ot <- otu_table(
        cbind(matrix(1:12, 6, 2), 0),
        taxa_are_rows = TRUE
    )
    test_sa <- sample_data(
        data.frame(sample = c("sa1", "sa2", "sa3"))
    )
    test_ps <- phyloseq(test_ot, test_sa)

    test_ot2 <- otu_table(
        cbind(matrix(1:12, 6, 2)),
        taxa_are_rows = TRUE
    )
    expect_identical(check_samples(test_ps), "sa3")
    expect_null(check_samples(test_ot2))
})

test_that("remove samples with missing values in the specified var", {
    otu <- otu_table(
        data.frame(
            s1 = runif(10),
            s2 = runif(10),
            s3 = runif(10)
        ),
        taxa_are_rows = TRUE
    )
    tax <- tax_table(data.frame(feature = paste0("sp", 1:10)) %>% as.matrix())
    sam <- data.frame(group = c(NA, "A", "B"))
    rownames(sam) <- paste0("s", 1:3)
    test_ps <- phyloseq(otu, tax, sample_data(sam))

    new_samples <- sample_names(remove_na_samples(test_ps, "group"))
    expect_identical(new_samples, c("s2", "s3"))
})


# create contrast
test_that("create contrast", {
    expect_error(create_contrast("a"), "at least two groups")

    # multiple groups, pairwise comparisons
    groups <- factor(rep(c("a", "b", "c"), each = 3))
    mat <- matrix(c(-1, 1, 0, -1, 0, 1, 0, -1, 1), 3)
    row.names(mat) <- c("a", "b", "c")
    colnames(mat) <- c("b-a", "c-a", "c-b")
    expect_identical(create_pairwise_contrast(levels(groups)), mat)
    expect_error(create_contrast(groups, c("a", "d")), "contained in `groups`")
    expect_error(create_contrast(groups, "a"), "two length")
    expect_identical(
        create_contrast(groups, c("a", "b")),
        matrix(c(-1, 1, 0), dimnames = list(c("a", "b", "c"), "b-a"))
    )
    expect_identical(
        create_contrast(groups, c("b", "c")),
        matrix(c(0, -1, 1), dimnames = list(c("a", "b", "c"), "c-b"))
    )

    # create contrast
    groups_two <- factor(rep(c("a", "b"), each = 3))
    expect_identical(create_contrast(groups_two), c(-1, 1))
    expect_identical(create_contrast(groups), mat)
    expect_warning(
        ctra <- create_contrast(groups_two, c("a", "b")),
        "`contrast` is ignored"
    )
    expect_identical(ctra, c(-1, 1))
    expect_identical(ctra, create_contrast(groups_two))
})


# return marker
test_that("marker_table, if no significant marker return all the features", {
    sig_ft1 <- data.frame()
    ft <- data.frame(feature = letters[1:3], ef = runif(3))
    expect_warning(return_marker(sig_ft1, ft), "No significant feature")
    expect_equal(marker_table(ft), return_marker(sig_ft1, ft))

    sig_ft2 <- data.frame(feature = "a", ef = 1)
    expect_equal(marker_table(sig_ft2), return_marker(sig_ft2, ft))
})


# extract the specific taxa rank
test_that("extract the specific taxa rank", {
    otu <- otu_table(
        data.frame(
            s1 = runif(10),
            s2 = runif(10),
            s3 = runif(10)
        ),
        taxa_are_rows = TRUE
    )
    tax <- tax_table(
        data.frame(
            rank1 = paste0("rank1", 1:10),
            rank2 = paste0("rank2", 1:10)
        ) %>%
            as.matrix()
    )
    test_ps <- phyloseq(otu, tax)

    # taxa names keep inconsistent with the taxa_rank
    expect_identical(
        taxa_names(extract_rank(test_ps, "rank2")),
        paste0("rank2", 1:10)
    )
    expect_identical(
        taxa_names(extract_rank(test_ps, "none")),
        paste0("sp", 1:10)
    )
})
