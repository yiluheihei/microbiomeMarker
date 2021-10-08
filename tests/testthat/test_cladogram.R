context("test cladogram")

test_that("Generate unique id for short annotation label", {
    uid <- get_unique_id(500)
    expect_equal(uid[26], "z")
    expect_equal(uid[30], "ad")
    expect_equal(uid[260], "iz")
    expect_equal(uid[500], "sf")
})

test_that("drop the levels (no marker) if the enrich_group is a factor", {
    marker <- readRDS("data/generate_cladogram_annotation.rds")
    group_n <- length(unique(marker$enrich_group))
    expect_error(
        generate_cladogram_annotation(marker, color = colors()[1:group_n]),
        NA
    )
})

test_that("generate tree data from phyloseq", {
    data(pediatric_ibd)
    td <- get_treedata_phyloseq(pediatric_ibd)
    # all node classes must be in r, k,p, c,o, f, g, s
    expect_true(
        all(levels(td@data$node_class) %in%
            c("r", "k", "p", "c", "o", "f", "g", "s"))
    )
})
