context("microbiomeMarker class")

test_that("microbiomeMarker constructor", {
    marker1 <- marker_table(
        data.frame(
            feature = paste0("sp", 1:5),
            enrich_group = c("cr", "er", "cr", "cr", "er"),
            stringsAsFactors = FALSE
        )
    )
    marker2 <- marker_table(
        data.frame(
            feature = paste0("sp", c(1:5, 11)),
            enrich_group = c("cr", "er", "cr", "cr", "er", "cr"),
            stringsAsFactors = FALSE
        )
    )

    otu1 <- otu_table(
        data.frame(
            s1 = runif(10),
            s2 = runif(10)
        ),
        taxa_are_rows = TRUE
    )
    tax1 <- tax_table(data.frame(feature = paste0("sp", 1:10)) %>% as.matrix())
    otu2 <- otu1[1:3, ]
    tax2 <- tax1[1:3, ]

    # expect error message for microbiomeMarker constructor
    expect_microbiomeMarker_error <- function(message, ...) {
        expect_error(microbiomeMarker(...), message, fixed = TRUE)
    }

    # otu_table is required
    expect_microbiomeMarker_error("otu_table is required", marker1)
    expect_microbiomeMarker_error(
        "otu_table is required",
        marker1,
        tax_table = tax1
    )
    # tax_table is required
    expect_microbiomeMarker_error(
        "tax_table is required",
        marker1,
        otu_table = otu1
    )

    expect_silent(microbiomeMarker(marker1, otu_table = otu1, tax_table = tax1))

    msg1 <- paste0(
        "The number of different feature must be smaller than the",
        " total number of feature"
    )
    expect_microbiomeMarker_error(
        msg1,
        marker1,
        tax_table = tax2,
        otu_table = otu2
    )

    msg2 <- "marker in marker_table must be contained in tax"
    expect_microbiomeMarker_error(
        msg2,
        marker2,
        otu_table = otu1,
        tax_table = tax1
    )

    msg4 <- "nrow of `otu_table` must be equal to the length of `tax_table()`"
    expect_microbiomeMarker_error(
        msg4,
        marker1,
        tax_table = tax1,
        otu_table = otu2,
    )

    expect_equal(
        is(microbiomeMarker(marker1, otu_table = otu1, tax_table = tax1)),
        c("microbiomeMarker", "phyloseq")
    )
})
