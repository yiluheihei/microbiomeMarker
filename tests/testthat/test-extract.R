test_that("extract methods", {
    marker <- marker_table(
        data.frame(
            feature = paste0("sp", 1:5),
            enrich_group = c("cr", "er", "cr", "cr", "er"),
            stringsAsFactors = FALSE
        )
    )

    expect_s4_class(marker[1], "marker_table")
    expect_s4_class(marker[1, ], "marker_table")
})
