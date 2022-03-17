test_that("confounder analysis", {
    data("caporaso")
    res <- withr::with_seed(
        2022,
        confounder(caporaso, 
            "SampleType", 
            confounders = "ReportedAntibioticUsage"
        )
    )
    expect_identical(res$confounder, "ReportedAntibioticUsage")
    expect_equal(res$pvalue, 0.239)
})
