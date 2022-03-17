data("caporaso")
data("pediatric_ibd")
test_that("check confounding variables", { 
    expect_error(
        check_confounder(caporaso, "target_var"),
        "`target_var` must be contained in the meta data")
    
    expect_error(
        check_confounder(pediatric_ibd, target_var = "Class"),
        "No confounding var"
    )
    
    expect_error(
        check_confounder(caporaso, "SampleType", c("Year", "test")),
        "`test` not be contained"
    )
    
    expect_identical(
        check_confounder(caporaso, "SampleType", c("Year", "Day")),
        c("Year", "Day")
    )
    
    vars <- names(sample_data(caporaso))
    confounders <- setdiff(vars, "SampleType")
    expect_identical(
        check_confounder(caporaso, "SampleType"),
        confounders
    )
    
    
})

test_that("confounder analysis", {
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
