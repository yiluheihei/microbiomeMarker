context("`abundances` methods")

test_that("abundances",{
  abd <- abundances(sample_ps, "log10p")
  abd_t <- transform_abundances(sample_ps, "log10p")

  # expect_warning(abundances(sample_ps, norm = TRUE), "not been normalized")
  expect_identical(abd, as(otu_table(abd_t), "matrix"))
})

test_that("abundances normalization", {
  abd <- abundances(pediatric_ibd)
  abd_norm <- normalize(pediatric_ibd, "CSS")
  nf <- get_norm_factors(abd_norm)

  expect_identical(
    abundances(abd_norm, norm = TRUE),
    sweep(abd, 2, nf, "/")
  )
  expect_identical(abd, abundances(abd_norm))
})
