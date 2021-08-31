context("`abundances` methods")

ps <- phyloseq::phyloseq(
  otu_table = otu_table(
    matrix(
      sample(100, 40), nrow = 2,
      dimnames = list(
        c("feature1", "feature2"),
        paste0("sample", 1:20))),
    taxa_are_rows = TRUE),
  tax_table = tax_table(
    matrix(
      c("taxa1", "taxa2"), nrow = 2,
      dimnames = list(c("feature1", "feature2"), c("Species")))),
  sam_data = sample_data(
    data.frame(
      group = rep(c("group1", "group2"), 10),
      row.names = paste0("sample", 1:20)))
)

test_that("abundances",{
  abd <- abundances(ps, "log10p")
  abd_t <- transform_abundances(ps, "log10p")

  # expect_warning(abundances(sample_ps, norm = TRUE), "not been normalized")
  expect_identical(abd, as(otu_table(abd_t), "matrix"))
})

test_that("abundances normalization", {
  abd <- abundances(ps)
  abd_norm <- normalize(ps, "CSS")
  nf <- get_norm_factors(abd_norm)

  expect_identical(
    abundances(abd_norm, norm = TRUE),
    sweep(abd, 2, nf, "/")
  )
  expect_identical(abd, abundances(abd_norm))
})
