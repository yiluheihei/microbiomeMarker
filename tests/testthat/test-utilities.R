context("test utilities")

test_that("check upper first letter", {
  expect_equal(
    upper_firstletter(c("abc", "ABC", "Abc")),
    c("Abc", "Abc", "Abc")
  )
})

test_that(" check whether all names of taxonomic ranks include in available_ranks", {
  expect_true(check_rank_names(oxygen))
  expect_true(check_rank_names(caporaso))
})

test_that("taxa prefix",
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
