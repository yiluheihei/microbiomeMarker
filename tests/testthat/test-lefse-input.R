context("lefse input preprocessing")

test_that("add missing levels: keep abundance is lower than 1", {
  oxygen_feature <- otu_table(oxygen)
  feature <- add_missing_levels(oxygen_feature)

  expect_true(max(feature) <= 1)
  expect_true(min(feature) >= 0)
})

test_that("check whether taxa has level prefix", {
  prefix <- paste0(c("k", "p", "c", "o", "f", "g", "s"), "__")
  tax_nms1 <- c("Bacteria|Verrucomicrobia|Verrucomicrobiae|")
  tax_nms2 <- paste0(prefix, tax_nms1)
  check1 <- check_tax_prefix(tax_nms1)
  check2 <- purrr::map_lgl(tax_nms2, check_tax_prefix)

  expect_false(check1)
  expect_true(all(check2))
})
