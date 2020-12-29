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
  msg1 <- "slot `otu_table` and `tax_table` are required"
  msg2 <-  paste0(
    "The number of different feature must be smaller than the",
    " total number of feature"
  )
  msg3 <- "marker in marker_table must be contained in `taxa_names`"

  expect_microbiomeMarker_error(msg1, marker1)
  expect_microbiomeMarker_error(msg1, marker1, otu1)
  expect_microbiomeMarker_error(msg1, marker1, tax1)
  expect_microbiomeMarker_error(msg2, marker1, tax2, otu2, tax2)
  expect_microbiomeMarker_error(msg3, marker2, tax2, otu1, tax1)

  # not reaise error, since the phyloseq does not requires otu_table
  # and tax_table have the same row
  # expect_error(microbiomeMarker(diff_feature, otu1, tax1[1:8, ]))

  expect_equal(
    is(microbiomeMarker(marker1, tax2, otu1, tax1)),
    c("microbiomeMarker", "phyloseq")
  )
})

# test_that("marker are contained in the original taxa", {
#   skip_on_cran()
#   skip_on_bioc()
#
#   expect_error(
#     lefse(caporaso,
#       norm = 1e6,
#       class = "SampleType",
#       multicls_strat = TRUE
#     ),
#     NA
#   )
# })
