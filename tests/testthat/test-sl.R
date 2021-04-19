test_that("supervised machine learning method workds properly", {
  skip_on_bioc()
  mm_lr <- run_sl(enterotypes_arumugam, "Gender", method = "LR")
  expect_identical(nmarker(mm_lr), 10L)

  mm_svm <- run_sl(
    enterotypes_arumugam,
    nfolds = 2,
    nrepeats = 2,
    "Gender",
    tune_length = 2,
    method = "SVM"
  )
  expect_identical(nmarker(mm_svm), 10L)

  mm_rf <- run_sl(
    enterotypes_arumugam,
    nfolds = 2,
    nrepeats = 2,
    "Gender",
    tune_length = 2,
    method = "RF",
    importance = "impurity"
  )
  expect_identical(nmarker(mm_rf), 10L)
})
