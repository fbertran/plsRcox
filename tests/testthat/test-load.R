test_that("package loads quietly", {
  expect_true(requireNamespace("plsRcox", quietly = TRUE))
  expect_no_condition(suppressPackageStartupMessages(library(plsRcox)))
})
