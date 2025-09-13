test_that("datasets load", {
  expect_no_error(data("Xmicro.censure_compl_imp", package = "plsRcox", envir = environment()))
  expect_true(exists("Xmicro.censure_compl_imp", inherits = FALSE))
  expect_no_error(data("micro.censure", package = "plsRcox", envir = environment()))
  expect_true(exists("micro.censure", inherits = FALSE))
})
