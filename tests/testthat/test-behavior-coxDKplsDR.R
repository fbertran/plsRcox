
test_that("coxDKplsDR runs (if deps installed) and returns sensible shapes", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("kernlab")
  dat <- sim_surv_data(n = 100, p = 6, beta = c(1.0, -1.1, rep(0, 4)), seed = 456)
  fit_all <- coxDKplsDR(Xplan = dat$X, time = dat$time, event = dat$event,
                        type = "right", ncomp = 2, allres = TRUE, plot = FALSE, scaleX = TRUE,
                        verbose = FALSE)
  expect_type(fit_all, "list")
  expect_true(all(c("tt_DKplsDR","cox_DKplsDR") %in% names(fit_all)))
  expect_equal(ncol(fit_all$tt_DKplsDR), 2)
  expect_s3_class(fit_all$cox_DKplsDR, "coxph")
  eta_hat <- as.numeric(fit_all$cox_DKplsDR$linear.predictors)
  rho <- suppressWarnings(stats::cor(eta_hat, dat$eta, method = "spearman", use = "complete.obs"))
  expect_gte(rho, 0.3)
})
