
test_that("coxsplsDR builds components and predicts risk sensibly", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("spls")
  dat <- sim_surv_data(n = 120, p = 8, beta = c(1.2, -1.0, 0.8, rep(0, 5)), seed = 321)
  fit_all <- coxsplsDR(Xplan = dat$X, time = dat$time, event = dat$event,
                       type = "right", ncomp = 2, allres = TRUE, plot = FALSE,
                       eta = 0.6, scaleX = TRUE)
  expect_type(fit_all, "list")
  expect_true(all(c("tt_splsDR","cox_splsDR","splsDR_mod") %in% names(fit_all)))
  expect_equal(ncol(fit_all$tt_splsDR), 2)
  expect_s3_class(fit_all$cox_splsDR, "coxph")
  eta_hat <- as.numeric(fit_all$cox_splsDR$linear.predictors)
  rho <- suppressWarnings(stats::cor(eta_hat, dat$eta, method = "spearman", use = "complete.obs"))
  expect_gte(rho, 0.3)
  # Direct coxph object when allres = FALSE
  fit <- coxsplsDR(Xplan = dat$X, time = dat$time, event = dat$event,
                   type = "right", ncomp = 2, allres = FALSE, plot = FALSE, eta = 0.6)
  expect_s3_class(fit, "coxph")
})
