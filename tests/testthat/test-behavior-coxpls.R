
test_that("coxpls recovers risk ranking and returns expected structure", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")
  dat <- sim_surv_data(n = 120, p = 6, beta = c(1.3, -1.1, 0, 0, 0, 0), seed = 123)
  # allres = TRUE for richer object
  fit_all <- coxpls(Xplan = dat$X, time = dat$time, event = dat$event,
                    type = "right", ncomp = 2, allres = TRUE, plot = FALSE, scaleX = TRUE)
  expect_type(fit_all, "list")
  expect_true(all(c("tt_pls","cox_pls","pls_mod") %in% names(fit_all)))
  expect_s3_class(fit_all$cox_pls, "coxph")
  expect_equal(ncol(fit_all$tt_pls), 2)
  eta_hat <- as.numeric(fit_all$cox_pls$linear.predictors)
  expect_length(eta_hat, nrow(dat$X))
  rho <- suppressWarnings(stats::cor(eta_hat, dat$eta, method = "spearman", use = "complete.obs"))
  expect_gte(rho, 0.4)
  # allres = FALSE should return the coxph fit directly
  fit <- coxpls(Xplan = dat$X, time = dat$time, event = dat$event,
                type = "right", ncomp = 2, allres = FALSE, plot = FALSE, scaleX = TRUE)
  expect_s3_class(fit, "coxph")
})
