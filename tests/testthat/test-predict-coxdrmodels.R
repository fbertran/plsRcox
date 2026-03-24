test_that("coxplsDR allres objects predict from original covariates", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")

  dat <- sim_surv_data(n = 50, p = 5, beta = c(1.1, -0.9, 0, 0, 0), seed = 31)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- coxplsDR(dat$X, time = dat$time, event = dat$event, ncomp = 2, allres = TRUE, plot = FALSE)

  direct_lp <- predict(fit$cox_plsDR, newdata = fit$tt_plsDR[1:5, , drop = FALSE], type = "lp")
  wrapped_lp <- predict(fit, newdata = dat$X[1:5, , drop = FALSE], type = "lp")

  expect_equal(as.numeric(wrapped_lp), as.numeric(direct_lp), tolerance = 1e-8)

  expected_full <- predict(fit, type = "expected")
  new_x <- dat$X[1:5, , drop = FALSE]
  rownames(new_x) <- paste0("new", seq_len(nrow(new_x)))
  expected_subset <- predict(fit, newdata = new_x, type = "expected")

  expect_equal(as.numeric(expected_subset), as.numeric(expected_full[1:5]), tolerance = 1e-8)
})

test_that("coxsplsDR allres objects predict from original covariates", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("spls")

  dat <- sim_surv_data(n = 50, p = 5, beta = c(1.1, -0.9, 0, 0, 0), seed = 32)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- coxsplsDR(dat$X, time = dat$time, event = dat$event, ncomp = 2, allres = TRUE, plot = FALSE, eta = 0.6)

  direct_lp <- predict(fit$cox_splsDR, newdata = fit$tt_splsDR[1:5, , drop = FALSE], type = "lp")
  wrapped_lp <- predict(fit, newdata = dat$X[1:5, , drop = FALSE], type = "lp")

  expect_equal(as.numeric(wrapped_lp), as.numeric(direct_lp), tolerance = 1e-8)
})

test_that("coxDKsplsDR allres objects predict from original covariates", {
  skip_if_not_installed("survival")
  skip_if_not_installed("mixOmics")
  skip_if_not_installed("spls")

  dat <- sim_surv_data(n = 50, p = 5, beta = c(1.1, -0.9, 0, 0, 0), seed = 33)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- coxDKsplsDR(dat$X, time = dat$time, event = dat$event, ncomp = 2, allres = TRUE, plot = FALSE, eta = 0.6, verbose = FALSE)

  direct_lp <- predict(fit$cox_DKsplsDR, newdata = fit$tt_DKsplsDR[1:5, , drop = FALSE], type = "lp")
  wrapped_lp <- predict(fit, newdata = dat$X[1:5, , drop = FALSE], type = "lp")

  expect_equal(as.numeric(wrapped_lp), as.numeric(direct_lp), tolerance = 1e-8)
})
