test_that("predict.plsRcoxmodel supports expected predictions on training subsets", {
  skip_if_not_installed("survival")

  dat <- sim_surv_data(n = 50, p = 4, beta = c(1.1, -0.9, 0, 0), seed = 11)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- plsRcox(dat$X, time = dat$time, event = dat$event, nt = 2, verbose = FALSE)

  expected_full <- predict(fit, type = "expected")
  expected_subset <- predict(fit, newdata = dat$X[1:5, , drop = FALSE], type = "expected")

  expect_equal(as.numeric(expected_subset), as.numeric(expected_full[1:5]), tolerance = 1e-8)
})

test_that("predict.plsRcoxmodel matches unnamed training subsets for expected predictions", {
  skip_if_not_installed("survival")

  dat <- sim_surv_data(n = 50, p = 4, beta = c(1.1, -0.9, 0, 0), seed = 12)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- plsRcox(dat$X, time = dat$time, event = dat$event, nt = 2, verbose = FALSE)
  expected_full <- predict(fit, type = "expected")

  unnamed_subset <- dat$X[1:5, , drop = FALSE]
  rownames(unnamed_subset) <- paste0("new", seq_len(nrow(unnamed_subset)))

  expected_subset <- predict(fit, newdata = unnamed_subset, type = "expected")

  expect_equal(as.numeric(expected_subset), as.numeric(expected_full[1:5]), tolerance = 1e-8)
})

test_that("predict.plsRcoxmodel can use an explicit Surv response for new expected predictions", {
  skip_if_not_installed("survival")

  dat <- sim_surv_data(n = 50, p = 4, beta = c(1.1, -0.9, 0, 0), seed = 12)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- plsRcox(dat$X, time = dat$time, event = dat$event, nt = 2, verbose = FALSE)

  new_x <- dat$X[1:5, , drop = FALSE]
  rownames(new_x) <- paste0("new", seq_len(nrow(new_x)))
  new_x[, 1] <- new_x[, 1] + 0.1

  expect_error(predict(fit, newdata = new_x, type = "expected"), "requires follow-up information")

  new_y <- survival::Surv(dat$time[1:5], dat$event[1:5])
  expected_new <- predict(fit, newdata = new_x, type = "expected", y = new_y)
  scaled_y <- new_y
  scaled_y[, 1] <- (scaled_y[, 1] - attr(fit$RepY, "scaled:center")) / attr(fit$RepY, "scaled:scale")
  score_frame <- data.frame(
    YwotNA = scaled_y,
    tt = unname(as.matrix(predict(fit, newdata = new_x, type = "scores")))
  )
  direct_expected <- predict(fit$FinalModel, newdata = score_frame, type = "expected")

  expect_equal(as.numeric(expected_new), as.numeric(direct_expected), tolerance = 1e-8)
})

test_that("predict.plsRcoxmodel forwards reference to predict.coxph", {
  skip_if_not_installed("survival")

  dat <- sim_surv_data(n = 50, p = 4, beta = c(1.1, -0.9, 0, 0), seed = 13)
  rownames(dat$X) <- seq_len(nrow(dat$X))

  fit <- plsRcox(dat$X, time = dat$time, event = dat$event, nt = 2, verbose = FALSE)

  new_x <- dat$X[1:5, , drop = FALSE]
  score_frame <- data.frame(tt = unname(as.matrix(predict(fit, newdata = new_x, type = "scores"))))

  direct_lp <- predict(fit$FinalModel, newdata = score_frame, type = "lp", reference = "zero")
  wrapped_lp <- predict(fit, newdata = new_x, type = "lp", reference = "zero")

  expect_equal(as.numeric(wrapped_lp), as.numeric(direct_lp), tolerance = 1e-8)
})
