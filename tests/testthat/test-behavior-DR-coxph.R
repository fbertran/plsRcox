
test_that("DR_coxph returns finite residuals with right-censoring", {
  skip_if_not_installed("survival")
  dat <- sim_surv_data(n = 60, p = 4, beta = c(1.1, 0, 0, 0), seed = 99)
  res <- DR_coxph(time = dat$time, event = dat$event, type = "right", typeres = "deviance")
  expect_type(res, "double")
  expect_length(res, length(dat$time))
  expect_true(all(is.finite(res)))
})
