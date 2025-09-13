
# Helper to simulate simple right-censored survival data with a Cox signal
sim_surv_data <- function(n = 120, p = 6, beta = c(1.3, -1.1, rep(0, p-2)), seed = 42) {
  stopifnot(length(beta) == p)
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", seq_len(p))
  eta <- as.vector(X %*% beta)
  # Event times: exponential with rate proportional to exp(eta)
  rate <- exp(scale(eta, center = TRUE, scale = FALSE))
  T <- rexp(n, rate = as.numeric(rate))
  # Censoring times, chosen to give a reasonable (~20-40%) censoring fraction
  C <- rexp(n, rate = 0.5)
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  list(X = X, time = time, event = event, eta = eta)
}
