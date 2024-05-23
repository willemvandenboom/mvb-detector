## Helper function for generating results on simulated data
run_inference <- function (
  data, cens_level = 0, cens_time = NULL, n_iter = 1e5L, verbose = FALSE,
  cp_constraint = TRUE
) {
  n_risk <- nlevels(data$event) - 1L
  n <- nrow(data)
  p <- ncol(data) - 2L
  X <- as.matrix(data[, 2L + seq_len(p)])
  n_cens <- as.integer(cens_level * n)
  
  if (n_cens > 0L) {
    if (missing(cens_time)) stop("`cens_time` is required for censoring.")
    tmp <- which(!is.na(cens_time))
    if (length(tmp) < n_cens) stop("Too few available for censoring")
    if (tmp[n_cens] != n_cens) warning("Skipping over `cens_time = NA`")
    cens_ind <- tmp[1:n_cens]
    data$time[cens_ind] <- cens_time[cens_ind]
    data$event[cens_ind] <- "censored"
  }
  
  ## Run the multivariate Bernoulli detector.
  set.seed(1L)
  if (verbose) print("Running the multivariate Bernoulli detector...")
  t_start <- Sys.time()
  
  mvbd_fit <- mvb.detector::run_mvb_detector(
    data = data, t_max = t_max, n_iter = n_iter, show_progress_bar = verbose,
    cp_constraint = cp_constraint
  )
  
  t_end <- Sys.time()
  
  
  ## Comparison with maximum likelihood estimation
  # Transform the data to long format.
  data_long <- discSurv::dataLongCompRisks(
    dataShort = data, timeColumn = "time", eventColumn = "event",
    eventColumnsAsFactor = TRUE
  )
  
  X_long <- as.matrix(data_long[, 5L + n_risk + seq_len(p)])
  data_long$y <- 0L
  tmp <- which(data_long$e0 == 0)
  data_long$y[tmp] <- as.integer(data_long$event[tmp]) - 1L
  
  nnet_fit <- nnet::multinom(formula = if (p == 0L) {
    y ~ as.factor(timeInt)
  } else {
    y ~ as.factor(timeInt) + X_long
  }, data = data_long, maxit = 1e3L, trace = FALSE)
  
  
  ## Comparision with `brea`
  # Semiparametric model by King and Weiss (2021)
  set.seed(1L)
  if (verbose) print("Running the model by King and Weiss (2021)...")
  
  brea_fit <- brea::brea_mcmc(
    x = cbind(data_long$timeInt, X_long + 1L),
    y = as.matrix(data_long[, 3L + 1:n_risk]),
    priors = c(list(list("gmrf", 3, .01)), rep(list(list("cat", 4)), p)),
    S = n_iter, B = n_iter %/% 10L, K = c(t_max, rep(2L, p))
  )
  
  
  ## Estimate alpha and beta
  alpha_mvbd <- colMeans(mvbd_fit$alpha_MCMC)
  beta_mvbd <- colMeans(mvbd_fit$beta_MCMC)
  
  alpha_nnet <- coef(nnet_fit)[, 1:t_max]
  alpha_nnet[, 2:t_max] <- alpha_nnet[, 1] + alpha_nnet[, 2:t_max]
  beta_nnet <- t(coef(nnet_fit)[, t_max + seq_len(p)])
  
  alpha_brea <- matrix(data = NA_real_, nrow = n_risk, ncol = t_max)
  beta_brea <- matrix(data = NA_real_, nrow = p, ncol = n_risk)
  offset <- rowMeans(brea_fit$b_0_s)
  
  for (r in 1:n_risk) {
    alpha_brea[r, ] <- offset[r] + rowMeans(brea_fit$b_m_s[[1]][r, , ])
    
    for (j in seq_len(p)) {
      tmp <- brea_fit$b_m_s[[1L + j]][r, , ]
      beta_brea[j, r] <- mean(tmp[2, ] - tmp[1, ])
    }
  }
  
  return (list(
    mvbd_fit = mvbd_fit, nnet_fit = nnet_fit, brea_fit = brea_fit,
    alpha_mvbd = alpha_mvbd, beta_mvbd = beta_mvbd, alpha_nnet = alpha_nnet,
    beta_nnet = beta_nnet, alpha_brea = alpha_brea, beta_brea = beta_brea,
    t_start = t_start, t_end = t_end
  ))
}
