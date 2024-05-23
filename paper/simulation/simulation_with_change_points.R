### Simulation study with change points

# This code produces the simulation study with censoring in Supplementary
# Materials of the paper "The Multivariate Bernoulli detector: Change point
# estimation in discrete survival analysis" by Willem van den Boom, Maria De
# Iorio, Fang Qian and Alessandra Guglielmi.


## Simulate data

n_risk <- 3L  # Number of competing risks

# Set change points.
n_vec_true <- c(5L, 7L, 8L)
t_max <- sum(n_vec_true)  # Number of time points

# No. of regimes is no. of change points + 1.
n_regimes_true <- length(n_vec_true)

alpha_star_true <- matrix(NA_real_, nrow = n_risk, ncol = n_regimes_true)
alpha_star_true[1, ] <- c(-9, -4, -2)
alpha_star_true[2, ] <- c(-9, -3, -2)
alpha_star_true[3, ] <- c(-9, -3, -3)

z_star_true <- matrix(data = TRUE, nrow = n_risk, ncol = n_regimes_true)

if (n_regimes_true > 1L) for (k in 2:n_regimes_true) {
  z_star_true[, k] <- alpha_star_true[, k] != alpha_star_true[, k - 1L]
}

alpha_true <- matrix(NA_real_, nrow = n_risk, ncol = t_max)
for (r in 1:n_risk) alpha_true[r, ] <- rep(alpha_star_true[r, ], n_vec_true)

# Generate survival times
n <- 300L  # Number of observations
set.seed(1L)
data <- mvb.detector::simulate_data(n, alpha_true)

# Censoring time *if* the corresponding individual is censored
cens_time <- integer(n)
set.seed(1L)

for (i in 1:n) {
  cens_time[i] <- if (with(data, time[i] == 1 && event[i] != "censored")) {
    NA_integer_
  } else sample.int(
    n = with(data, time[i] - (event[i] != "censored")), size = 1L
  )
}

cens_level_vec <- c(0, .1, .5)
n_cens_level <- length(cens_level_vec)



## Run the multivariate Bernoulli detector for all three scenarios.

source("simulation_util.R")
res_list <- vector(mode = "list", length = n_cens_level)

for (cens_ind in 1:n_cens_level) {
  cens_level <- cens_level_vec[cens_ind]
  print(paste("Working on censoring level", cens_level, "..."))
  
  res_list[[cens_ind]] <- run_inference(
    data = data, cens_level = cens_level, cens_time = cens_time, n_iter = 1e5L,
    verbose = TRUE
  )
}



## Plot the results.

# Plot the change points.
cens_lab <- c("No censoring", paste(c("10%", "50%"), "censoring"))

for (BF in c(FALSE, TRUE)) {
  pdf(
    file = paste0("simulation_with_censoring", if (BF) "_BF", ".pdf"),
    width = 8, height = 4.5
  )
  
  par(mfrow = c(n_cens_level, n_risk + 1L), mai = c(0.55, 0.6, 0.3, 0.2))
  
  for (cens_ind in 1:n_cens_level) {
    tmp_res <- res_list[[cens_ind]]
    main_vec <- c(cens_lab[cens_ind], rep("", n_risk))
    if (cens_ind == 1L) main_vec[1L + 1:n_risk] <- paste("Risk", 1:n_risk)
    
    mvb.detector::plot_change_points(
      mvbd_fit = tmp_res$mvbd_fit, n_vec_true = n_vec_true,
      z_star_true = z_star_true, BF = BF, main_vec = main_vec, BF_max = 50
    )
  }
  
  dev.off()
}


for (cens_ind in 1:n_cens_level) {
  # Compute the Bayes factor.
  tmp_res <- res_list[[cens_ind]]
  
  cp_candidates <- mvb.detector:::get_cp_candidates(
    tmp_res$mvbd_fit$data, t_max, TRUE
  )
  
  n_cp_candidates <- length(cp_candidates)
  
  prior_prob_no_cp <- exp(
    mvb.detector:::log_prior_n_regimes(1L) - mvb.detector:::logSumExp(
      mvb.detector:::log_prior_n_regimes(1:n_cp_candidates)
    )
  )
  
  print("Posterior probability of no change points:")
  print(mean(tmp_res$mvbd_fit$n_regimes_MCMC == 1))
  
  print("Bayes factor:")
  print(mean(tmp_res$mvbd_fit$n_regimes_MCMC == 1) / prior_prob_no_cp)
}


# Plot inference on alpha.
pdf_height <- 4.5
ylim <- c(-12, 0)

pdf(
  file = "simulation_with_censoring_alpha.pdf", width = 7, height = pdf_height
)

par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  main_vec <- character(n_risk)
  
  for (r in 1:n_risk) main_vec[r] <- paste0(
    if (r == 1L) paste0(cens_lab[cens_ind], ". ") else "", "Risk ", r
  )
  
  mvb.detector::plot_baseline_hazards(
    mvbd_fit = res_list[[cens_ind]]$mvbd_fit, alpha_true = alpha_true,
    ylim = ylim, main_vec = main_vec
  )
}

dev.off()


## Plot the comparison results.
xlab <- expression(paste("Time ", italic(t)))

# Semiparametric model by King and Weiss (2021)
pdf(file = "simulation_with_censoring_brea.pdf", width = 7, height = pdf_height)
par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  brea_fit <- res_list[[cens_ind]]$brea_fit
  alpha_mean <- res_list[[cens_ind]]$alpha_brea
  alpha_CI <- array(data = NA_real_, dim = c(2L, n_risk, t_max))
  offset <- rowMeans(brea_fit$b_0_s)
  
  for (r in 1:n_risk) for (i in 1:2) alpha_CI[i, r, ] <- offset[r] + apply(
    X = brea_fit$b_m_s[[1]][r, , ], MARGIN = 1L, FUN = quantile,
    probs = c(.025, .975)[i]
  )
  
  
  for (r in 1:n_risk) {
    
    plot(
    alpha_mean[r, ], type = "n", lwd = 2, ylim = ylim, ylab = c(
      expression(italic(alpha[1][t])), expression(italic(alpha[2][t])),
      expression(italic(alpha[3][t]))
    )[r], xlab = xlab, xaxs = "i", main = paste0(
      if (r == 1L) paste0(cens_lab[cens_ind], ". ") else "", "Risk ", r
    ))
    
    polygon(
      x = c(1:t_max, t_max:1), y = c(alpha_CI[1, r, ], rev(alpha_CI[2, r, ])),
      border = NA, col = "darkgrey"
    )
    
    lines(alpha_mean[r, ], lwd = 2)
    lines(alpha_true[r, ], col = "#d62728", lty = 2, lwd = 2)
    
  }

}

dev.off()


# Maximum likelihood estimation
pdf(file = "simulation_with_censoring_nnet.pdf", width = 7, height = pdf_height)
par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  nnet_fit <- res_list[[cens_ind]]$nnet_fit
  nnet_fit_summary <- summary(nnet_fit)
  se_cov <- nnet:::vcov.multinom(nnet_fit)
  
  alpha_ind <- 1:ncol(nnet_fit_summary$coefficients)
  alpha_mean <- res_list[[cens_ind]]$alpha_nnet
  alpha_se <- matrix(data = NA_real_, nrow = n_risk, ncol = t_max)
  alpha_se[1, ] <- diag(se_cov)[alpha_ind]
  alpha_se[2, ] <- diag(se_cov)[alpha_ind[t_max] + alpha_ind]
  
  alpha_se[1, -1] <-
    alpha_se[1, -1] + alpha_se[1, 1] + 2 * se_cov[1, alpha_ind[-1]]
  
  alpha_se[2, -1] <-
    alpha_se[2, -1] + alpha_se[2, 1] +
    2 * se_cov[alpha_ind[t_max] + 1L, alpha_ind[t_max] + alpha_ind[-1]]
  
  alpha_se <- sqrt(alpha_se)
  alpha_CI <- array(data = NA_real_, dim = c(2L, n_risk, t_max))
  
  for (r in 1:n_risk) for (t in 1:t_max) {
    alpha_CI[, r, t] <- alpha_mean[r, t] + c(-1, 1) * 1.96 * alpha_se[r, t]
  }
  
  for (r in 1:n_risk) {
    plot(alpha_mean[r, ], type = "n", lwd = 2, ylim = ylim, ylab = c(
      expression(italic(alpha[1][t])), expression(italic(alpha[2][t])),
      expression(italic(alpha[3][t]))
    )[r], xlab = xlab, xaxs = "i", main = paste0(
      if (r == 1L) paste0(cens_lab[cens_ind], ". ") else "", "Risk ", r
    ))
    
    polygon(
      x = c(1:t_max, t_max:1), y = c(alpha_CI[1, r, ], rev(alpha_CI[2, r, ])),
      border = NA, col = "darkgrey"
    )
    
    lines(alpha_mean[r, ], lwd = 2)
    lines(alpha_true[r, ], col = "#d62728", lty = 2, lwd = 2)
  }
}

dev.off()
