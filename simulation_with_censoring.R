### Simulation study with censoring

# This code produces the simulation study with censoring in Section 3 of the
# paper "The Multivariate Bernoulli detector: Change point estimation in
# discrete survival analysis" by Willem van den Boom, Maria De Iorio, Fang Qian
# and Alessandra Guglielmi.


source("mvb_detector.R")



## Simulate data


simulate_data <- function (n, alpha_true) {
  n_risk <- nrow(alpha_true)
  t_max <- ncol(alpha_true)
  time_to_event <- integer(n)
  event <- integer(n)  # Record which of the competing events happens
  
  for (i in seq_len(n)) {
    
    for (t in 1:t_max) {
      
      lam <- exp(log_hazard_rate_compete(alpha_true[, t]))
      U <- runif(1)
      
      for (r in 1:n_risk) if ({U <- U - lam[r]} < 0) {
        event[i] <- r
        break
      }
      
      if (event[i] != 0L) break
      
    }
    
    time_to_event[i] <- if (event[i] != 0L) t else t + 1L
    
  }
  
  cens <- time_to_event > t_max
  time_to_event[cens] <- t_max
  
  return (list(
    time_to_event = time_to_event, event = event, cens = cens
  ))
}


n_risk <- 3L  # Number of competing risks
t_max <- 20L

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
tmp <- simulate_data(n, alpha_true)
time_to_event_uncensored <- tmp$time_to_event
cens_uncensored <- tmp$cens
event <- tmp$event

# Censoring time *if* the corresponding individual is censored
cens_time <- integer(n)
set.seed(1L)

for (i in 1:n) cens_time[i] <- sample.int(
  n = time_to_event_uncensored[i] - !cens_uncensored[i], size = 1L
)

cens_level_vec <- c(0, .1, .5)
n_cens_level <- length(cens_level_vec)

# Regression terms
p <- 0L
X <- matrix(data = NA_real_, nrow = n, ncol = 0L)



## Run the multivariate Bernoulli detector for all three scenarios.

n_iter <- 1e5L  # Number of MCMC iterations
res_list <- list()

for (cens_ind in 1:n_cens_level) {
  cens_level <- cens_level_vec[cens_ind]
  print(paste("Working on censoring level", cens_level, "..."))
  n_cens <- as.integer(cens_level * n)
  time_to_event <- time_to_event_uncensored
  cens <- cens_uncensored
  time_to_event[1:n_cens] <- cens_time[1:n_cens]
  cens[1:n_cens] <- TRUE
  observed_times <- unique(time_to_event[!cens])
  
  # Specify the prior on change points:
  # Only allow a start of a regime at or right after an observed time, or time
  # 1. The last time point t_max cannot be the start of a regime, which is a
  # standard restriction in change point modelling, e.g.,
  # doi:10.1080/02664763.2011.559209.
  cp_candidates <- unique(c(
    1L, observed_times[observed_times < t_max],
    observed_times[observed_times < t_max - 1L] + 1L
  ))
  
  # Do not allow a start of a regime at time t right after a time t - 1 with an
  # observed event if time t + 1 also has an observed event.
  if (t_max > 2L) for (t in 2:(t_max - 1L)) {
    if (!(t %in% cp_candidates) | (t %in% observed_times)) next
    
    if (((t - 1L) %in% observed_times) & ((t + 1L) %in% observed_times)) {
      cp_candidates <- cp_candidates[-which(cp_candidates == t)]
    }
  }
  
  n_cp_candidates <- length(cp_candidates)
  
  
  log_prior_n_vec <- function (n_vec) {
    # Truncated prior that puts zero mass on change points at time points with
    # no observed events
    n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
    if (!all((cumsum(n_vec[-n_regimes]) + 1L) %in% cp_candidates)) return (-Inf)
    
    return (log_prior_n_regimes(n_regimes) - lchoose(
      n_cp_candidates - 1L, n_regimes - 1L
    ))
  }
  
  
  # Transform the data to long format.
  tmp <- discSurv::dataLong(dataShort = data.frame(
    time_to_event = time_to_event, eventColumn = ifelse(cens, 0L, 1L),
    event = event, X = X
  ), timeColumn = "time_to_event", eventColumn = "eventColumn")
  
  X_long <- as.matrix(tmp[, 6L + seq_len(p)])
  data_long <- tmp[, c("timeInt", "y", "event")]
  data_long$y <- as.integer(data_long$y)
  tmp <- which(data_long$y == 1L)
  data_long$y[tmp] <- data_long$event[tmp]
  data_long$event <- NULL
  
  # Run the MCMC.
  set.seed(1L)
  res <- run_MCMC(n_iter)
  
  res_list[[cens_ind]] <- list(
    cp_candidates = cp_candidates, res = res, cens = cens,
    time_to_event = time_to_event
  )
}



## Plot the results.

# Plot with the change points
cens_lab <- c("No censoring", paste(c("10%", "50%"), "censoring"))
xlab = expression(paste("Time ", italic(t)))
BF_max <- 50
pdf(file = "simulation_with_censoring.pdf", width = 8, height = 4.5)
par(mfrow = c(n_cens_level, n_risk + 1L), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  tmp_res <- res_list[[cens_ind]]
  cp_candidates <- tmp_res$cp_candidates
  res <- tmp_res$res
  cens <- tmp_res$cens
  time_to_event <- tmp_res$time_to_event
  n_cp_candidates <- length(cp_candidates)
  
  # Get rid of 10% of iterations as burnin.
  alpha_post <- res$alpha_MCMC[(n_iter %/% 10L + 1L):n_iter, , , drop = FALSE]
  
  # Plot change points
  gam_mean <- numeric(t_max)
  z_mean <- matrix(data = 1, nrow = n_risk, ncol = t_max)
  gam_mean[1] <- 1
  n_iter_post <- dim(alpha_post)[1]
  z_t_MCMC <- matrix(data = NA, nrow = n_risk, ncol = n_iter_post)
  
  for (t in 2:t_max) {
    for (r in 1:n_risk) for (s in 1:n_iter_post) {
      z_t_MCMC[r, s] <- alpha_post[s, r, t] != alpha_post[s, r, t - 1L]
    }
    
    z_mean[, t] <- rowMeans(z_t_MCMC)
    gam_mean[t] <- mean(colSums(z_t_MCMC) > 0)
  }
  
  # Compute the prior marginal probability of a change point.
  log_prior_n_regimes_normalizing_constant <- matrixStats::logSumExp(
    log_prior_n_regimes(1:n_cp_candidates)
  )
  
  log_prior_expec_n_regimes <- sum((1:n_cp_candidates) * exp(
    log_prior_n_regimes(1:n_cp_candidates) -
      log_prior_n_regimes_normalizing_constant)
  )
  
  prior_prob_cp <- (log_prior_expec_n_regimes - 1) / (n_cp_candidates - 1)
  
  Bayes_factor <- (gam_mean / (1 - gam_mean)) / (
    prior_prob_cp / (1 - prior_prob_cp)
  )
  
  Bayes_factor[is.infinite(Bayes_factor)] <- 2 * BF_max
  Bayes_factor[-cp_candidates] <- NA_real_
  Bayes_factor[1] <- NA_real_
  
  plot(
    x = 1:t_max, y = Bayes_factor, col = "darkgrey", type = "h", xlab = xlab,
    ylab = expression(paste("Bayes factor of ", italic(gamma[t]))), lwd = 3,
    ylim = c(0, BF_max), main = cens_lab[cens_ind]
  )
  
  points(x = 1:t_max, y = Bayes_factor, pch = 20, col = "black", cex = 0.75)
  
  abline(
    v = cumsum(n_vec_true[-n_regimes_true]) + 1L, lty = "dashed", lwd = 2,
    col = "red"
  )
  
  prior_prob_cp_single_risk <- prior_prob_cp * 2^(n_risk - 1L) / (2^n_risk - 1L)
  
  for (r in 1:n_risk) {
    Bayes_factor <- (z_mean[r, ] / (1 - z_mean[r, ])) / (
      prior_prob_cp_single_risk / (1 - prior_prob_cp_single_risk)
    )
    
    Bayes_factor[is.infinite(Bayes_factor)] <- 2 * BF_max
    Bayes_factor[-cp_candidates] <- NA_real_
    Bayes_factor[1] <- NA_real_
    
    plot(
      x = 1:t_max, y = Bayes_factor, col = "darkgrey", type = "h", xlab = xlab,
      ylab = c(
        expression(paste("Bayes factor of ", italic(z[1][t]))),
        expression(paste("Bayes factor of ", italic(z[2][t]))),
        expression(paste("Bayes factor of ", italic(z[3][t])))
      )[r], lwd = 3, ylim = c(0, BF_max),
      main = if (cens_ind == 1L) paste("Competing risk", r) else ""
    )
    
    points(x = 1:t_max, y = Bayes_factor, pch = 20, col = "black", cex = 0.75)
    
    abline(
      v = cumsum(n_vec_true[-n_regimes_true])[z_star_true[r, -1]] + 1L,
      lty = "dashed", lwd = 2, col = "red"
    )
  }
  
  # Compute the Bayes factor.
  n_regimes_MCMC <- res$n_regimes_MCMC
  n_regimes_post <- tail(n_regimes_MCMC, 0.9 * n_iter)
  post_prob_no_cp <- mean(n_regimes_post == 1)
  
  prior_prob_no_cp <- exp(
    log_prior_n_regimes(1L) - log_prior_n_regimes_normalizing_constant
  )
  
  print("Bayes factor:")
  print(post_prob_no_cp / prior_prob_no_cp)
}

dev.off()

# Plot with inference on alpha
pdf_height <- 4.5
ylim <- c(-12, 0)
pdf(file = "simulation_with_censoring_alpha.pdf", width = 7, height = pdf_height)
par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  tmp_res <- res_list[[cens_ind]]
  res <- tmp_res$res
  cens <- tmp_res$cens
  time_to_event <- tmp_res$time_to_event
  
  # Get rid of 10% of iterations as burnin.
  alpha_post <- res$alpha_MCMC[(n_iter %/% 10L + 1L):n_iter, , , drop = FALSE]
  alpha_mean <- colMeans(alpha_post)
  alpha_CI <- array(NA_real_, dim = c(2, n_risk, t_max))
  
  for (r in 1:n_risk) for (t in 1:t_max) {
    alpha_CI[, r, t] <- quantile(alpha_post[, r, t], c(0.025, 0.975))
  }
  
  
  for (r in 1:n_risk) {
    
    plot(
      alpha_mean[r, ], type = "n", lwd = 2, ylim = ylim, ylab = c(
        expression(italic(alpha[1][t])), expression(italic(alpha[2][t])),
        expression(italic(alpha[3][t]))
      )[r], xlab = xlab, xaxs = "i", main = paste0(
        if (r == 1L) paste0(cens_lab[cens_ind], ". ") else "", "Risk ", r
      )
    )
    
    polygon(
      x = c(1:t_max, t_max:1), y = c(alpha_CI[1, r, ], rev(alpha_CI[2, r, ])),
      border = NA, col = "darkgrey"
    )
    
    lines(alpha_mean[r, ], lwd = 2)
    lines(alpha_true[r, ], col = "red", lty = 2, lwd = 2)
  }
  
}

dev.off()



## Fit the models for the comparisons.

brea_mcmc2 <- function (x, y, priors = NULL, S = 1000, B = 100, n = NULL, K = NULL, 
    store_re = FALSE) 
{   # Modification of `brea::brea_mcmc` to show a progress bar
    xK <- brea:::check_set_predictors(x, K)
    x <- xK[[1]]
    K <- xK[[2]]
    M <- ncol(x)
    N <- nrow(x)
    yn <- brea:::check_set_outcome(y, n, N)
    y <- yn[[1]]
    n <- yn[[2]]
    R <- ncol(y)
    priors <- brea:::check_set_priors(priors, M, R)
    if (!(is.numeric(S) && length(S) == 1 && brea:::all_whole(S) && 
        S > 0.5)) {
        stop("the number of iterations to store S must be a positive whole number")
    }
    if (!is.integer(S)) 
        S <- as.integer(round(S))
    if (!(is.numeric(B) && length(B) == 1 && brea:::all_whole(B) && 
        B > -0.5)) {
        stop("the number of burn-in iterations B must be a nonnegative integer")
    }
    if (!is.integer(B)) 
        B <- as.integer(round(B))
    b_0 <- log(colSums(y)/(sum(n) - sum(y)))
    b_0_s <- matrix(0, R, S)
    b_m <- vector("list", M)
    b_m_s <- vector("list", M)
    b_m_a <- vector("list", M)
    for (m in seq_len(M)) {
        b_m[[m]] <- matrix(0, R, K[m])
        if (priors[[m]]$type %in% c("cat", "gmrf") || (priors[[m]]$type == 
            "re" && store_re)) {
            b_m_s[[m]] <- array(0, c(R, K[m], S))
        }
        b_m_a[[m]] <- rep(0L, K[m])
    }
    s_m <- vector("list", M)
    s_m_s <- vector("list", M)
    prec_m <- vector("list", M)
    for (m in seq_len(M)) {
        if (priors[[m]]$type == "gmrf") {
            s_m[[m]] <- rep(1, R)
            s_m_s[[m]] <- matrix(0, R, S)
            prec_m[[m]] <- 1/s_m[[m]]^2
        }
        else if (priors[[m]]$type == "re") {
            s_m[[m]] <- diag(R)
            s_m_s[[m]] <- array(0, c(R, R, S))
            prec_m[[m]] <- solve(s_m[[m]])
        }
    }
    eta <- matrix(0, R, N)
    eta <- eta + b_0
    for (m in seq_len(M)) eta <- eta + b_m[[m]][, x[, m]]
    nlpartfun <- -log1p(colSums(exp(eta)))
    subindex <- vector("list", M)
    y_sum <- vector("list", M)
    for (m in seq_len(M)) {
        subindex[[m]] <- vector("list", K[m])
        y_sum[[m]] <- matrix(0, R, K[m])
        for (k in seq_len(K[m])) {
            subindex[[m]][[k]] <- which(x[, m] == k)
            y_sum[[m]][, k] <- colSums(y[subindex[[m]][[k]], 
                , drop = FALSE])
        }
    }
    
    pb <- txtProgressBar(max = S + B, style = 3L)
    
    for (s in seq_len(S + B)) {
        for (m in seq_len(M)) {
            if (priors[[m]]$type == "cat") {
                for (k in seq_len(K[m])) {
                  prop_step <- rnorm(R, 0, 2.4/sqrt(R * (priors[[m]]$tau + 
                    y_sum[[m]][, k])))
                  b_star <- b_m[[m]][, k] + prop_step
                  cm <- rowMeans(b_m[[m]][, -k, drop = FALSE])
                  lpr <- sum(dnorm(b_star, cm, priors[[m]]$sd, 
                    TRUE) - dnorm(b_m[[m]][, k], cm, priors[[m]]$sd, 
                    TRUE))
                  subi <- subindex[[m]][[k]]
                  if (length(subi) == 0L) {
                    if (log(runif(1)) < lpr) {
                      b_m[[m]][, k] <- b_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                  else {
                    eta_star <- eta[, subi] + prop_step
                    nlpartfun_star <- -log1p(.colSums(exp(eta_star), 
                      R, length(subi)))
                    llr <- (crossprod(y_sum[[m]][, k], prop_step) + 
                      crossprod(n[subi], nlpartfun_star - nlpartfun[subi]))
                    if (log(runif(1)) < (llr + lpr)) {
                      b_m[[m]][, k] <- b_star
                      eta[, subi] <- eta_star
                      nlpartfun[subi] <- nlpartfun_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                }
                rms <- rowMeans(b_m[[m]])
                b_m[[m]] <- b_m[[m]] - rms
                b_0 <- b_0 + rms
            }
            else if (priors[[m]]$type == "gmrf") {
                for (k in seq_len(K[m])) {
                  fc_prec <- ifelse(k == 1L || k == K[m], 1, 
                    2) * prec_m[[m]] + y_sum[[m]][, k]
                  prop_step <- rnorm(R, 0, 2.4/sqrt(R * fc_prec))
                  b_star <- b_m[[m]][, k] + prop_step
                  if (k == 1L) {
                    lpr <- sum(dnorm(b_star, b_m[[m]][, 2L], 
                      s_m[[m]], TRUE) - dnorm(b_m[[m]][, k], 
                      b_m[[m]][, 2L], s_m[[m]], TRUE))
                  }
                  else if (k == K[m]) {
                    lpr <- sum(dnorm(b_star, b_m[[m]][, k - 1L], 
                      s_m[[m]], TRUE) - dnorm(b_m[[m]][, k], 
                      b_m[[m]][, k - 1L], s_m[[m]], TRUE))
                  }
                  else {
                    b_near <- (b_m[[m]][, k - 1L] + b_m[[m]][, 
                      k + 1L])/2
                    lpr <- sum(dnorm(b_star, b_near, s_m[[m]]/sqrt(2), 
                      TRUE) - dnorm(b_m[[m]][, k], b_near, s_m[[m]]/sqrt(2), 
                      TRUE))
                  }
                  subi <- subindex[[m]][[k]]
                  if (length(subi) == 0L) {
                    if (log(runif(1)) < lpr) {
                      b_m[[m]][, k] <- b_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                  else {
                    eta_star <- eta[, subi] + prop_step
                    nlpartfun_star <- -log1p(.colSums(exp(eta_star), 
                      R, length(subi)))
                    llr <- (crossprod(y_sum[[m]][, k], prop_step) + 
                      crossprod(n[subi], nlpartfun_star - nlpartfun[subi]))
                    if (log(runif(1)) < (llr + lpr)) {
                      b_m[[m]][, k] <- b_star
                      eta[, subi] <- eta_star
                      nlpartfun[subi] <- nlpartfun_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                }
                rms <- rowMeans(b_m[[m]])
                b_m[[m]] <- b_m[[m]] - rms
                b_0 <- b_0 + rms
                SS <- apply(b_m[[m]], 1, function(v) crossprod(diff(v)))
                prec_m[[m]] <- rgamma(R, priors[[m]]$a + 0.5 * 
                  (K[m] - 1), priors[[m]]$b + 0.5 * SS)
                s_m[[m]] <- sqrt(1/prec_m[[m]])
            }
            else if (priors[[m]]$type == "re") {
                for (k in seq_len(K[m])) {
                  b_old <- b_m[[m]][, k]
                  fc_prec <- prec_m[[m]] + diag(y_sum[[m]][, 
                    k])
                  prop_step <- 2.4/sqrt(R) * backsolve(chol(fc_prec), 
                    rnorm(R))
                  b_star <- b_old + prop_step
                  lpr <- -0.5 * ((b_star %*% prec_m[[m]] %*% 
                    b_star) - (b_old %*% prec_m[[m]] %*% b_old))
                  subi <- subindex[[m]][[k]]
                  if (length(subi) == 0L) {
                    if (log(runif(1)) < lpr) {
                      b_m[[m]][, k] <- b_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                  else {
                    eta_star <- eta[, subi] + prop_step
                    nlpartfun_star <- -log1p(.colSums(exp(eta_star), 
                      R, length(subi)))
                    llr <- (crossprod(y_sum[[m]][, k], prop_step) + 
                      crossprod(n[subi], nlpartfun_star - nlpartfun[subi]))
                    if (log(runif(1)) < (llr + lpr)) {
                      b_m[[m]][, k] <- b_star
                      eta[, subi] <- eta_star
                      nlpartfun[subi] <- nlpartfun_star
                      if (s > B) 
                        b_m_a[[m]][k] <- b_m_a[[m]][k] + 1L
                    }
                  }
                }
                prec_m[[m]] <- rWishart(1, priors[[m]]$nu + K[m], 
                  solve(priors[[m]]$s + tcrossprod(b_m[[m]])))[, 
                  , 1]
            }
        }
        if (s > B) {
            b_0_s[, s - B] <- b_0
            for (m in seq_len(M)) {
                if (priors[[m]]$type %in% c("cat", "gmrf") || 
                  (priors[[m]]$type == "re" && store_re)) {
                  b_m_s[[m]][, , s - B] <- b_m[[m]]
                }
                if (priors[[m]]$type == "gmrf") 
                  s_m_s[[m]][, s - B] <- s_m[[m]]
                if (priors[[m]]$type == "re") 
                  s_m_s[[m]][, , s - B] <- solve(prec_m[[m]])
            }
        }
      
      setTxtProgressBar(pb, s)
    }
    
    close(pb)
    list(b_0_s = b_0_s, b_m_s = b_m_s, s_m_s = s_m_s, b_m_a = b_m_a)
}


set.seed(1L)
res_list_brea <- list()
res_list_nnet <- list()

for (cens_ind in 1:n_cens_level) {
  print(paste("Working on cens_ind =", cens_ind))
  time_to_event <- res_list[[cens_ind]]$time_to_event
  cens <- res_list[[cens_ind]]$cens
  
  tmp <- discSurv::dataLong(dataShort = data.frame(
    time_to_event = time_to_event, eventColumn = ifelse(cens, 0L, 1L),
    event = event
  ), timeColumn = "time_to_event", eventColumn = "eventColumn")
  
  data_long <- tmp[, c("timeInt", "y", "event")]
  data_long$y <- as.integer(data_long$y)
  tmp <- which(data_long$y == 1L)
  data_long$y[tmp] <- data_long$event[tmp]
  data_long$event <- NULL
  
  # Semiparametric model by King and Weiss (2021)
  y_mat <- matrix(data = 0L, nrow = nrow(data_long), ncol = n_risk)
  for (r in 1:n_risk) y_mat[, r] <- data_long$y == r
  
  res_list_brea[[cens_ind]] <- brea_mcmc2(
    x = as.matrix(data_long$timeInt), y = y_mat,
    priors = list(list("gmrf", 3, .01)), S = n_iter, B = n_iter %/% 10L,
    K = max(data_long$timeInt)
  )
  
  # Maximum likelihood estimation
  res_list_nnet[[cens_ind]] <- nnet::multinom(
    formula = y ~ as.factor(timeInt), data = data_long, maxit = 1e3L
  )
}



## Plot the comparison results.

# Semiparametric model by King and Weiss (2021)
pdf(file = "simulation_with_censoring_brea.pdf", width = 7, height = pdf_height)
par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  time_to_event <- res_list[[cens_ind]]$time_to_event
  cens <- res_list[[cens_ind]]$cens
  brea_fit <- res_list_brea[[cens_ind]]
  alpha_mean <- matrix(data = NA_real_, nrow = n_risk, ncol = t_max)
  alpha_CI <- array(data = NA_real_, dim = c(2L, n_risk, t_max))
  offset <- rowMeans(brea_fit$b_0_s)
  
  for (r in 1:n_risk) {
    alpha_mean[r, ] <- offset[r] + rowMeans(brea_fit$b_m_s[[1]][r, , ])
  }
  
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
    lines(alpha_true[r, ], col = "red", lty = 2, lwd = 2)
    
  }

}

dev.off()


# Maximum likelihood estimation
pdf(file = "simulation_with_censoring_nnet.pdf", width = 7, height = pdf_height)
par(mfrow = c(n_cens_level, n_risk), mai = c(0.55, 0.6, 0.3, 0.2))

for (cens_ind in 1:n_cens_level) {
  time_to_event <- res_list[[cens_ind]]$time_to_event
  cens <- res_list[[cens_ind]]$cens
  nnet_fit <- res_list_nnet[[cens_ind]]
  nnet_fit_summary <- summary(nnet_fit)
  se_cov <- vcov(nnet_fit)
  
  alpha_ind <- 1:ncol(nnet_fit_summary$coefficients)
  alpha_mean <- nnet_fit_summary$coefficients[, alpha_ind]
  alpha_mean[, 2:t_max] <- alpha_mean[, 1] + alpha_mean[, 2:t_max]
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
    lines(alpha_true[r, ], col = "red", lty = 2, lwd = 2)
    
  }
  
}

dev.off()
