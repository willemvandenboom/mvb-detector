### Multivariate Bernoulli detector

# This code implements the Markov chain Monte Carlo described in the paper
# "The Multivariate Bernoulli detector: Change point estimation in discrete
# survival analysis" by Willem van den Boom, Maria De Iorio, Fang Qian and
# Alessandra Guglielmi.


## Move the following to the correct spot outside the R package.
# # Install the required packages.
# for (tmp in c("brea", "matrixStats", "nnet")) {
#   if(!tmp %in% rownames(installed.packages())) install.packages(
#     pkgs = tmp, repos = "https://cloud.r-project.org", dependencies = TRUE
#   )
# }



## Code for the augmented updates, i.e. local MCMC


# Code that provides indices that assign the discrete survival data in long
# format to their respective regimes which are marked by the change points:
get_ind_list <- function (n_vec, z_star, data_long) {
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  regime <- .bincode(x = data_long$timeInt, breaks = c(0L, cumsum(n_vec)))
  # ind_list_tmp <- list()
  # for (k in 1:n_regimes) ind_list_tmp[[k]] <- which(regime == k)
  ind_list_tmp <- get_ind_list_Rcpp(regime, n_regimes)
  
  ind_list <- list()
  n_risk <- nrow(z_star)
  
  for (r in 1:n_risk) {
    ind_list[[r]] <- list(ind_list_tmp[[1]])
    ind <- 1L
    
    if (n_regimes > 1L) for (k in 2:n_regimes) if (z_star[r, k]) {
      ind <- ind + 1L
      ind_list[[r]][[ind]] <- ind_list_tmp[[k]]
    } else ind_list[[r]][[ind]] <- c(ind_list[[r]][[ind]], ind_list_tmp[[k]])
  }
  
  ind_list[[n_risk + 1L]] <- ind_list_tmp
  return (ind_list)
}


sample_augmentation <- function (
  alpha_star, beta, ind_list, data_long, X_long
) {
  n_long <- nrow(data_long)
  n_risk <- nrow(alpha_star)
  n_regimes <- ncol(alpha_star)  # No. of regimes is no. of change points + 1.
  z <- matrix(data = NA_real_, nrow = n_risk, ncol = n_long)
  m <- z
  phi <- z
  
  for (k in 1:n_regimes) {
    ind <- ind_list[[n_risk + 1L]][[k]]
    
    # Sample the data augmentation.
    # This follows step (b) in Section 3.2 of doi:10.1016/j.csda.2006.10.006.
    n_long_k <- length(ind)
    eta <- matrix(NA_real_, nrow = n_risk, ncol = n_long_k)
    
    for (r in 1:n_risk) eta[r, ] <- alpha_star[r, k] + X_long[
      ind, , drop = FALSE
    ] %*% beta[, r, drop = FALSE]
    
    lam <- exp(eta)
    U_i <- stats::runif(n_long_k)
    colSums_lam <- colSums(lam)
    
    for (r in 1:n_risk) z[r, ind] <- -log(-log(U_i) / (1 + colSums_lam) - log(
      stats::runif(n_long_k)
    ) / lam[r, ] * (data_long$y[ind] != r))
    
    # Table 1 of doi:10.1016/j.csda.2006.10.006
    w_mix <- c(.00397, .0396, .168, .147, .125, .101, .104, .116, .107, .088)
    m_mix <- c(5.09, 3.29, 1.82, 1.24, .764, .391, .0431, -.306, -.673, -1.06)
    phi_mix <- c(4.50, 2.02, 1.10, .422, .198, .107, .0778, .0766, .0947, .146)
    
    tmp1 <- log(w_mix) - .5 * log(phi_mix)
    tmp2 <- .5 / phi_mix
    
    U <- matrix(
      data = stats::runif(n_risk * n_long_k), nrow = n_risk, ncol = n_long_k
    )
    
    r_i_mat <- sample_r_i_mat(
      U, z[, ind, drop = FALSE] - eta, m_mix, tmp1, tmp2
    )
    
    for (r in 1:n_risk) {
      m[r, ind] <- m_mix[r_i_mat[r, ]]
      phi[r, ind] <- phi_mix[r_i_mat[r, ]]
    }
  }
  
  return (list(z = z - m, phi = phi))
}


update_alpha <- function (z, phi, z_star, ind_list) {
  n_risk <- nrow(z)
  n_regimes <- ncol(z_star)  # No. of regimes is no. of change points + 1.
  alpha_star <- matrix(data = NA_real_, nrow = n_risk, ncol = n_regimes)
  
  for (r in 1:n_risk) for (k in 1:n_regimes) {
    if (!z_star[r, k]) {
      alpha_star[r, k] <- alpha_star[r, k - 1L]
      next
    }
    
    ind <- ind_list[[r]][[sum(z_star[r, 1:k])]]
    
    # Sample alpha.
    prec_post_alpha <- 1 / s2_alpha + sum(1 / phi[r, ind])
    
    mean_post_alpha <- (
      mu_alpha / s2_alpha + sum(z[r, ind] / phi[r, ind])
    ) / prec_post_alpha
    
    alpha_star[r, k] <- mean_post_alpha + stats::rnorm(1L) / sqrt(
      prec_post_alpha
    )
  }
  
  return (alpha_star)
}


compute_log_marginal_likelihood <- function (ind_list, z, phi) {
  # Marginal likelihood with the data augmentation
  n_risk <- nrow(z)
  n_long <- ncol(z)
  phi_inv <- 1 / phi
  z_div_phi <- z / phi
  log_lik <- 0
  
  for (r in 1:n_risk) {
    n_regimes <- length(ind_list[[r]])
    
    # Create "assignment" matrix that links each observation with their regime.
    # tau_mat <- matrix(data = FALSE, nrow = n_long, ncol = n_regimes)
    # for (k in 1:n_regimes) tau_mat[ind_list[[r]][[k]], k] <- TRUE
    
    # Likelihood follows from
    # z ~ N(mu_alpha, s2_alpha * tau_mat %*% t(tau_mat) + diag(phi))
    
    # for (r in 1:n_risk) log_lik <- log_lik + mvtnorm::dmvnorm(
    #   x = z[r, ] - mu_alpha,
    #   sigma = s2_alpha * tcrossprod(tau_mat) + diag(phi[r, ]),
    #   log = TRUE
    # )
    
    # Gaussian density evaluation using the matrix determinant lemma and the
    # Woodbury matrix identity, and dropping constant terms (those that do not
    # vary with `ind_list`):
    
    # tmp_mat <-
    #   diag(n_regimes) + s2_alpha * crossprod(tau_mat / phi[r, ], tau_mat)
    # 
    # tmp <- crossprod(tau_mat, (z[r, ] - mu_alpha) / phi[r, ])
    
    tmp_vec <- numeric(n_regimes)  # tmp_mat = diag(tmp_vec)
    tmp <- numeric(n_regimes)
    
    for (k in 1:n_regimes) {
      ind <- ind_list[[r]][[k]]
      
      # tmp_vec[k] <- 1 + s2_alpha * sum(1 / phi[r, ind])
      # tmp[k] <- sum((z[r, ind] - mu_alpha) / phi[r, ind])
      
      sum_phi_inv <- sum(phi_inv[r, ind])
      tmp_vec[k] <- 1 + s2_alpha * sum_phi_inv
      tmp[k] <- sum(z_div_phi[r, ind]) - mu_alpha * sum_phi_inv
    }
    
    # log_lik <- log_lik - 0.5 * determinant(tmp_mat)$modulus +
    #   0.5 * s2_alpha * crossprod(tmp, solve(tmp_mat, tmp))
    
    log_lik <-
      log_lik - 0.5 * sum(log(tmp_vec)) + 0.5 * s2_alpha * sum(tmp^2 / tmp_vec)
  }
  
  return (log_lik)
}


log_posterior_augmented <- function (
  n_vec, z_star, z, phi, data_long, log_prior_n_vec
) {
  # Here, the prior on cause-specific change points p(z_t | gamma_t) is not
  # included in `log_prior_n_vec`, and thus not included in this log-posterior
  # computation.
  ind_list <- get_ind_list(n_vec, z_star, data_long)
  
  return (
    log_prior_n_vec(n_vec) + compute_log_marginal_likelihood(ind_list, z, phi)
  )
}


update_n_vec_local <- function (
  n_vec, z_star, z, phi, data_long, n_risk, t_max, log_prior_n_vec
) {
  q <- 0.5  # Split proposal probability
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  
  log_post <- log_posterior_augmented(
    n_vec, z_star, z, phi, data_long, log_prior_n_vec
  )
  
  if (n_regimes == 1 || (n_regimes < t_max && stats::runif(1) < q)) {  # Split
    splittable <- (1:n_regimes)[n_vec > 1L]
    
    if (length(splittable) > 1L) {
      j <- sample(x = (1:n_regimes)[n_vec > 1L], size = 1L)
    } else {
      j <- splittable[1]
    }
    
    l <- sample.int(n = n_vec[j] - 1L, size = 1L)
    
    n_vec_prop <- c(
      utils::head(n_vec, j - 1L), l, n_vec[j] - l,
      utils::tail(n_vec, n_regimes - j)
    )
    
    z_new <- sample_z(n_risk)
    
    z_star_prop <- cbind(
      z_star[, seq_len(j), drop = FALSE],
      z_new,
      z_star[, j + seq_len(n_regimes - j), drop = FALSE]
    )
    
    log_post_prop <- log_posterior_augmented(
      n_vec_prop, z_star_prop, z, phi, data_long, log_prior_n_vec
    )
    
    if (stats::runif(1) < exp(
      (n_regimes < t_max - 1L) * log(1 - q) - (n_regimes > 1) * log(q) +
        log(length(splittable)) + log(n_vec[j] - 1L) - log(n_regimes) +
        log_post_prop - log_post
    )) {
      n_vec <- n_vec_prop
      z_star <- z_star_prop
      log_post <- log_post_prop
    }
  } else {  # Merge
    j <- sample.int(n = n_regimes - 1L, size = 1L)
    n_merge <- n_vec[j] + n_vec[j + 1L]
    
    n_vec_prop <- c(
      utils::head(n_vec, j - 1L), n_merge,
      utils::tail(n_vec, n_regimes - j - 1L)
    )
    
    z_star_prop <- cbind(
      z_star[, seq_len(j), drop = FALSE],
      z_star[, j + 1L + seq_len(n_regimes - j - 1L), drop = FALSE]
    )
    
    log_post_prop <- log_posterior_augmented(
      n_vec_prop, z_star_prop, z, phi, data_long, log_prior_n_vec
    )
    
    if (stats::runif(1) < exp(
      (n_regimes > 2L) * log(q) - (n_regimes < t_max) * log(1 - q) +
        log(n_regimes - 1L) - log(sum(n_vec_prop > 1L)) - log(n_merge - 1L) +
        log_post_prop - log_post
    )) {
      n_vec <- n_vec_prop
      z_star <- z_star_prop
      log_post <- log_post_prop
    }
  }
  
  n_regimes <- length(n_vec)
  
  if (n_regimes > 1) {  # Shuffle while keeping `z_star` fixed
    i <- sample.int(n = n_regimes - 1L, size = 1L)
    j <- sample.int(n = n_vec[i] + n_vec[i + 1L] - 1L, size = 1L)
    n_vec_prop <- n_vec
    n_vec_prop[i + 1L] <- n_vec[i] + n_vec[i + 1L] - j
    n_vec_prop[i] <- j
    
    log_post_prop <- log_posterior_augmented(
      n_vec_prop, z_star, z, phi, data_long, log_prior_n_vec
    )
    
    if (stats::runif(1) < exp(log_post_prop - log_post)) {
      n_vec <- n_vec_prop
      log_post <- log_post_prop
    }
  }
  
  return (list(log_post = log_post, n_vec = n_vec, z_star = z_star))
}


update_z_local <- function (
  log_post, n_vec, z_star, z, phi, data_long, n_risk, log_prior_n_vec
) {
  # Resample a randomly selected z_star.
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.

  if (n_regimes > 1L) {
    j <- if (n_regimes == 2L) 2L else sample(x = 2:n_regimes, size = 1L)
    z_new <- sample_z(n_risk)
    
    if (all(z_new == z_star[, j, drop = FALSE])) return (
      list(log_post = log_post, z_star = z_star)
    )

    z_star_prop <- z_star
    z_star_prop[, j] <- z_new
    
    log_post_prop <- log_posterior_augmented(
      n_vec, z_star_prop, z, phi, data_long, log_prior_n_vec
    )
    
    if (stats::runif(1) < exp(log_post_prop - log_post)) z_star <- z_star_prop
  }

  return (list(log_post = log_post, z_star = z_star))
}



update_beta <- function (
  beta, z, phi, alpha_star, ind_list, X_long, pred_group
) {
  p <- length(pred_group)
  if (p == 0L) return (beta)
  n_risk <- nrow(alpha_star)
  n_regimes <- ncol(alpha_star)  # No. of regimes is no. of change points + 1.
  n_long <- ncol(z)
  z_min_alpha <- matrix(NA_real_, nrow = n_risk, ncol = n_long)
  
  for (k in 1:n_regimes) {
    ind <- ind_list[[n_risk + 1L]][[k]]
    z_min_alpha[, ind] <- z[, ind] - alpha_star[, k]
  }
  
  groups <- unique(pred_group)
  
  # Predictor group inclusion indicators
  inc_mat <- matrix(nrow = length(groups), ncol = n_risk)
  
  for (r in 1:n_risk) for (g in groups) {
    inc_mat[g, r] <- all(beta[pred_group == g, r] == 0)
  }
  
  # Sample the prior inclusion probability with a Beta(1, 1) (= uniform) prior.
  inc_prob <- stats::rbeta(
    n = 1L, shape1 = 1 + sum(inc_mat), shape2 = 1 + sum(!inc_mat)
  )
  
  for (r in 1:n_risk) {
    prec_post_beta_all <- diag(p) / s2_beta + crossprod(X_long / sqrt(phi[r, ]))
    Xt.z <- crossprod(X_long, z_min_alpha[r, ] / phi[r, ])
    inc_vec <- pred_group %in% which(inc_mat[, r])  # Predictor inclusion vector
    p_inc <- sum(inc_vec)  # Number of included predictors
    prec_post_beta <- prec_post_beta_all[inc_vec, inc_vec, drop = FALSE]
    
    mean_post_beta <- if (p_inc == 0L) numeric(0L) else solve(
      a = prec_post_beta, b = Xt.z[inc_vec]
    )
    
    # Metropolis-within-Gibbs for the covariate inclusion indicators that does
    # not condition on `beta`
    for (g in groups) {
      add <- !inc_mat[g, r]  # Whether we are adding a covariate
      inc_vec_prop <- inc_vec
      inc_vec_prop[pred_group == g] <- add
      p_inc_prop <- sum(inc_vec_prop)
      
      prec_post_beta_prop <-
        prec_post_beta_all[inc_vec_prop, inc_vec_prop, drop = FALSE]
      
      mean_post_beta_prop <- if (p_inc_prop == 0L) numeric(0L) else solve(
        a = prec_post_beta_prop, b = Xt.z[inc_vec_prop]
      )
      
      # Equation 12 of doi:10.1214/06-BA105:
      if (stats::runif(1) < exp(.5 * (
        -determinant(prec_post_beta_prop)$modulus + p_inc * log(s2_beta)
          + crossprod(
            mean_post_beta_prop, prec_post_beta_prop %*% mean_post_beta_prop
          ) + determinant(prec_post_beta)$modulus - p_inc_prop * log(s2_beta)
          - crossprod(mean_post_beta, prec_post_beta %*% mean_post_beta)
      ) + (2 * add - 1) * (log(inc_prob) - log1p(-inc_prob)))) {
        inc_vec <- inc_vec_prop
        p_inc <- p_inc_prop
        prec_post_beta <- prec_post_beta_prop
        mean_post_beta <- mean_post_beta_prop
      }
    }
    
    # Sample `beta` conditionally on `inc_vec`.
    beta[, r] <- 0
    if (p_inc == 0L) next
    
    beta[inc_vec, r] <- mvtnorm::rmvnorm(
      n = 1L, mean = mean_post_beta, sigma = solve(prec_post_beta)
    )
  }
  
  return (beta)
}


#' Log of the prior on `n_regimes` up to proportionality
#' 
#'  This implements the Geometric distribution with success probability 0.5 used
#'  in the paper.
#' 
#' @param n_regimes No. of regimes, which is the no. of change points + 1
#' 
#' @keywords internal
log_prior_n_regimes <- function (n_regimes) {
  # Log of the prior on `n_regimes` up to proportionality
  # No. of regimes is no. of change points + 1.
  return (stats::dgeom(x = n_regimes - 1L, prob = 0.5, log = TRUE))
}


get_cp_candidates <- function (data, t_max, cp_constraint) {
  if (!cp_constraint) return (1:t_max)
  event <- as.integer(data$event) - 1L
  cens <- event == 0L
  observed_times <- unique(data$time[!cens])
  
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
  
  return (cp_candidates)
}


# Put a Gaussian prior on alpha.
mu_alpha <- -9
s2_alpha <- 3

s2_beta <- 1  # Gaussian prior on the nonzero beta, i.e. the slab variance


#' Run the multivariate Bernoulli detector
#' 
#' Function that runs the MCMC algorithm for the multivariate Bernoulli detector
#' 
#' @param data Data frame with as first variable an integer vector `time` with
#'  the time-to-event data, as second variable a factor `event` of event
#'  types/causes with the first level indicating that the time was censored, and
#'  optionally additional variables that are treated as predictor values.
#' @param t_max The integer \eqn{t_{\max}}
#' @param n_iter Number of MCMC iterations, including burn-in
#' @param pred_group Group indicator for predictors that should be jointly added
#'  or removed when it comes to variable selection
#' @param burnin Number of burn-in MCMC iterations
#' @param show_progress_bar Whether to show a progress bar
#' @param cp_constraint Whether to limit where change points are allowed in the
#'  prior
#' 
#' @returns A list with the following items:
#' \tabular{ll}{
#'   \code{alpha_MCMC} \tab Array with MCMC samples of alpha \cr
#'   \code{beta_MCMC} \tab Array with MCMC samples of beta \cr
#'   \code{n_regimes_MCMC} \tab Vector with MCMC samples of the no. of regimes,
#'     which is the no. of change points + 1 \cr
#'   \code{log_lik_MCMC} \tab Vector with log-likelihood values of the MCMC
#'     samples. \cr
#'   \code{data} \tab The value of the function argument \code{data} \cr
#'   \code{cp_constraint} \tab The value of the function argument
#'     \code{cp_constraint} \cr
#' }
#' 
#' @export
run_mvb_detector <- function (
  data, t_max = NULL, n_iter = 1e4L, pred_group = NULL, burnin = n_iter %/% 10L,
  show_progress_bar = TRUE, cp_constraint = TRUE
) {
  stopifnot(all(colnames(data[, 1:2]) == c("time", "event")))
  n_risk <- nlevels(data$event) - 1L
  p <- ncol(data) - 2L
  
  # Transform the data to long format.
  if (n_risk > 1L) {
    # The next line also checks whether `data` is in the correct format.
    data_long <- discSurv::dataLongCompRisks(
      dataShort = data, timeColumn = "time", eventColumns = "event",
      eventColumnsAsFactor = TRUE
    )
    
    data_long$y <- 0L
    tmp <- which(data_long$e0 == 0)
    data_long$y[tmp] <- as.integer(data_long$event[tmp]) - 1L
  } else {
    data_tmp <- data
    data_tmp$event <- as.integer(data$event != "censored")
    
    data_long <- discSurv::dataLong(
      dataShort = data_tmp, timeColumn = "time", eventColumn = "event"
    )
  }
  
  X_long <- as.matrix(data_long[, 5L + n_risk + seq_len(p)])
  data_long <- data_long[, c("timeInt", "y")]
  time <- data$time
  event <- as.integer(data$event) - 1L
  cens <- event == 0L
  if (missing(t_max)) t_max <- max(time)
  X <- as.matrix(data[, 2L + seq_len(p)])
  
  stopifnot(all(time <= t_max))
  if (n_risk > 64L) stop("The code only supports `n_risk` <= 64.")
  n <- length(time)
  if (missing(X)) X <- matrix(data = 0, nrow = n, ncol = 0L)
  p <- NCOL(X)
  if (missing(pred_group)) pred_group <- seq_len(p)
  
  if (length(pred_group) != p) stop(
    "The provided `pred_group` does not match the number of columns of `X`."
  )
  
  groups <- unique(pred_group)
  
  if (!all(sort(groups) == seq_along(groups))) stop(
    "The groups in `pred_group` need to be numbered consecutively starting 1."
  )
  
  # Specify the prior on change points:
  # Only allow a start of a regime at or right after an observed time, or time
  # 1. The last time point t_max cannot be the start of a regime, which is a
  # standard restriction in change point modelling, e.g.,
  # doi:10.1080/02664763.2011.559209.
  cp_candidates <- get_cp_candidates(data, t_max, cp_constraint)
  n_cp_candidates <- length(cp_candidates)
  
  # Function that computes the log of the prior on change points
  log_prior_n_vec <- function (n_vec) {
    # Truncated prior that puts zero mass on change points at time points that
    # are not in `cp_candidates`, which is a subset of times if `cp_constraint`
    # is true. `log_prior_n_vec` does not include the prior on cause-specific
    # change points p(z_t | gamma_t).
    n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
    
    if (!all((cumsum(n_vec[-n_regimes]) + 1L) %in% cp_candidates)) return (-Inf)
    
    return (log_prior_n_regimes(n_regimes) - lchoose(
      n_cp_candidates - 1L, n_regimes - 1L
    ))
  }
  
  # Initialize objects to store MCMC results.
  n_recorded <- n_iter - burnin
  alpha_MCMC <- array(data = NA_real_, dim = c(n_recorded, n_risk, t_max))
  beta_MCMC <- array(data = NA_real_, dim = c(n_recorded, p, n_risk))
  n_regimes_MCMC <- integer(n_recorded)
  
  n_regimes_r_MCMC <- matrix(
    data = NA_integer_, nrow = n_recorded, ncol = n_risk
  )
  
  log_lik_MCMC <- numeric(n_recorded)
  
  # Initialize with n_regimes = 1.
  n_vec <- t_max
  
  # Crude estimate of instantaneous death probability
  p_hat <- 1 / mean(time[!cens])
  
  alpha_star <- matrix(
    data = log(p_hat) - log1p(-p_hat) - log(n_risk), nrow = n_risk, ncol = 1L
  )
  
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  z_star <- matrix(TRUE, nrow = n_risk, ncol = n_regimes)
  
  if (n_regimes > 1L) for (k in 2:n_regimes) {
    z_star[, k] <- alpha_star[, k] != alpha_star[, k - 1L]
  }
  
  beta <- matrix(data = stats::rnorm(
    n = p * n_risk, sd = sqrt(s2_beta)
  ), nrow = p, ncol = n_risk)
  
  for (r in 1:n_risk) for (g in groups) {
    if (stats::runif(1L) < 0.5) beta[pred_group == g, r] <- 0
  }
  
  if (show_progress_bar) pb <- utils::txtProgressBar(max = n_iter, style = 3L)
  
  for (s in 1:n_iter) {
    # Local step
    ind_list <- get_ind_list(n_vec, z_star, data_long)
    tmp <- sample_augmentation(alpha_star, beta, ind_list, data_long, X_long)
    z <- tmp$z
    phi <- tmp$phi
    
    z_min_mu <- z - t(X_long %*% beta)
    
    tmp <- update_n_vec_local(
      n_vec, z_star, z_min_mu, phi, data_long, n_risk, t_max, log_prior_n_vec
    )
    
    n_vec <- tmp$n_vec
    z_star <- tmp$z_star
    
    tmp <- update_z_local(
      tmp$log_post, n_vec, z_star, z_min_mu, phi, data_long, n_risk,
      log_prior_n_vec
    )
    
    z_star <- tmp$z_star
    
    ind_list <- get_ind_list(n_vec, z_star, data_long)
    alpha_star <- update_alpha(z_min_mu, phi, z_star, ind_list)
    
    beta <- update_beta(beta, z, phi, alpha_star, ind_list, X_long, pred_group)
    
    # Global step
    tmp <- update_n_vec_global(
      n_vec, alpha_star, z_star, beta, time, event, cens, X, n_risk,
      t_max, log_prior_n_vec, mu_alpha, s2_alpha
    )
    
    n_vec <- tmp$n_vec
    alpha_star <- tmp$alpha_star
    z_star <- tmp$z_star
    
    tmp <- update_z_global(
      tmp$log_post, n_vec, alpha_star, z_star, beta, time, event, cens,
      X, n_risk, log_prior_n_vec, mu_alpha, s2_alpha
    )
    
    alpha_star <- tmp$alpha_star
    z_star <- tmp$z_star
    
    if (s > burnin) {
      alpha_MCMC[s - burnin, , ] <-
        alpha_star[, rep(x = 1:length(n_vec), times = n_vec)]
      
      beta_MCMC[s - burnin, , ] <- beta
      n_regimes_MCMC[s - burnin] <- length(n_vec)
      n_regimes_r_MCMC[s - burnin, ] <- rowSums(z_star)
      
      log_lik_MCMC[s - burnin] <- log_likelihood(
        n_vec, alpha_star, event, time, cens, mu = X %*% beta
      )
    }
    
    if (show_progress_bar) utils::setTxtProgressBar(pb, s)
  }
  
  if (show_progress_bar) close(pb)
  
  return (list(
    alpha_MCMC = alpha_MCMC, beta_MCMC = beta_MCMC,
    n_regimes_MCMC = n_regimes_MCMC, n_regimes_r_MCMC = n_regimes_r_MCMC,
    log_lik_MCMC = log_lik_MCMC, data = data, cp_constraint = cp_constraint
  ))
}


#' Simulate data
#' 
#' Simulate competing risks discrete survival data from the multinomial logit
#' model for the cause-specific hazard functions.
#' 
#' @param n Number of observations
#' @param alpha Matrix with baseline hazard parameters
#' @param X Matrix with predictors
#' @param beta Matrix with cause-specific regression coefficients
#' 
#' @returns Data frame with an integer vector `time` with the time-to-event
#'  data, a factor `event` of event types/causes with the first level indicating
#'  that the time was censored, and optionally additional variables that are the
#'  predictor values provided by the argument `X`.
#' 
#' @export
simulate_data <- function (n, alpha, X = NULL, beta = 0) {
  n_risk <- nrow(alpha)
  t_max <- ncol(alpha)
  if (missing(X)) X <- matrix(data = 0, nrow = n, ncol = 0L)
  if (missing(beta)) beta <- matrix(data = 0, nrow = NCOL(X), ncol = n_risk)
  mu <- X %*% beta
  
  time <- integer(n)
  event_int <- integer(n)  # Record which of the competing events happens
  
  for (i in seq_len(n)) {
    
    for (t in 1:t_max) {
      
      lam <- exp(log_hazard_rate_compete(alpha[, t] + mu[i, ]))
      U <- stats::runif(1)
      
      for (r in 1:n_risk) if ({U <- U - lam[r]} < 0) {
        event_int[i] <- r
        break
      }
      
      if (event_int[i] != 0L) break
      
    }
    
    time[i] <- if (event_int[i] != 0L) t else t_max + 1L
    
  }
  
  cens <- time > t_max
  time[cens] <- t_max
  event <- factor(x = event_int, levels = 0:n_risk)
  levels(event) <- c("censored", 1:n_risk)
  return (data.frame(time = time, event = event, X = X))
}
