### Multivariate Bernoulli detector

# This code implements the Markov chain Monte Carlo described in the paper
# "The Multivariate Bernoulli detector: Change point estimation in discrete
# survival analysis" by Willem van den Boom, Maria De Iorio, Fang Qian and
# Alessandra Guglielmi.


# Install the required packages.
for (tmp in c(
  "brea", "discSurv", "matrixStats", "mvtnorm", "nnet", "Rcpp", "RcppParallel",
  "R.utils"
)) {
  if(!tmp %in% rownames(installed.packages())) install.packages(
    pkgs = tmp, repos = "https://cloud.r-project.org", dependencies = TRUE
  )
}



## Code for the augmented updates, i.e. local MCMC


# Code that provides indices that assign the discrete survival data in long
# format to their respective regimes which are marked by the change points:

Rcpp::sourceCpp(code="
#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::List get_ind_list_Rcpp(
  const Rcpp::IntegerVector regime, const int n_regimes
) {
  int n_long = regime.length();
  std::vector<std::vector<int> > ind_list(n_regimes);
  for (int i = 0; i < n_long; i++) ind_list[regime[i] - 1].push_back(i + 1);
  Rcpp::List R_list(n_regimes);
  for (int k = 0; k < n_regimes; k++) R_list[k] = ind_list[k];
  return R_list;
}
"
)


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


Rcpp::sourceCpp(code="
// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>


int sample_r_i(
  double U, const double z_min_eta,
  const RcppParallel::RVector<double> m_mix,
  const RcppParallel::RVector<double> tmp1,
  const RcppParallel::RVector<double> tmp2
) {
  int n_mix = m_mix.length();
  std::vector<double> prob(n_mix);
  
  for (int i = 0; i < n_mix; i++) {
    prob[i] = z_min_eta - m_mix[i];
    prob[i] = tmp1[i] - tmp2[i] * prob[i] * prob[i];
  }
  
  double prob_sum = 0.0, prob_max = *max_element(prob.begin(), prob.end());
  
  for (int i = 0; i < n_mix; i++) {
    prob[i] = std::exp(prob[i] - prob_max);
    prob_sum += prob[i];
  }
  
  int r_i = 0;
  
  while (r_i < n_mix) {
    prob[r_i] /= prob_sum;
    if (U < prob[r_i]) break;
    U -= prob[r_i++];
  }
  
  return r_i + 1;
}


struct sample_r_i_worker : public RcppParallel::Worker
{
   // Parameters
   const RcppParallel::RMatrix<double> U;
   const RcppParallel::RMatrix<double> z_min_eta;
   const RcppParallel::RVector<double> m_mix;
   const RcppParallel::RVector<double> tmp1;
   const RcppParallel::RVector<double> tmp2;
   
   // destination matrix
   RcppParallel::RMatrix<int> r_i_mat;
   
   // initialize with parameters and destination
   sample_r_i_worker(
    const Rcpp::NumericMatrix U, const Rcpp::NumericMatrix z_min_eta,
    const Rcpp::NumericVector m_mix, const Rcpp::NumericVector tmp1,
    const Rcpp::NumericVector tmp2, Rcpp::IntegerMatrix r_i_mat
  ) 
    : U(U), z_min_eta(z_min_eta), m_mix(m_mix), tmp1(tmp1), tmp2(tmp2),
      r_i_mat(r_i_mat) {}
   
   void operator()(std::size_t begin, std::size_t end) {
     std::transform(
      U.begin() + begin,
      U.begin() + end,
      z_min_eta.begin() + begin,
      r_i_mat.begin() + begin,
      [=](double U, double z_min_eta){
        return sample_r_i(U, z_min_eta, m_mix, tmp1, tmp2);
      }
     );
  }
};


// [[Rcpp::export]]
Rcpp::IntegerMatrix sample_r_i_mat(
  const Rcpp::NumericMatrix U, const Rcpp::NumericMatrix z_min_eta,
  const Rcpp::NumericVector m_mix, const Rcpp::NumericVector tmp1,
  const Rcpp::NumericVector tmp2
) {
  int n_risk = z_min_eta.nrow(), n_long_k = z_min_eta.ncol();
  Rcpp::IntegerMatrix r_i_mat(n_risk, n_long_k);  // Allocate the output matrix.
  
  // `sample_r_i` functor (Pass input and output matrices.)
  sample_r_i_worker worker(U, z_min_eta, m_mix, tmp1, tmp2, r_i_mat);
  
  // Call `parallelFor` to do the work.
  RcppParallel::parallelFor(0, U.length(), worker);
  
  return r_i_mat;  // Return the output matrix.
}
"
)


sample_augmentation <- function (alpha_star, beta, ind_list, data_long) {
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
    U_i <- runif(n_long_k)
    colSums_lam <- colSums(lam)
    
    for (r in 1:n_risk) z[r, ind] <- -log(-log(U_i) / (1 + colSums_lam) - log(
      runif(n_long_k)
    ) / lam[r, ] * (data_long$y[ind] != r))
    
    # Table 1 of doi:10.1016/j.csda.2006.10.006
    w_mix <- c(.00397, .0396, .168, .147, .125, .101, .104, .116, .107, .088)
    m_mix <- c(5.09, 3.29, 1.82, 1.24, .764, .391, .0431, -.306, -.673, -1.06)
    phi_mix <- c(4.50, 2.02, 1.10, .422, .198, .107, .0778, .0766, .0947, .146)
    
    tmp1 <- log(w_mix) - .5 * log(phi_mix)
    tmp2 <- .5 / phi_mix
    U <- matrix(data = runif(n_risk * n_long_k), nrow = n_risk, ncol = n_long_k)
    
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
    
    alpha_star[r, k] <- mean_post_alpha + rnorm(1L) / sqrt(prec_post_alpha)
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


log_posterior_augmented <- function (n_vec, z_star, z, phi, data_long) {
  ind_list <- get_ind_list(n_vec, z_star, data_long)
  
  return (
    log_prior_n_vec(n_vec) + compute_log_marginal_likelihood(ind_list, z, phi)
  )
}


sample_z <- function () {
  n_z <- as.integer(2^n_risk) - 1L
  z <- rep(FALSE, n_risk)
  k <- sample.int(n = n_z, size = 1L)  # Uniform distribution
  
  z[min(n_risk - ceiling(log2(k + 1L) - 1L), n_risk):n_risk] <- as.logical(
    as.numeric(strsplit(R.utils::intToBin(k), split = "")[[1]])
  )
  
  return (z)
}


update_n_vec_local <- function (n_vec, z_star, z, phi, data_long) {
  q <- 0.5  # Split proposal probability
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  log_post <- log_posterior_augmented(n_vec, z_star, z, phi, data_long)
  
  if (n_regimes == 1 || (n_regimes < t_max && runif(1) < q)) {  # Split
    splittable <- (1:n_regimes)[n_vec > 1L]
    
    if (length(splittable) > 1L) {
      j <- sample(x = (1:n_regimes)[n_vec > 1L], size = 1L)
    } else {
      j <- splittable[1]
    }
    
    l <- sample.int(n = n_vec[j] - 1L, size = 1L)
    
    n_vec_prop <- c(
      head(n_vec, j - 1L), l, n_vec[j] - l, tail(n_vec, n_regimes - j)
    )
    
    z_new <- if (j == 1L) rep(TRUE, n_risk) else sample_z()
    
    z_star_prop <- cbind(
      z_star[, seq_len(j), drop = FALSE],
      z_new,
      z_star[, j + seq_len(n_regimes - j), drop = FALSE]
    )
    
    log_post_prop <- log_posterior_augmented(
      n_vec_prop, z_star_prop, z, phi, data_long
    )
    
    if (runif(1) < exp(
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
      head(n_vec, j - 1L), n_merge, tail(n_vec, n_regimes - j - 1L)
    )
    
    z_star_prop <- cbind(
      z_star[, seq_len(j), drop = FALSE],
      z_star[, j + 1L + seq_len(n_regimes - j - 1L), drop = FALSE]
    )
    
    log_post_prop <- log_posterior_augmented(
      n_vec_prop, z_star_prop, z, phi, data_long
    )
    
    if (runif(1) < exp(
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
      n_vec_prop, z_star, z, phi, data_long
    )
    
    if (runif(1) < exp(log_post_prop - log_post)) {
      n_vec <- n_vec_prop
      log_post <- log_post_prop
    }
  }
  
  return (list(log_post = log_post, n_vec = n_vec, z_star = z_star))
}


update_z_local <- function (log_post, n_vec, z_star, z, phi, data_long) {
  # Resample a randomly selected z_star.
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.

  if (n_regimes > 1L) {
    j <- if (n_regimes == 2L) 2L else sample(x = 2:n_regimes, size = 1L)
    z_new <- sample_z()
    
    if (all(z_new == z_star[, j, drop = FALSE])) return (
      list(log_post = log_post, z_star = z_star)
    )

    z_star_prop <- z_star
    z_star_prop[, j] <- z_new
    
    log_post_prop <- log_posterior_augmented(
      n_vec, z_star_prop, z, phi, data_long
    )
    
    if (runif(1) < exp(log_post_prop - log_post)) z_star <- z_star_prop
  }

  return (list(log_post = log_post, z_star = z_star))
}



update_beta <- function (beta, z, phi, alpha_star, ind_list) {
  if (p == 0L) return (beta)
  
  # Sample the prior inclusion probability with a Beta(1, 1) prior.
  inc_prob <- rbeta(
    n = 1, shape1 = 1 + sum(beta != 0), shape2 = 1 + sum(beta == 0)
  )
  
  n_risk <- nrow(alpha_star)
  n_regimes <- ncol(alpha_star)  # No. of regimes is no. of change points + 1.
  n_long <- ncol(z)
  z_min_alpha <- matrix(NA_real_, nrow = n_risk, ncol = n_long)
  
  for (k in 1:n_regimes) {
    ind <- ind_list[[n_risk + 1L]][[k]]
    z_min_alpha[, ind] <- z[, ind] - alpha_star[, k]
  }
  
  for (r in 1:n_risk) {
    inc_vec <- beta[, r] != 0  # Covariate inclusion vector
    p_inc <- sum(inc_vec)  # Number of included covariates
    X_long_inc <- X_long[, inc_vec, drop = FALSE]
    
    prec_post_beta <- diag(p_inc) / s2_beta + crossprod(
      X_long_inc / phi[r, ], X_long_inc
    )
    
    mean_post_beta <- if (p_inc == 0L) numeric(0) else solve(
      a = prec_post_beta, b = crossprod(X_long_inc, z_min_alpha[r, ] / phi[r, ])
    )
    
    # Metropolis-within-Gibbs for the covariate inclusion indicators that does
    # not condition on `beta`
    for (j in 1:p) {
      inc_vec_prop <- inc_vec
      inc_vec_prop[j] <- !inc_vec[j]
      p_inc_prop <- sum(inc_vec_prop)
      X_long_inc_prop <- X_long[, inc_vec_prop, drop = FALSE]
      
      prec_post_beta_prop <- diag(p_inc_prop) / s2_beta + crossprod(
        X_long_inc_prop / phi[r, ], X_long_inc_prop
      )
      
      mean_post_beta_prop <- if (p_inc_prop == 0L) numeric(0) else solve(
        a = prec_post_beta_prop,
        b = crossprod(X_long_inc_prop, z_min_alpha[r, ] / phi[r, ])
      )
      
      # Equation 12 of doi:10.1214/06-BA105:
      if (runif(1) < exp(.5 * (
        -determinant(prec_post_beta_prop)$modulus + p_inc * log(s2_beta)
          + crossprod(
            mean_post_beta_prop, prec_post_beta_prop %*% mean_post_beta_prop
          ) + determinant(prec_post_beta)$modulus - p_inc_prop * log(s2_beta)
          - crossprod(mean_post_beta, prec_post_beta %*% mean_post_beta)
      ) + (p_inc_prop - p_inc) * (log(inc_prob) - log1p(-inc_prob)))) {
        inc_vec <- inc_vec_prop
        p_inc <- p_inc_prop
        X_long_inc <- X_long_inc_prop
        prec_post_beta <- prec_post_beta_prop
        mean_post_beta <- mean_post_beta_prop
      }
    }
    
    # Sample `beta` conditional on `inc_vec`.
    beta[, r] <- 0
    if (p_inc == 0L) next
    
    beta[inc_vec, r] <- mvtnorm::rmvnorm(
      n = 1, mean = mean_post_beta, sigma = solve(prec_post_beta)
    )
  }
  
  return (beta)
}



## Code for the unaugmented updates, i.e. global MCMC


Rcpp::sourceCpp(code="
#include <Rcpp.h>
#include <cmath>       /* log1p */


// [[Rcpp::export]]
Rcpp::NumericVector log_hazard_rate_compete(Rcpp::NumericVector alpha) {
  // This function assumes that `alpha` is an `n_risk`-dimensional vector.
  return alpha - std::log1p(Rcpp::sum(Rcpp::exp(alpha)));
}


// [[Rcpp::export]]
double log_likelihood(
  Rcpp::IntegerVector n_vec, Rcpp::NumericMatrix alpha_star,
  Rcpp::IntegerVector event, Rcpp::IntegerVector time_to_event,
  Rcpp::LogicalVector cens, Rcpp::NumericMatrix mu
) {
  int n_regimes = n_vec.length(), n = event.length();
  Rcpp::IntegerVector n_vec_sum = Rcpp::cumsum(n_vec);
  double log_lik = 0.0;
  
  for (int i = 0; i < n; i++) for (int k = 0; k < n_regimes; k++) {
    Rcpp::NumericVector log_rate = log_hazard_rate_compete(
      alpha_star(Rcpp::_, k) + mu(i, Rcpp::_)
    );
    
    double log_surv_prob = std::log1p(-Rcpp::sum(Rcpp::exp(log_rate)));
    
    if (time_to_event[i] <= n_vec_sum[k]) {
      log_lik += (
        time_to_event[i] - (k == 0 ? 0 : n_vec_sum[k - 1]) - !cens[i]
      ) * log_surv_prob;
      
      if (not cens[i]) log_lik += log_rate[event[i] - 1];
      break;  // Go to next individual.
    }
    
    log_lik += n_vec[k] * log_surv_prob;
  }
  
  return log_lik;
}
"
)


log_prior_alpha <- function (alpha_star) {
  n_risk <- nrow(alpha_star)
  ret <- 0
  
  for (r in 1:n_risk) ret <- ret + sum(dnorm(
    x = unique(alpha_star[r, ]), mean = mu_alpha, sd = sqrt(s2_alpha),
    log = TRUE
  ))
  
  return (ret)
}


log_posterior_unaugmented <- function (n_vec, alpha_star, beta, X) {
  # Log of posterior for discrete-time survival with `n_risk` risks.
  return (log_prior_n_vec(n_vec) + log_prior_alpha(alpha_star) + log_likelihood(
    n_vec, alpha_star, event, time_to_event, cens, mu = X %*% beta
  ))
}


# Standard deviation of the Gaussian random walk for new alpha
s_prop_alpha <- 1


update_n_vec_global <- function (n_vec, alpha_star, z_star, beta) {
  q <- 0.5  # Split proposal probability
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  log_post <- log_posterior_unaugmented(n_vec, alpha_star, beta, X)
  
  if (n_regimes == 1 || (n_regimes < t_max && runif(1) < q)) {  # Split
    splittable <- (1:n_regimes)[n_vec > 1L]
    
    if (length(splittable) > 1L) {
      j <- sample(x = (1:n_regimes)[n_vec > 1L], size = 1L)
    } else {
      j <- splittable[1]
    }
    
    l <- sample.int(n = n_vec[j] - 1L, size = 1L)
    
    n_vec_prop <- c(
      head(n_vec, j - 1L), l, n_vec[j] - l, tail(n_vec, n_regimes - j)
    )
    
    z_new <- if (j == 1L) rep(TRUE, n_risk) else sample_z()
    alpha_new <- alpha_star[, j]
    
    alpha_new[z_new] <- rnorm(
      n = sum(z_new), mean = alpha_star[z_new, j, drop = FALSE],
      sd = s_prop_alpha
    )
    
    alpha_star_prop <- cbind(
      alpha_star[, seq_len(j), drop = FALSE], alpha_new,
      alpha_star[, j + seq_len(n_regimes - j), drop = FALSE]
    )
    
    for (r in 1:n_risk) {
      tmp <- z_star[r, j + seq_len(n_regimes - j), drop = FALSE]
      ind <- ifelse(!any(tmp), n_regimes - j, which.max(tmp) - 1L)
      for (t in seq_len(ind)) alpha_star_prop[r, j + 1L + t] <- alpha_new[r]
    }
    
    log_post_prop <- log_posterior_unaugmented(
      n_vec_prop, alpha_star_prop, beta, X
    )
    
    if (runif(1) < exp(
      (n_regimes < t_max - 1L) * log(1 - q) - (n_regimes > 1) * log(q) +
        log(length(splittable)) + log(n_vec[j] - 1L) - log(n_regimes) +
        log_post_prop - log_post - sum(dnorm(
          x = alpha_new[z_new], mean = alpha_star[z_new, j, drop = FALSE],
          sd = s_prop_alpha, log = TRUE
        ))
    )) {
      log_post <- log_post_prop
      n_vec <- n_vec_prop
      alpha_star <- alpha_star_prop
      z_star <- cbind(
        z_star[, seq_len(j), drop = FALSE], z_new,
        z_star[, j + seq_len(n_regimes - j), drop = FALSE]
      )
    }
  } else {  # Merge
    j <- sample.int(n = n_regimes - 1L, size = 1L)
    n_merge <- n_vec[j] + n_vec[j + 1L]
    
    n_vec_prop <- c(
      head(n_vec, j - 1L), n_merge, tail(n_vec, n_regimes - j - 1L)
    )
    
    alpha_star_prop <- cbind(
      alpha_star[, seq_len(j), drop = FALSE],
      alpha_star[, j + 1L + seq_len(n_regimes - j - 1L), drop = FALSE]
    )
    
    for (r in 1:n_risk) {
      tmp <- z_star[r, j + 1L + seq_len(n_regimes - j - 1L), drop = FALSE]
      ind <- ifelse(!any(tmp), n_regimes - j - 1L, which.max(tmp) - 1L)
      for (t in seq_len(ind)) alpha_star_prop[r, j + t] <- alpha_star[r, j]
    }
    
    log_post_prop <- log_posterior_unaugmented(
      n_vec_prop, alpha_star_prop, beta, X
    )
    
    if (runif(1) < exp(
      (n_regimes > 2L) * log(q) - (n_regimes < t_max) * log(1 - q) +
        log(n_regimes - 1L) - log(sum(n_vec_prop > 1L)) - log(n_merge - 1L) +
        log_post_prop - log_post + sum(dnorm(
          x = alpha_star[z_star[, j + 1L, drop = FALSE], j + 1L],
          mean = alpha_star_prop[z_star[, j + 1L, drop = FALSE], j],
          sd = s_prop_alpha, log = TRUE
        ))
    )) {
      log_post <- log_post_prop
      n_vec <- n_vec_prop
      alpha_star <- alpha_star_prop
      
      z_star <- cbind(
        z_star[, seq_len(j), drop = FALSE],
        z_star[, j + 1L + seq_len(n_regimes - j - 1L), drop = FALSE]
      )
    }
  }
  
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  
  if (n_regimes > 1) {  # Shuffle
    i <- sample.int(n = n_regimes - 1L, size = 1L)
    j <- sample.int(n = n_vec[i] + n_vec[i + 1L] - 1L, size = 1L)
    n_vec_prop <- n_vec
    n_vec_prop[i + 1L] <- n_vec[i] + n_vec[i + 1L] - j
    n_vec_prop[i] <- j
    log_post_prop <- log_posterior_unaugmented(n_vec_prop, alpha_star, beta, X)
    
    if (runif(1) < exp(log_post_prop - log_post)) {
      log_post <- log_post_prop
      n_vec <- n_vec_prop
    }
  }
  
  return (list(
    log_post=log_post, n_vec = n_vec, alpha_star = alpha_star, z_star = z_star
  ))
}


update_z_global <- function (log_post, n_vec, alpha_star, z_star, beta) {
  # Resample a randomly selected z_star.
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  
  if (n_regimes > 1L) {
    j <- if (n_regimes == 2L) 2L else sample(x = 2:n_regimes, size = 1L)
    z_new <- sample_z()
    
    if (all(z_new == z_star[, j, drop = FALSE])) {
      return (list(alpha_star = alpha_star, z_star = z_star))
    }
    
    alpha_star_prop <- alpha_star
    log_prop_alpha <- 0
    
    for (r in 1:n_risk) {
      
      if (z_new[r] && !z_star[r, j]) {
        
        alpha_new <- rnorm(n = 1L, mean = alpha_star[r, j], sd = s_prop_alpha)
        alpha_star_prop[r, j] <- alpha_new
        
        log_prop_alpha <- log_prop_alpha + dnorm(
          x = alpha_new, mean = alpha_star[r, j], sd = s_prop_alpha, log = TRUE
        )
        
      } else if (!z_new[r] && z_star[r, j]) {
        
        alpha_star_prop[r, j] <- alpha_star[r, j - 1L]
        
        log_prop_alpha <- log_prop_alpha - dnorm(
          x = alpha_star[r, j], mean = alpha_star_prop[r, j], sd = s_prop_alpha,
          log = TRUE
        )
        
      }
      
      tmp <- z_star[r, j + seq_len(n_regimes - j), drop = FALSE]
      ind <- ifelse(!any(tmp), n_regimes - j, which.max(tmp) - 1L)
      for (t in seq_len(ind)) alpha_star_prop[r, j + t] <- alpha_star_prop[r, j]
      
    }
    
    log_post_prop <- log_posterior_unaugmented(n_vec, alpha_star_prop, beta, X)
    
    if (runif(1) < exp(log_post_prop - log_post - log_prop_alpha)) {
      log_post <- log_post_prop
      z_star[, j] <- z_new
      alpha_star <- alpha_star_prop
    }
    
  }
  
  return (list(alpha_star = alpha_star, z_star = z_star))
}


log_prior_n_regimes <- function (n_regimes) {
  # Log of the prior on `n_regimes` up to proportionality
  # No. of regimes is no. of change points + 1.
  return (dgeom(x = n_regimes - 1L, prob = 0.5, log = TRUE))
}


# Put a Gaussian prior on alpha.
mu_alpha <- -9
s2_alpha <- 3

s2_beta <- 1  # Gaussian prior on the nonzero beta, i.e. the slab variance


run_MCMC <- function (n_iter = 1e4L, initialize_at_true = FALSE) {
  # Function that runs the main MCMC algorithm
  
  alpha_MCMC <- array(data = NA_real_, dim = c(n_iter, n_risk, t_max))
  beta_MCMC <- array(data = NA_real_, dim = c(n_iter, p, n_risk))
  n_regimes_MCMC <- integer(n_iter)
  log_lik_MCMC <- numeric(n_iter)
  
  if (initialize_at_true) {
    n_vec <- n_vec_true
    alpha_star <- alpha_star_true
  } else {
    # Initialize with n_regimes = 1.
    n_vec <- t_max
    
    # Crude estimate of instantaneous death probability
    p_hat <- 1 / mean(time_to_event[!cens])
    
    alpha_star <- matrix(
      data = log(p_hat) - log1p(-p_hat) - log(n_risk), nrow = n_risk, ncol = 1L
    )
  }
  
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  z_star <- matrix(TRUE, nrow = n_risk, ncol = n_regimes)
  
  if (n_regimes > 1L) for (k in 2:n_regimes) {
    z_star[, k] <- alpha_star[, k] != alpha_star[, k - 1L]
  }
  
  beta <- matrix(data = 0.5, nrow = p, ncol = n_risk)
  pb <- txtProgressBar(max = n_iter, style = 3L)
  
  for (s in 1:n_iter) {
    # Local step
    ind_list <- get_ind_list(n_vec, z_star, data_long)
    tmp <- sample_augmentation(alpha_star, beta, ind_list, data_long)
    z <- tmp$z
    phi <- tmp$phi
    
    z_min_mu <- z - t(X_long %*% beta)
    
    tmp <- update_n_vec_local(n_vec, z_star, z_min_mu, phi, data_long)
    n_vec <- tmp$n_vec
    z_star <- tmp$z_star
    tmp <- update_z_local(tmp$log_post, n_vec, z_star, z_min_mu, phi, data_long)
    z_star <- tmp$z_star
    
    ind_list <- get_ind_list(n_vec, z_star, data_long)
    alpha_star <- update_alpha(z_min_mu, phi, z_star, ind_list)
    
    beta <- update_beta(beta, z, phi, alpha_star, ind_list)
    
    # Global step
    tmp <- update_n_vec_global(n_vec, alpha_star, z_star, beta)
    n_vec <- tmp$n_vec
    alpha_star <- tmp$alpha_star
    z_star <- tmp$z_star
    
    tmp <- update_z_global(tmp$log_post, n_vec, alpha_star, z_star, beta)
    alpha_star <- tmp$alpha_star
    z_star <- tmp$z_star
    
    alpha_MCMC[s, , ] <- alpha_star[, rep(x = 1:length(n_vec), times = n_vec)]
    beta_MCMC[s, , ] <- beta
    n_regimes_MCMC[s] <- length(n_vec)
    
    log_lik_MCMC[s] <- log_likelihood(
      n_vec, alpha_star, event, time_to_event, cens, mu = X %*% beta
    )
    
    if(!isTRUE(getOption('knitr.in.progress'))) setTxtProgressBar(pb, s)
  }
  
  close(pb)
  
  return (list(
    alpha_MCMC = alpha_MCMC, beta_MCMC = beta_MCMC,
    n_regimes_MCMC = n_regimes_MCMC, log_lik_MCMC = log_lik_MCMC
  ))
}
