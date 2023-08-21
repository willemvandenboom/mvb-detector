### Simulation study without change points

# This code produces the simulation study without change points in the
# appendices of the paper "The Multivariate Bernoulli detector: Change point
# estimation in discrete survival analysis" by Willem van den Boom, Maria De
# Iorio, Fang Qian and Alessandra Guglielmi.


source("mvb_detector.R")



## Simulate data

set.seed(1L)
n_risk <- 3L  # Number of competing risks

# Cause-specific hazard
lam <- exp(log_hazard_rate_compete(c(-2, -3, -4)))

# Generate survival times
n <- 1e2L  # Number of observations
time_to_event <- integer(n)
event <- integer(n)  # Record which of the competing events happened
cens <- logical(n)  # Whether the time is censored


for (i in seq_len(n)) {
  
  for (t in 1:1e2L) {
    
    U <- runif(1)
    
    for (r in 1:n_risk) if ({U <- U - lam[r]} < 0) {
      event[i] <- r
      break
    }
    
    if (event[i] != 0L) break
    
  }
  
  time_to_event[i] <- t
  
}



## Specify the prior on change points.

t_max <- max(time_to_event)
observed_times <- unique(time_to_event[!cens])

# Only allow a start of a regime at or right after an observed time, or time 1.
# The last time point t_max cannot be the start of a regime, which is a standard
# restriction in change point modelling, e.g., doi:10.1080/02664763.2011.559209.
cp_candidates <- unique(c(
  1L, observed_times[observed_times < t_max],
  observed_times[observed_times < t_max - 1L] + 1L
))

# Do not allow a start of a regime at time t right after a time t - 1 with an
# observed event if time t + 1 has no observed events.
if (t_max > 2L) for (t in 2:(t_max - 1L)) {
  if (!(t %in% cp_candidates) | (t %in% observed_times)) next
  
  if (((t - 1L) %in% observed_times) & ((t + 1L) %in% observed_times)) {
    cp_candidates <- cp_candidates[-which(cp_candidates == t)]
  }
}

n_cp_candidates <- length(cp_candidates)


log_prior_n_vec <- function (n_vec) {
  # Truncated prior that puts zero mass on change points at time points with no
  # observed events.
  n_regimes <- length(n_vec)  # No. of regimes is no. of change points + 1.
  if (!all((cumsum(n_vec[-n_regimes]) + 1L) %in% cp_candidates)) return (-Inf)
  
  return (log_prior_n_regimes(n_regimes) - lchoose(
    n_cp_candidates - 1L, n_regimes - 1L
  ))
}



## Transform the data to long format.

# Regression terms
p <- 0L
X <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
beta_true <- matrix(data = 0, nrow = p, ncol = n_risk)

tmp <- discSurv::dataLong(
  dataShort = data.frame(
    time_to_event = time_to_event, eventColumn = ifelse(cens, 0L, 1L),
    event = event, X = X
  ),
  timeColumn = "time_to_event", eventColumn = "eventColumn"
)

X_long <- as.matrix(tmp[, 6L + seq_len(p)])
data_long <- tmp[, c("timeInt", "y", "event")]
data_long$y <- as.integer(data_long$y)
tmp <- which(data_long$y == 1L)
data_long$y[tmp] <- data_long$event[tmp]
data_long$event <- NULL



## Run the MCMC.
n_iter <- 1e5L  # Number of MCMC iterations
res <- run_MCMC(n_iter)



## Compute the Bayes factor.

n_regimes_MCMC <- res$n_regimes_MCMC
n_regimes_post <- tail(n_regimes_MCMC, 0.9 * n_iter)
post_prob_no_cp <- mean(n_regimes_post == 1)

# Compute prior marginal probability of a change point
log_prior_n_regimes_normalizing_constant <- matrixStats::logSumExp(
  log_prior_n_regimes(1:n_cp_candidates)
)

prior_prob_no_cp <- exp(
  log_prior_n_regimes(1L) - log_prior_n_regimes_normalizing_constant
)

print("Bayes factor:")
print(post_prob_no_cp / prior_prob_no_cp)



### Compute the confidence interval for alpha.

# Get rid of 10% of iterations as burnin.
alpha_post <- res$alpha_MCMC[(n_iter %/% 10L + 1L):n_iter, , , drop = FALSE]
alpha_mean <- colMeans(alpha_post)
alpha_CI <- array(NA_real_, dim = c(2, n_risk, t_max))

for (r in 1:n_risk) for (t in 1:t_max) {
  alpha_CI[, r, t] <- quantile(alpha_post[, r, t], c(0.025, 0.975))
}

print("Confidence intervals:")
print(alpha_CI)
