### Simulation study without change points

# This code produces the simulation study without change points in Supplementary
# Materials of the paper "The Multivariate Bernoulli detector: Change point
# estimation in discrete survival analysis" by Willem van den Boom, Maria De
# Iorio, Fang Qian and Alessandra Guglielmi.


## Simulate data.

set.seed(1L)
n_risk <- 3L  # Number of competing risks

# Generate survival times
data <- mvb.detector::simulate_data(
  n = 200L, alpha = matrix(data = c(-2, -3, -4), nrow = n_risk, ncol = 100L)
)


## Run the MCMC.
mvbd_fit <- mvb.detector::run_mvb_detector(data = data, n_iter = 1e5L)


## Compute the Bayes factor.
t_max <- dim(mvbd_fit$alpha_MCMC)[3]

cp_candidates <- mvb.detector:::get_cp_candidates(
  data, t_max, mvbd_fit$cp_constraint
)

n_cp_candidates <- length(cp_candidates)

prior_prob_no_cp <- exp(
  mvb.detector:::log_prior_n_regimes(1L) - mvb.detector:::logSumExp(
    mvb.detector:::log_prior_n_regimes(1:n_cp_candidates)
  )
)

post_prob <- mean(mvbd_fit$n_regimes_MCMC == 1)
print("Posterior probability:")
print(post_prob)

print("Bayes factor:")
print(post_prob / prior_prob_no_cp)


## Compute the confidence interval for alpha.
alpha_CI <- array(NA_real_, dim = c(2, n_risk, t_max))

for (r in 1:n_risk) for (t in 1:t_max) {
  alpha_CI[, r, t] <- quantile(mvbd_fit$alpha_MCMC[, r, t], c(0.025, 0.975))
}

print("Confidence intervals:")
print(alpha_CI)
for (r in 1:n_risk) print(range(alpha_CI[, r, ]))
