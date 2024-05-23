### Simulation study on restricting change point locations

# This code produces the simulation study on restricting change point locations
# in Supplementary Materials of the paper "The Multivariate Bernoulli detector:
# Change point estimation in discrete survival analysis" by Willem van den Boom,
# Maria De Iorio, Fang Qian and Alessandra Guglielmi.


n_risk <- 1L  # Number of competing risks
n_vec_true <- c(5L, 7L, 8L)  # Specify change points.
t_max <- sum(n_vec_true)  # Number of time points

# No. of regimes is no. of change points + 1.
n_regimes_true <- length(n_vec_true)

alpha_star_true <- matrix(NA_real_, nrow = n_risk, ncol = n_regimes_true)
alpha_star_true[1, ] <- c(-4, -9, -2)
alpha_true <- matrix(NA_real_, nrow = n_risk, ncol = t_max)
for (r in 1:n_risk) alpha_true[r, ] <- rep(alpha_star_true[r, ], n_vec_true)

# Seeds 22 and 35 are chosen to have simulated data with no observations near
# the first change point.
seed_vec <- c(22L, 35L, 1L, 2L)
pdf(file = "simulation_change_point_location.pdf", width = 7, height = 7)
par(mfcol = c(3L, length(seed_vec) %/% 2L), mai = c(0.7, 0.5, 0.3, 0.5))

for (seed in seed_vec) {
  set.seed(seed)
  data <- mvb.detector::simulate_data(100L, alpha_true)
  
  midpoints <- barplot(
    height = hist(
      x = data$time[data$event != "censored"], breaks = 0:t_max + 0.5,
      plot = FALSE
    )$counts,
    names.arg = 1:t_max, main = "Histogram",
    xlab = expression(paste("Time ", italic(t))), ylab = "Frequency",
    ylim = c(0, 15)
  )
  
  for (t in cumsum(n_vec_true[-n_regimes_true])) {
    abline(v = mean(midpoints[t + 0:1]), lty = "dashed", col = "#d62728")
  }
  
  plot_res <- function (cp_constraint) {
    mvbd_fit <- mvb.detector::run_mvb_detector(
      data = data, n_iter = 1e5L, cp_constraint = cp_constraint
    )

    mvb.detector::plot_change_points(
      mvbd_fit = mvbd_fit, n_vec_true = n_vec_true,
      main_vec = paste(
        if (cp_constraint) "Restricted" else "Unrestricted", "change points"
      ), BF_max = 50
    )
  }

  plot_res(TRUE)
  plot_res(FALSE)
}

dev.off()
