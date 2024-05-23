### Simulation study with replicates

# This code produces the simulation study with replicates in Supplementary
# Materials of the paper "The Multivariate Bernoulli detector: Change point
# estimation in discrete survival analysis" by Willem van den Boom, Maria De
# Iorio, Fang Qian and Alessandra Guglielmi.


source("simulation_util.R")


generate_data <- function (alpha_star, n, p) {
  n_risk <- nrow(alpha_star)
  alpha <- matrix(NA_real_, nrow = n_risk, ncol = t_max)
  for (r in 1:n_risk) alpha[r, ] <- rep(alpha_star[r, ], n_vec_true)
  
  beta <- matrix(data = sample(
    x = c(-0.1, 0, 0.1), size = p * n_risk, replace = TRUE
  ), nrow = p, ncol = n_risk)
  
  X <- matrix(data = as.integer(runif(n * p) < 0.5), nrow = n, ncol = p)
  data <- mvb.detector::simulate_data(n, alpha, X, beta)
  
  # Censoring time *if* the corresponding individual is censored
  cens_time <- integer(n)
  
  for (i in 1:n) {
    cens_time[i] <- if (with(data, time[i] == 1 && event[i] != "censored")) {
      NA_integer_
    } else sample.int(
      n = with(data, time[i] - (event[i] != "censored")), size = 1L
    )
  }
  
  return (list(alpha = alpha, beta = beta, data = data, cens_time = cens_time))
}


## Simulate data

# Set change points.
n_vec_true <- c(5L, 7L, 8L)
t_max <- sum(n_vec_true)  # Number of time points

# No. of regimes is no. of change points + 1.
n_regimes_true <- length(n_vec_true)

max_n_risk <- 4L
alpha_star_true <- matrix(NA_real_, nrow = max_n_risk, ncol = n_regimes_true)
alpha_star_true[1, ] <- c(-5, -4, -2)
alpha_star_true[2, ] <- c(-5, -3, -3)
alpha_star_true[3, ] <- c(-5, -3, -2)
alpha_star_true[4, ] <- c(-4, -4, -4)
n_regimes_r_true <- c(3L, 2L, 3L, 1L)

n_vec <- c(200L, 500L, 1000L)
cens_level_vec <- c(0, .1, .5)
n_risk_vec <- 2:max_n_risk

# Generate survival times
n_rep <- 128L
set.seed(1L)

replica_list <- replicate(n = n_rep, expr = generate_data(
  alpha_star = alpha_star_true[1:2, ], n = max(n_vec), p = 3L
), simplify = FALSE)

replica_list_p5 <- replicate(n = n_rep, expr = generate_data(
  alpha_star = alpha_star_true[1:2, ], n = n_vec[2], p = 5L
), simplify = FALSE)

replica_list_n_risk <- c(list(replica_list), lapply(
  X = n_risk_vec[-1],
  FUN = function (n_risk) replicate(n = n_rep, expr = generate_data(
    alpha_star = alpha_star_true[1:n_risk, ], n = max(n_vec), p = 3L
  ), simplify = FALSE)
))

RcppParallel::setThreadOptions(numThreads = 1L)  # Avoid nested parallelization
n <- n_vec[2]
n_cens_level <- length(cens_level_vec)
result_list_cens <- vector(mode = "list", length = n_cens_level)

for (cens_ind in 1:n_cens_level) {
  cens_level <- cens_level_vec[cens_ind]
  print(paste("Working on censoring level", cens_level, "..."))
  cl <- parallel::makeForkCluster(parallel::detectCores(logical = FALSE))
  
  result_list_cens[[cens_ind]] <- pbapply::pblapply(
    X = 1:n_rep, FUN = function (i) run_inference(
      replica_list[[i]]$data[1:n, ], cens_level,
      replica_list[[i]]$cens_time[1:n]
    ), cl = cl
  )
  
  parallel::stopCluster(cl)
}

n_n <- length(n_vec)
result_list_n <- vector(mode = "list", length = n_n)

for (n_ind in 1:n_n) {
  if (n_ind == 2L) {
    result_list_n[[n_ind]] <- result_list_cens[[2]]
    next
  }
  
  n <- n_vec[n_ind]
  print(paste("Working on n =", n, "..."))
  cl <- parallel::makeForkCluster(parallel::detectCores(logical = FALSE))
  
  result_list_n[[n_ind]] <- pbapply::pblapply(
    X = 1:n_rep, FUN = function (i) run_inference(
      replica_list[[i]]$data[1:n, ], cens_level_vec[2],
      replica_list[[i]]$cens_time[1:n]
    ), cl = cl
  )
  
  parallel::stopCluster(cl)
}

result_list_p <- list(result_list_cens[[2]], NULL)
p <- 5L
print(paste("Working on p =", p, "..."))
cl <- parallel::makeForkCluster(parallel::detectCores(logical = FALSE))

result_list_p[[2]] <- pbapply::pblapply(
  X = 1:n_rep, FUN = function (i) run_inference(
    replica_list_p5[[i]]$data, cens_level_vec[2], replica_list_p5[[i]]$cens_time
  ), cl = cl
)

parallel::stopCluster(cl)

n_n_risk <- length(n_risk_vec)
print("Varying `n_risk` with n = 500...")
result_list_n_risk <- vector(mode = "list", length = n_n_risk)

for (n_risk_ind in 1:n_n_risk) {
  if (n_risk_ind == 1L) {
    result_list_n_risk[[n_risk_ind]] <- result_list_cens[[2]]
    next
  }
  
  n_risk <- n_risk_vec[n_risk_ind]
  print(paste("Working on n_risk =", n_risk, "..."))
  cl <- parallel::makeForkCluster(parallel::detectCores(logical = FALSE))
  
  result_list_n_risk[[n_risk_ind]] <- pbapply::pblapply(
    X = 1:n_rep, FUN = function (i) run_inference(
      replica_list_n_risk[[n_risk_ind]][[i]]$data[1:n_vec[2], ],
      cens_level_vec[2],
      replica_list_n_risk[[n_risk_ind]][[i]]$cens_time[1:n_vec[2]]
    ), cl = cl
  )
  
  parallel::stopCluster(cl)
}

print("Varying `n_risk` with n = 1000...")
result_list_n_risk_n1000 <- vector(mode = "list", length = n_n_risk)

for (n_risk_ind in 1:n_n_risk) {
  if (n_risk_ind == 1L) {
    result_list_n_risk_n1000[[n_risk_ind]] <- result_list_n[[3]]
    next
  }
  
  n_risk <- n_risk_vec[n_risk_ind]
  print(paste("Working on n_risk =", n_risk, "..."))
  cl <- parallel::makeForkCluster(parallel::detectCores(logical = FALSE))
  
  result_list_n_risk_n1000[[n_risk_ind]] <- pbapply::pblapply(
    X = 1:n_rep, FUN = function (i) run_inference(
      replica_list_n_risk[[n_risk_ind]][[i]]$data, cens_level_vec[2],
      replica_list_n_risk[[n_risk_ind]][[i]]$cens_time
    ), cl = cl
  )
  
  parallel::stopCluster(cl)
}


methods <- c("mvbd", "nnet", "brea")
n_methods <- length(methods)


compute_bias_MSE <- function (result_list, replica_list) {
  n_setups <- length(result_list)

  if (length(replica_list[[1]]) == 4L) {
    if (n_rep == 4L) stop("Structure of `replica_list` is ambiguous.")
    replica_list <- rep(list(replica_list), n_setups)
  }

  MSE <- array(
    data = NA_real_, dim = c(n_setups, n_rep, 2L * n_methods + 2L),
    dimnames = list(NULL, NULL, c(
      paste0("alpha_", methods), paste0("beta_", methods), "n_regimes",
      "n_regimes_r"
    ))
  )

  bias <- MSE

  for (ind in 1:n_setups) for (i in 1:n_rep) {
    replica <- replica_list[[ind]][[i]]
    result <- result_list[[ind]][[i]]
    n_risk <- nrow(replica$alpha)

    # Inference on baseline hazards:
    for (alpha in c("alpha_mvbd", "alpha_nnet", "alpha_brea")) {
      error <- result[[alpha]] - replica$alpha
      bias[ind, i, alpha] <- mean(error)
      MSE[ind, i, alpha] <- mean(error^2)
    }

    # Inference on regression coefficients:
    for (beta in c("beta_mvbd", "beta_nnet", "beta_brea")) {
      error <- result[[beta]] - replica$beta
      bias[ind, i, beta] <- mean(error)
      MSE[ind, i, beta] <- mean(error^2)
    }

    error <- mean(result$mvbd_fit$n_regimes_MCMC) - n_regimes_true
    bias[ind, i, "n_regimes"] <- error
    MSE[ind, i, "n_regimes"] <- error^2

    error <-
      colMeans(result$mvbd_fit$n_regimes_r_MCMC) - n_regimes_r_true[1:n_risk]

    bias[ind, i, "n_regimes_r"] <- mean(error)
    MSE[ind, i, "n_regimes_r"] <- mean(error^2)
  }

  return (list(bias = bias, MSE = MSE))
}


get_MSE <- function (result_list, replica_list, ind_vec, outcome) {
  compute_bias_MSE(result_list, replica_list)[[outcome]][, , ind_vec]
}


for (outcome in c("bias", "MSE")) for (param in c("alpha", "beta", "K")) {
  ind_vec <- switch(
    EXPR = param, alpha = 1:n_methods, beta = n_methods + 1:n_methods,
    K = 2L * n_methods + 1:2
  )
  
  MSE_arr <- abind::abind(
    get_MSE(result_list_cens, replica_list, ind_vec, outcome),
    get_MSE(result_list_n, replica_list, ind_vec, outcome),
    get_MSE(
      result_list_p, list(replica_list, replica_list_p5), ind_vec, outcome
    ),
    get_MSE(result_list_n_risk, replica_list_n_risk, ind_vec, outcome),
    get_MSE(result_list_n_risk_n1000, replica_list_n_risk, ind_vec, outcome),
    along = 1L
  )
  
  log_scale <- outcome == "MSE"
  n_setups <- nrow(MSE_arr)
  pch_vec <- c(19L, 4L, 2L)
  lwd <- 1.5
  col_vec <- c("#1f77b4", "#ff7f0e", if (param != "K") "#2ca02c")
  fill_vec <- c("#D3E4F0", "#FFE5CF", if (param != "K") "#D3E4F0")
  y_vec <- numeric()
  MSE_mat <- matrix(data = NA_real_, nrow = 0L, ncol = n_rep)
  
  for (ind in if (param == "K") 1:2 else 1:n_methods) {
    y_vec <- c(y_vec, 1:n_setups + (ind - 2L) * 0.19)
    MSE_mat <- rbind(MSE_mat, MSE_arr[, , ind])
  }
  
  # We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
  cairo_pdf(
    file = paste0("simulation_", outcome, "_", param, ".pdf"), width = 8,
    height = 9
  )
  
  par(mar = c(5, 9, 4, 2) + 0.1)
  
  boxplot(
    x = t(MSE_mat), border = rep(x = col_vec, each = n_setups),
    col = rep(x = fill_vec, each = n_setups), log = if (log_scale) "x" else "",
    boxwex = 0.12, staplewex = 1, horizontal = TRUE, at = y_vec, pch = 4L,
    cex = 0.5, pt.lwd = 0.1, lwd = lwd, lty = 1L,
    outcol = rep(
      x = adjustcolor(col = col_vec, alpha.f = 0.7), each = n_setups
    ),
    yaxt = "n", xlim = c(n_setups + 1.5, 1),
    ylim = if (outcome == "bias" && param == "alpha") c(-0.8, 0.4) else NULL,
    xlab = parse(text = paste0(outcome, "[italic(", switch(
      EXPR = param, alpha = "\u03B1", beta = "\u03B2",
      K = paste0("K)]~or~", outcome, "[italic(K[r]")
    ), ")]"))
  )
  
  abline(v = 0, lty = 3L)
  title(ylab = "Simulation scenario", mgp = c(8, 1, 0))
  abline(h = 1:(n_setups - 1L) + 0.5, lty = 3L)
  abline(h = cumsum(c(3, 3, 2, 3)) + 0.5)
  axis(1)  # Default x-axis
  
  axis(
    side = 2L, at = 1:n_setups,
    labels = parse(text = c(
      paste0(c("No", "10*'%'", "50*'%'"), "~censoring"),
      paste0("italic(n)~'='~", n_vec), paste0("italic(p)~'='~", c(3, 5)),
      paste0("(italic(n)*','~italic(m))~'='~(", c(
        paste0("500*','~", n_risk_vec), paste0("1000*','~", n_risk_vec)
      ), ")")
    )),
    las = 1L
  )
  
  legend(
    x = "bottomleft",
    legend = if (param == "K") {
      c(expression(italic(K)), expression(italic(K[r])))
    } else {
      c("MVB detector", "MLE", "King and Weiss")
    },
    pt.bg = fill_vec, col = col_vec, pt.lwd = lwd, pch = 22L, pt.cex = 2
  )
  
  box()
  dev.off()
}
