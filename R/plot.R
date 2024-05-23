### Code for plotting
K_lab <- expression(paste("Number of change points ", italic(K)))


get_risk_names <- function (mvbd_fit) {
  risk_names <- levels(mvbd_fit$data$event)[-1]
  n_risk <- length(risk_names)
  
  return (if (all(risk_names == 1:n_risk)) {
    if (n_risk == 1L) "" else paste("Risk", 1:n_risk)
  } else {
    risk_names
  })
}


plot_CI <- function (vec, CI, ylim, xlab, ylab, vec_true = NULL, main = "") {
  t_max <- length(vec)
  
  plot(
    x = vec, type = "n", lwd = 2, ylim = ylim, ylab = ylab, xlab = xlab,
    xaxs = "i", yaxs = "i", main = main
  )
  
  graphics::polygon(
    x = c(1:t_max, t_max:1), y = c(CI[1, ], rev(CI[2, ])), border = NA,
    col = "darkgrey"
  )
  
  graphics::lines(vec, lwd = 2)
  
  if (!is.null(vec_true)) {
    graphics::lines(vec_true, col = "#d62728", lty = 2, lwd = 2)
  }
}


#' Plot baseline hazards
#' 
#' Plot the posterior inference on the baseline hazards
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param alpha_true Matrix with true population parameters
#' @param xlab x-axis label
#' @param ylim Numeric vector of length 2, giving the y coordinates range
#' @param main_vec Plot titles
#'  
#' @export
plot_baseline_hazards <- function (
  mvbd_fit, alpha_true = NULL, xlab = expression(paste("Time ", italic(t))),
  ylim = NULL, main_vec = NULL
) {
  alpha_MCMC <- mvbd_fit$alpha_MCMC
  alpha_mean <- colMeans(alpha_MCMC)
  n_risk <- nrow(alpha_mean)
  t_max <- ncol(alpha_mean)
  if (is.null(main_vec)) main_vec <- get_risk_names(mvbd_fit)
  alpha_CI <- array(NA_real_, dim = c(2, n_risk, t_max))
  
  for (r in 1:n_risk) for (t in 1:t_max) {
    alpha_CI[, r, t] <- stats::quantile(alpha_MCMC[, r, t], c(0.025, 0.975))
  }
  
  for (r in 1:n_risk) plot_CI(
    vec = alpha_mean[r, ], CI = alpha_CI[, r, ],
    ylim = if (is.null(ylim)) range(alpha_CI) else ylim, xlab = xlab,
    ylab = parse(text = paste("italic(alpha[", r, "][t])", sep="")),
    vec_true = if (!is.null(alpha_true)) alpha_true[r, ] else NULL,
    main = main_vec[r]
  )
}


#' Plot cumulative hazard function
#' 
#' Plot the posterior inference on the cumulative hazard function for a subject
#' with all covariates equal to zero.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param xlab x-axis label
#'  
#' @export
plot_CHF <- function (
  mvbd_fit, xlab = expression(paste("Time ", italic(t)))
) {
  alpha_MCMC <- mvbd_fit$alpha_MCMC
  n_iter_post <- dim(alpha_MCMC)[1]
  t_max <- dim(alpha_MCMC)[3]
  CHF_post <- array(data = NA_real_, dim = c(n_iter_post, t_max))
  
  for (s in 1:n_iter_post) {
    CHF_post[s, 1] <- sum(exp(log_hazard_rate_compete(alpha_MCMC[s, , 1])))
    
    for (t in 2:t_max) {
      CHF_post[s, t] <- CHF_post[s, t - 1L] + (1 - CHF_post[s, t - 1L]) * sum(
        exp(log_hazard_rate_compete(alpha_MCMC[s, , t]))
      )
    }
  }
  
  CHF_CI <- array(NA_real_, dim = c(2, t_max))
  
  for (t in 1:t_max) {
    CHF_CI[, t] <- stats::quantile(CHF_post[, t], c(0.025, 0.975))
  }
  
  plot_CI(
    vec = colMeans(CHF_post), CI = CHF_CI, ylim = 0:1, xlab = xlab,
    ylab = "Cumulative hazard function"
  )
}


plot_cp_single <- function (
  BF, prob_vec, BF_vec, t_max, xlab, cp_true, main, BF_max, prob_lab, BF_lab
) {
  # Plot posterior probability and Bayes factor for the overall change points or
  # a single risk only.
  # The code is based on https://stackoverflow.com/a/52702095.
  both <- is.null(BF)
  
  plot_vec <- function (prob) {
    x <- 1:t_max + if (!both) 0 else (1L - 2L * prob) / 5
    vec <- if (prob) prob_vec else BF_vec
    col <- if (!both) "black" else if (prob) "#1f77b4" else "#ff7f0e"
    side <- if (prob || !both) 2L else 4L
    
    graphics::plot.window(
      xlim = c(1, t_max), ylim = c(0, if (prob) 1 else BF_max)
    )
    
    graphics::segments(
      x0 = x, y0 = 0, y1 = vec, col = if (both) col else "darkgrey",
      lwd = if (both) 1.5 else 3, lend = 1L
    )
    
    if(!is.null(cp_true)) graphics::abline(
      v = cp_true + 1L, lty = "dashed", col = "#d62728",
      lwd = if (both) 1 else 2
    )
    
    graphics::points(
      x = x, y = vec, pch = if (prob) 4L else 20L, col = col, cex = 0.75
    )
    
    graphics::axis(side = side, col.axis = col)
    
    graphics::mtext(
      text = if (prob) prob_lab else BF_lab, side = side, las = 3L, line = 2,
      col = col, cex = 0.7
    )
  }
  
  graphics::plot.new()
  graphics::abline(h = 0, lwd = 0.5)
  if (!isTRUE(BF)) plot_vec(TRUE)
  if (!isFALSE(BF)) plot_vec(FALSE)
  graphics::axis(1L)
  graphics::box()
  graphics::title(xlab = xlab, main = main)
}


logSumExp <- function (x) {
  x_max <- max(x)
  return (x_max + log(sum(exp(x - x_max))))
}


#' Plot posterior inference on change points
#' 
#' Plot the posterior inference on the location of change points using Bayes
#' factors or posterior probabilities.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param n_vec_true Vector of true population regime durations
#' @param z_star_true Array specifying which risk are affected at change points
#' @param BF Whether to plot Bayes factors (TRUE), posterior probabilities
#'  (FALSE) or both (NULL, default)
#' @param xlab x-axis label
#' @param main_vec Plot titles
#' @param BF_max Upper limit of the y-axis range for the Bayes factors
#'  
#' @export
plot_change_points <- function (
  mvbd_fit, n_vec_true = NULL, z_star_true = NULL, BF = NULL,
  xlab = expression(paste("Time ", italic(t))), main_vec = NULL, BF_max = 100
) {
  n_risk <- nlevels(mvbd_fit$data$event) - 1L
  n_iter_post <- dim(mvbd_fit$alpha_MCMC)[1]
  t_max <- dim(mvbd_fit$alpha_MCMC)[3]
  
  if(n_risk > 1L & !is.null(n_vec_true) & is.null(z_star_true)) stop(
    "If `n_vec_true` is specified, then `z_star_true` should be specified too."
  )
  
  if (is.null(main_vec)) {
    main_vec <- if (n_risk == 1L) "" else c("Overall", get_risk_names(mvbd_fit))
  }
  
  gam_mean <- numeric(t_max)
  z_mean <- matrix(data = 1, nrow = n_risk, ncol = t_max)
  gam_mean[1] <- 1
  z_t_MCMC <- matrix(data = NA, nrow = n_risk, ncol = n_iter_post)
  
  for (t in 2:t_max) {
    for (r in 1:n_risk) for (s in 1:n_iter_post) z_t_MCMC[r, s] <-
      mvbd_fit$alpha_MCMC[s, r, t] != mvbd_fit$alpha_MCMC[s, r, t - 1L]
    
    z_mean[, t] <- rowMeans(z_t_MCMC)
    gam_mean[t] <- mean(colSums(z_t_MCMC) > 0)
  }
  
  cp_candidates <- get_cp_candidates(
    mvbd_fit$data, t_max, mvbd_fit$cp_constraint
  )
  
  gam_mean[-cp_candidates] <- NA_real_
  z_mean[, -cp_candidates] <- NA_real_
  gam_mean[1] <- NA_real_
  z_mean[, 1] <- NA_real_
  
  # Compute prior marginal probability of a change point
  n_cp_candidates <- length(cp_candidates)
  
  log_prior_n_regimes_normalizing_constant <-
    logSumExp(log_prior_n_regimes(1:n_cp_candidates))
  
  log_prior_expec_n_regimes <- sum((1:n_cp_candidates) * exp(
    log_prior_n_regimes(1:n_cp_candidates) -
      log_prior_n_regimes_normalizing_constant
  ))
  
  prior_prob_cp <- (log_prior_expec_n_regimes - 1) / (n_cp_candidates - 1)
  prior_prob_cp_single_risk <- prior_prob_cp * 2^(n_risk - 1L) / (2^n_risk - 1L)
  
  Bayes_factor <-
    (gam_mean / (1 - gam_mean)) / (prior_prob_cp / (1 - prior_prob_cp))
  
  Bayes_factor[is.infinite(Bayes_factor)] <- 2 * BF_max
  
  plot_cp_single(
    BF = BF, prob_vec = gam_mean, BF_vec = Bayes_factor, t_max = t_max,
    xlab = xlab,
    cp_true = if(is.null(n_vec_true)) NULL else {
      cumsum(n_vec_true[-length(n_vec_true)])
    },
    main = main_vec[1], BF_max = BF_max,
    prob_lab = parse(
      text = "italic(p)(italic(gamma)[italic(t)] != 0~'|'~data)"
    ),
    BF_lab = expression(paste("Bayes factor of ", italic(gamma[t])))
  )
  
  if (n_risk > 1L) for (r in 1:n_risk) {
    Bayes_factor <- (z_mean[r, ] / (1 - z_mean[r, ])) / (
      prior_prob_cp_single_risk / (1 - prior_prob_cp_single_risk)
    )
    
    Bayes_factor[is.infinite(Bayes_factor)] <- 2 * BF_max
    Bayes_factor[-cp_candidates] <- NA_real_
    Bayes_factor[1] <- NA_real_
    
    plot_cp_single(
      BF = BF, prob_vec = z_mean[r, ], BF_vec = Bayes_factor, t_max = t_max,
      xlab = xlab,
      cp_true = if(is.null(n_vec_true)) NULL else {
        cumsum(n_vec_true[-length(n_vec_true)])[z_star_true[r, -1]]
      },
      main = main_vec[r + 1L], BF_max = BF_max,
      prob_lab = parse(text = paste(
        "italic(p)(italic(z)[", r, "][italic(t)] != 0~'|'~data)",
        sep=""
      )),
      BF_lab = parse(
        text = paste("Bayes~factor~of~italic(z[", r, "][t])", sep="")
      )
    )
  }
}


#' Plot MCMC trace plots
#' 
#' Plot trace plots of the log-likelihood and/or the number of change points.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param quantity Vector specifying which trace plots are generated
#'  
#' @export
plot_trace <- function (mvbd_fit, quantity = c("log_lik", "K")) {
  for (y in quantity) graphics::plot(
    x = mvbd_fit[[paste0(
      if (y == "log_lik") "log_lik" else "n_regimes", "_MCMC"
    )]] - (y == "K"),
    type = "l", xlab = "MCMC iteration",
    ylab = if (y == "log_lik") "Log-likelihood" else K_lab
  )
}


get_beta_CI <- function (mvbd_fit) {
  # Compute 95% confidence intervals of the regression coefficients.
  beta_MCMC <- mvbd_fit$beta_MCMC
  p <- dim(beta_MCMC)[2]
  n_risk <- dim(beta_MCMC)[3]
  beta_CI <- array(data = NA_real_, dim = c(p, n_risk, 3L))
  beta_CI[, , 1] <- colMeans(beta_MCMC)
  
  for (r in 1:n_risk) for (j in 1:p) beta_CI[j, r, 2:3] <- stats::quantile(
    x = beta_MCMC[, j, r], prob = c(0.025, 0.975)
  )
  
  return (beta_CI)
}


#' Compute posterior inclusion probabilities
#' 
#' Compute the posterior inclusion probabilities of the regression coefficients.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param risk_names Names of the competing risks
#' @param pred_names Names of the predictors
#' @param digits If provided, the result is rounded to this many digits.
#' 
#' @returns Matrix of posterior inclusion probabilities with each row
#'  corresponding to a predictor and each column to a risk. The row names of the
#'  matrix are `pred_names`, if provided. The column names or `risk_names` or,
#'  if multiple risks, "Risk 1", etc.
#'  
#' @export
compute_post_inc_prob <- function (
    mvbd_fit, risk_names = NULL, pred_names = NULL, digits = NULL
) {
  res <- colMeans(mvbd_fit$beta_MCMC != 0)
  if (!missing(pred_names)) rownames(res) <- pred_names
  
  colnames(res) <- if (missing(risk_names)) {
    get_risk_names(mvbd_fit)
  } else {
    risk_names
  }
  
  return (if (is.null(digits)) res else round(res, digits))
}


plot_single_CI <- function(x, CI, lty, x_offset) {
  # Draw a single credible interval.
  
  # Draw the credible interval.
  graphics::lines(y = rep(x, 2), x = CI[-1], lty = lty, lwd = 2)
  
  # Add horizontal delimiters to both ends of the credible interval.
  for (y in CI[-1]) graphics::lines(
    y = c(x - x_offset / 3, x + x_offset / 3), x = rep(y, 2L), lwd = 2
  )
  
  # Add a dot at the posterior mean.
  graphics::points(y = x, x = CI[1], pch = 19L, cex = 0.7)
}


#' Plot posterior inference on regression coefficients
#' 
#' Plot posterior means and 95\% credible intervals of the regression
#' coefficients.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#' @param risk_names Names of the competing risks for use in the legend
#' @param pred_names Names of the predictors shown next to the y-axis
#' @param pred_names_offset Vertical offset from the y-axis of the labels in
#'  `pred_names`
#' @param inc_prob If true, the posterior inclusion probabilities are plot
#'  instead of the posterior means and credible intervals.
#' @param pred_group Group indicator for predictors that should be jointly added
#'  or removed when it comes to variable selection. Only used if
#'  inc_prob = TRUE.
#'  
#' @export
plot_regression <- function (
  mvbd_fit, risk_names = NULL, pred_names = NULL, pred_names_offset = 0.07,
  inc_prob = FALSE, pred_group = NULL
) {
  beta_CI <- if (is.array(mvbd_fit)) mvbd_fit else get_beta_CI(mvbd_fit)
  p <- dim(beta_CI)[1]
  n_risk <- dim(beta_CI)[2]
  xlim = if (inc_prob) 0:1 else range(beta_CI)
  x_offset <- .2
  
  if (inc_prob) {
    if (is.null(pred_group)) pred_group <- 1:p
    groups <- unique(pred_group)
    p <- length(groups)
    inc_prob_mat <- compute_post_inc_prob(mvbd_fit)
    
    # Only keep one row of `inc_prob_mat` per covariate.
    inc_prob_mat <- inc_prob_mat[match(1:p, pred_group), ]
  }
  
  # Set up empty plot.
  plot(
    x = 0, type = "n", xlim = xlim, ylim = c(-1.9, p + 2 * x_offset), ylab = "",
    xlab = if (inc_prob) {
      "Posterior inclusion probability"
    } else {
      parse(text = "'Regression coefficient'~italic(\u03B2[rj])")
    },
    yaxt = if (missing(pred_names)) "s" else "n", yaxs = "i"
  )
  
  # Add x-axis tick mark text.
  if (!missing(pred_names)) graphics::text(
    y = p:1, x = xlim[1] - pred_names_offset, adj = 1, labels = pred_names,
    xpd = TRUE
  )
  
  # Add a line at one or zero.
  graphics::abline(v =  c(0, if (inc_prob) 1), col = "grey")
  
  # Plot the p credible intervals or inclusion probabilities.
  for (r in 1:n_risk) for (j in 1:p) plot_single_CI(
    x = p + 1 - j - (
      r - 0.5 * (n_risk + 1L)
    ) / (0.5 * (n_risk - 1L)) * x_offset,
    CI = if (inc_prob) c(NA_real_, 0, inc_prob_mat[j, r]) else beta_CI[j, r, ],
    lty = r, x_offset = x_offset
  )
  
  if (n_risk > 1L) {
    # Add horizontal lines for legibility.
    for (j in seq_len(p - 1L)) {
      graphics::abline(h = j + 0.5, lwd = 0.5, col = "grey")
    }
  
    graphics::legend(
      x = "bottomleft",
      legend = if (missing(risk_names)) {
        get_risk_names(mvbd_fit)
      } else {
        risk_names
      },
      lty = 1:n_risk, lwd = 2
    )
  }
}


#' Plot posterior inference on number of change points
#' 
#' Plot posterior distribution of the number of overall change points.
#' 
#' @param mvbd_fit Output from the function \link{run_mvb_detector}
#'  
#' @export
plot_K <- function (mvbd_fit) {
  K_MCMC <- mvbd_fit$n_regimes_MCMC - 1L
  
  graphics::barplot(
    height = table(
      factor(x = K_MCMC, levels = min(K_MCMC):max(K_MCMC))
    ) / length(K_MCMC),
    xlab = K_lab, ylab = "Posterior probability"
  )
}
