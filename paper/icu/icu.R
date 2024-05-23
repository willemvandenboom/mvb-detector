### Application to ICU length of stay

# This code produces the results of the application to ICU length of stay in the
# paper "The Multivariate Bernoulli detector: Change point estimation in
# discrete survival analysis" by Willem van den Boom, Maria De Iorio, Fang Qian
# and Alessandra Guglielmi.


# Read in the CSV file with MIMIC-IV data created by the Python script icu.py.
df <- read.csv("mimic.csv")
time_to_event <- df$X
event <- df$J
cens <- df$J == 0L
t_max <- max(time_to_event[!cens])
time_to_event[time_to_event > t_max] <- t_max
n_risk <- max(event)
n <- nrow(df)
X <- df[, 2:37]
p <- ncol(X)
for (j in 1:p) if (is.character(X[[j]])) X[[j]] <- X[[j]] == "True"
X <- as.matrix(X)
tmp <- factor(x = event, levels = 0:n_risk)
levels(tmp) <- c("censored", 1:n_risk)
data <- data.frame(time = time_to_event, event = tmp, X = X)

# Grouping of predictors because of the categorical variables insurance, marital
# status, ethnicity and admissions count coded using dummy variables:
pred_group <- c(1:20, rep(21:24, c(2L, 3L, 4L, 2L)), 25:29)


## Run the MCMC of the multivariate Bernoulli detector.

n_iter <- 2e5L  # Number of MCMC iterations
burnin <- 5e4L  # Number of burn-in iterations
set.seed(1L)
print("Running the MCMC of the multivariate Bernoulli detector...")
t_start <- Sys.time()

res <- mvb.detector::run_mvb_detector(
  data = data, n_iter = n_iter, pred_group = pred_group, burnin = burnin
)

t_end <- Sys.time()
print(t_end - t_start)


## Create trace plots.

# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "mimic_trace.pdf", width = 9, height = 8)
par(mfrow = 2:1)
mvb.detector::plot_trace(res)
dev.off()


## Plot posterior inference on the hazards.

xlab <- expression(paste("ICU length of stay ", italic(t), " (days)"))
risk_names <- c("Home", "Transfer", "Death")
pdf(file = "mimic_hazard.pdf", width = 9, height = 5)
par(mfrow = c(2L, 1L + n_risk), mai = c(0.5, 0.6, 0.3, 0.2))

mvb.detector::plot_change_points(
  mvbd_fit = res, xlab = xlab, BF = FALSE, main_vec = c("Overall", risk_names)
)

mvb.detector::plot_CHF(mvbd_fit = res, xlab = xlab)

mvb.detector::plot_baseline_hazards(
  mvbd_fit = res, xlab = xlab, main_vec = character(n_risk)
)

dev.off()


# Plot Bayes factors as an alternative to posterior probabilities of change
# points.
pdf(file = "mimic_hazard_BF.pdf", width = 9, height = 2.5)
par(mfrow = c(1L, 1L + n_risk), mai = c(0.5, 0.6, 0.3, 0.2))

mvb.detector::plot_change_points(
  mvbd_fit = res, xlab = xlab, BF = TRUE, main_vec = c("Overall", risk_names)
)

dev.off()


## Plot posterior inference on the regression coefficients.

# Pretty predictor names for plotting
pred_names <- c(
  "Anion gap", "Bicarbonate", "Calcium total", "Chloride", "Creatinine",
  "Glucose", "Magnesium", "Phosphate", "Potassium", "Sodium", "Urea nitrogen",
  "Hematocrit", "Hemoglobin", "MCH", "MCH concentration", "MCV",
  "Platelet count", "RDW", "Red blood cells", "White blood cells",
  "Insurance: Medicare", "Insurance: other", "Marital status: married",
  "Marital status: single", "Marital status: widowed", "Ethnicity: black",
  "Ethnicity: Hispanic", "Ethnicity: other", "Ethnicity: white",
  "Admission number: 2", "Admission number: 3+", "Night admission",
  "Sex: female", "Direct emergency", "Recent admission", "Age"
)

# Names for each group of predictors:
cov_names <- pred_names[match(1:29, pred_group)]
cov_names[21] <- "Insurance"
cov_names[22] <- "Marital status"
cov_names[23] <- "Ethnicity"
cov_names[24] <- "Admission number"
cov_names[26] <- "Sex"

# Plot posterior inclusion probabilities.
# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "mimic_inc_prob.pdf", width = 6.5, height = 14)
par(mar = c(4, 9.3, 0, 0) + 0.1)  # Set margins to make room for x-axis labels.

mvb.detector::plot_regression(
  mvbd_fit = res, risk_names = risk_names, pred_names = cov_names,
  inc_prob = TRUE, pred_group = pred_group
)

dev.off()


plot_beta <- function (beta_CI, file, pred_names_offset = 0.07) {
  # We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
  cairo_pdf(file = file, width = 6.5, height = 14)
  
  # Set margins to make room for x-axis labels.
  par(mar = c(4, 9.3, 0, 0) + 0.1)
  
  mvb.detector::plot_regression(
    beta_CI, risk_names, pred_names, pred_names_offset
  )
  
  dev.off()
}


plot_beta(res, "mimic_regression.pdf")  # Plot posterior means and 95% CIs.


## Comparison with maximum likelihood estimation

# Transform the data into long format.
data_long <- discSurv::dataLongCompRisks(
  dataShort = data, timeColumn = "time", eventColumn = "event",
  eventColumnsAsFactor = TRUE
)

X_long <- as.matrix(data_long[, 5L + n_risk + seq_len(p)])
data_long$y <- 0L
tmp <- which(data_long$e0 == 0)
data_long$y[tmp] <- as.integer(data_long$event[tmp]) - 1L

nnet_fit <- nnet::multinom(
  formula = y ~ as.factor(timeInt) + X_long, data = data_long, maxit = 1e3L
)

se_cov <- vcov(nnet_fit)
alpha_hat <- coef(nnet_fit)[, 1:t_max]
alpha_hat[, 2:t_max] <- alpha_hat[, 1] + alpha_hat[, 2:t_max]
beta_hat <- t(coef(nnet_fit)[, t_max + 1:p])
alpha_se <- matrix(data = NA_real_, nrow = n_risk, ncol = t_max)
beta_CI <- array(data = NA_real_, dim = c(p, n_risk, 3L))

for (r in 1:n_risk) {
  tmp <- (r - 1L) * (t_max + p)
  alpha_se[r, ] <- diag(se_cov)[tmp + 1:t_max]
  
  alpha_se[r, -1] <-
    alpha_se[r, -1] + alpha_se[r, 1] + 2 * se_cov[tmp + 1L, tmp + 2:t_max]
  
  for (j in 1:p) beta_CI[j, r, ] <-
    beta_hat[j, r] + c(0, -1, 1) * 1.96 * sqrt(diag(se_cov)[tmp + j])
}

alpha_se <- sqrt(alpha_se)
alpha_CI <- array(data = NA_real_, dim = c(2L, n_risk, t_max))

for (r in 1:n_risk) for (t in 1:t_max) {
  alpha_CI[, r, t] <- alpha_hat[r, t] + c(-1, 1) * 1.96 * alpha_se[r, t]
}

# Plot baseline hazards.
# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "mimic_nnet_hazard.pdf", width = 8, height = 3.5)

# Set margins to make room for y-axis labels.
par(mfrow = c(1L, n_risk), mar = c(5, 5, 4, 2) + 0.1)

for (r in 1:n_risk) mvb.detector:::plot_CI(
  vec = alpha_hat[r, ], CI = alpha_CI[, r, ], ylim = range(alpha_CI),
  xlab = xlab, ylab = parse(text = paste("italic(\u03B1[", r, "][t])", sep="")),
  main = risk_names[r]
)

dev.off()

plot_beta(beta_CI, "mimic_nnet_regression.pdf", 0.1)  # Plot coefficients


## Comparison with the model by King & Weiss
## (2021, doi:10.1093/biostatistics/kxz029)


# Modify `brea::brea_mcmc` to show a progress bar.
brea_mcmc2 <- function (x, y, priors = NULL, S = 1000, B = 100, n = NULL, K = NULL, 
    store_re = FALSE) 
{
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


age_ordered <- sort(unique(data_long$X.standardized_age))
data_long$age_discretized <- match(data_long$X.standardized_age, age_ordered)
data_long$age_discretized[data_long$age_discretized == 73L] <- 74L
set.seed(1L)
print("Running the MCMC of the model by King & Weiss (2021)...")

brea_fit <- brea_mcmc2(
  x = cbind(data_long$timeInt, data_long$age_discretized, X_long[, -p] + 1L),
  y = as.matrix(data_long[, 3L + 1:n_risk]),
  priors = c(
    rep(list(list("gmrf", 3, .01)), 2L), rep(list(list("cat", 4)), p - 1L)
  ),
  S = n_iter, B = burnin,
  K = c(t_max, length(age_ordered) + 1L, rep(2L, p - 1L))
)

# Create trace plots.
# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "brea_trace.pdf", width = 9, height = 8)

# Set margins to make room for y-axis labels.
par(mfrow = c(n_risk, 1L), mar = c(5, 5, 4, 2) + 0.1)

# Trace plots of alpha
for (r in 1:n_risk) plot(
  x = brea_fit$b_0_s[r, -(1:burnin)], type = "l", xlab = "MCMC iteration",
  ylab = parse(text = paste("italic(\u03B2[0*", r, "])", sep="")),
  main = risk_names[r]
)

dev.off()

# Inference on the baseline hazards
alpha_brea_hat <- matrix(data = NA_real_, nrow = n_risk, ncol = t_max)
alpha_brea_CI <- array(data = NA_real_, dim = c(2L, n_risk, t_max))
offset <- rowMeans(brea_fit$b_0_s)

for (r in 1:n_risk) {
  alpha_brea_hat[r, ] <- offset[r] + rowMeans(brea_fit$b_m_s[[1]][r, , ])
}

for (r in 1:n_risk) for (i in 1:2) alpha_brea_CI[i, r, ] <- offset[r] + apply(
  X = brea_fit$b_m_s[[1]][r, , ], MARGIN = 1L, FUN = quantile,
  probs = c(.025, .975)[i]
)


## Plot baseline hazards.
# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "mimic_brea_hazard.pdf", width = 8, height = 3.5)

# Set margins to make room for y-axis labels.
par(mfrow = c(1L, n_risk), mar = c(5, 5, 4, 2) + 0.1)

for (r in 1:n_risk) mvb.detector:::plot_CI(
  vec = alpha_brea_hat[r, ], CI = alpha_brea_CI[, r, ],
  ylim = range(alpha_brea_CI), xlab = xlab,
  ylab = parse(text = paste("italic(\u03B1[", r, "][t])", sep = "")),
  main = risk_names[r]
)

dev.off()


## Plot inference on regression terms
p_binary <- p - 1L
beta_brea_CI <- array(data = NA_real_, dim = c(p_binary, n_risk, 3L))

for (r in 1:n_risk) for (j in 1:p_binary) {
  tmp <- brea_fit$b_m_s[[1L + p - p_binary + j]][r, , ]
  tmp <- tmp[2, ] - tmp[1, ]
  beta_brea_CI[j, r, 1] <- mean(tmp)
  beta_brea_CI[j, r, 2:3] <- quantile(x = tmp, probs = c(0.025, 0.975))
}

plot_beta(beta_brea_CI, "mimic_brea_regression.pdf")

# We use `cairo_pdf` instead of `pdf` to correctly display Unicode characters.
cairo_pdf(file = "mimic_brea_regression_age.pdf", width = 8, height = 3.5)

# Set margins to make room for y-axis labels.
par(mfrow = c(1L, n_risk), mar = c(5, 5, 4, 2) + 0.1)

for (r in 1:n_risk) mvb.detector:::plot_CI(
  vec = rowMeans(brea_fit$b_m_s[[2]][r, , ]), CI = apply(
    X = brea_fit$b_m_s[[2]][r, , ], MARGIN = 1L, FUN = quantile,
    probs = c(0.025, 0.975)
  ),
  ylim = c(-1.7, 1.7),
  xlab = expression(paste("Age (relative, starts at 1) ", italic(x[1]))),
  ylab = parse(text = paste("italic(f[\u03B2][", r, "1])(x[1])", sep="")),
  main = risk_names[r]
)

dev.off()
