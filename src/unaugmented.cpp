// Code for the unaugmented updates, i.e. global MCMC

// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>  // For `std::min`
#include <bitset>  // For `std::bitset`
#include <cmath>  // For `std::log1p`, `std::log2`, `std::pow`, etc.
#include <string>  // For `std::string`
#include <vector>  // For `std::vector`

#include <RcppArmadillo.h>


//' Log of the prior on alpha conditionally on the change points
//' 
//' @param alpha_star Matrix of unique values of alpha over each time interval
//' @param mu_alpha Mean of the marginal prior on alpha
//' @param s2_alpha Variance of the marginal prior on alpha
//' 
//' @keywords internal
//' @name log_prior_alpha
double log_prior_alpha(
  const Rcpp::NumericMatrix alpha_star, const double mu_alpha,
  const double s2_alpha
) {
  int n_risk = alpha_star.nrow();
  double ret = 0.0;
  
  for (int r = 0; r < n_risk; r++) ret += Rcpp::sum(Rcpp::dnorm(
    Rcpp::unique(alpha_star(r, Rcpp::_)), mu_alpha, std::sqrt(s2_alpha), true
  ));
  
  return ret;
}


//' Log of the hazard rate
//'
//' @param alpha Vector of baseline hazard parameters
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector log_hazard_rate_compete(const Rcpp::NumericVector alpha) {
  // This function assumes that `alpha` is an `n_risk`-dimensional vector.
  return alpha - std::log1p(Rcpp::sum(Rcpp::exp(alpha)));
}


// [[Rcpp::export]]
double log_likelihood(
  const Rcpp::IntegerVector n_vec, const Rcpp::NumericMatrix alpha_star,
  const Rcpp::IntegerVector event, const Rcpp::IntegerVector time_to_event,
  const Rcpp::LogicalVector cens, const Rcpp::NumericMatrix mu
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


//' Log of posterior for discrete-time survival with `n_risk` risks
//' 
//' Here, the prior on cause-specific change points p(z_t | gamma_t) is not
//' included in `log_prior_n_vec`, and thus not included in this log-posterior
//' computation.
//' 
//' @keywords internal
// [[Rcpp::export]]
double log_posterior_unaugmented(
  const Rcpp::IntegerVector n_vec, const Rcpp::NumericMatrix alpha_star,
  const arma::mat beta, const arma::mat X,
  const Rcpp::IntegerVector time_to_event, const Rcpp::IntegerVector event,
  const Rcpp::LogicalVector cens, const Rcpp::Function log_prior_n_vec,
  const double mu_alpha, const double s2_alpha
) {

  return Rcpp::as<double>(log_prior_n_vec(n_vec)) +
    log_prior_alpha(alpha_star, mu_alpha, s2_alpha) + log_likelihood(
        n_vec, alpha_star, event, time_to_event, cens,
        Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(X * beta))
    );
}


double runif() {
  return R::runif(0.0, 1.0);
}


int sample(const Rcpp::IntegerVector vec) {
  if (vec.length() == 1) return vec[0];
  Rcpp::Function sample_R("sample");
  PutRNGstate();
  
  int res = Rcpp::as<int>(
    sample_R(Rcpp::Named("x") = vec, Rcpp::Named("size") = 1)
  );
  
  GetRNGstate();
  return res;
}


// [[Rcpp::export]]
int sample_int(const int n) {
  Rcpp::Function sample_int_R("sample.int");
  PutRNGstate();
  
  int res = Rcpp::as<int>(
    sample_int_R(Rcpp::Named("n") = n, Rcpp::Named("size") = 1)
  );
  
  GetRNGstate();
  return res;
}


// [[Rcpp::export]]
Rcpp::LogicalVector sample_z(const int n_risk) {
  int n_z = int(std::pow(2, n_risk)) - 1;
  
  // Sample from the uniform distribution over the `n_z` possible `z`.
  int k = sample_int(n_z);
  
  // Transform the integer `k` to a corresponding Boolean vector `z`.
  // This code supports at most 64 competing risks.
  std::string bit_str = std::bitset<64>(k).to_string();
  
  std::vector<bool> bit_vec;
  bit_vec.reserve(n_risk);
  
  for (auto bit : bit_str.substr(bit_str.size() - n_risk)) {
    bit_vec.push_back(bit == '1');
  }
  
  return Rcpp::as<Rcpp::LogicalVector>(Rcpp::wrap(bit_vec));
}


// Standard deviation of the Gaussian random walk for new alpha
double s_prop_alpha = 1.0;


// [[Rcpp::export]]
Rcpp::List update_n_vec_global(
    Rcpp::IntegerVector n_vec, Rcpp::NumericMatrix alpha_star,
    Rcpp::LogicalMatrix z_star, const arma::mat beta,
    const Rcpp::IntegerVector time_to_event, const Rcpp::IntegerVector event,
    const Rcpp::LogicalVector cens, const arma::mat X, const int n_risk,
    const int t_max, const Rcpp::Function log_prior_n_vec,
    const double mu_alpha, const double s2_alpha
) {
  double q = 0.5;  // Split proposal probability
  
  // No. of regimes is no. of change points + 1.
  int n_regimes = n_vec.length();

  double log_post = log_posterior_unaugmented(
      n_vec, alpha_star, beta, X, time_to_event, event, cens, log_prior_n_vec,
      mu_alpha, s2_alpha
  );

  if (n_regimes == 1 || (n_regimes < t_max && runif() < q)) {  // Split
    Rcpp::IntegerVector splittable = Rcpp::seq(0, n_regimes - 1);
    splittable = splittable[n_vec > 1];
    int j = sample(splittable);
    int l = sample_int(n_vec[j] - 1);
    Rcpp::IntegerVector n_vec_prop = Rcpp::clone(n_vec);
    n_vec_prop.erase(j);
    n_vec_prop.insert(j, l);
    n_vec_prop.insert(j + 1, n_vec[j] - l);
    Rcpp::LogicalVector z_new = sample_z(n_risk);
    Rcpp::NumericVector alpha_new = alpha_star(Rcpp::_, j);
    
    Rcpp::NumericVector tmp_vec =
      Rcpp::as<Rcpp::NumericVector>(alpha_new[z_new]) +
      s_prop_alpha * Rcpp::rnorm(Rcpp::sum(z_new));
    
    alpha_new[z_new] = tmp_vec;
    
    Rcpp::NumericMatrix tmp_mat_N = alpha_star(Rcpp::_, Rcpp::seq(0, j)),
      alpha_star_prop = Rcpp::cbind(tmp_mat_N, alpha_new);
    
    if (n_regimes > j + 1) {
      tmp_mat_N = alpha_star(Rcpp::_, Rcpp::seq(j + 1, n_regimes - 1));
      alpha_star_prop = Rcpp::cbind(alpha_star_prop, tmp_mat_N);
    }
    
    for (int r = 0; r < n_risk; r++) for (int t = j + 1; t < n_regimes; t++) {
      if (z_star(r, t)) break;
      alpha_star_prop(r, t + 1) = alpha_new[r];
    }
    
    double log_post_prop = log_posterior_unaugmented(
      n_vec_prop, alpha_star_prop, beta, X, time_to_event, event, cens,
      log_prior_n_vec, mu_alpha, s2_alpha
    );
    
    tmp_vec = alpha_star(Rcpp::_, j);
    
    if (runif() < std::exp(
      (n_regimes < t_max - 1) * std::log(1 - q) -
        (n_regimes > 1) * std::log(q) + std::log(splittable.length()) +
        std::log(n_vec[j] - 1) - std::log(n_regimes) + log_post_prop -
        log_post - Rcpp::sum(Rcpp::dnorm(
            Rcpp::as<Rcpp::NumericVector>(alpha_new[z_new] - tmp_vec[z_new]),
            0.0, s_prop_alpha, true
        ))
    )) {
      log_post = log_post_prop;
      n_vec = n_vec_prop;
      alpha_star = alpha_star_prop;
      
      Rcpp::LogicalMatrix tmp_mat_L1 = z_star(Rcpp::_, Rcpp::seq(0, j));
      
      if (n_regimes > j + 1) {
        Rcpp::LogicalMatrix tmp_mat_L2 =
          z_star(Rcpp::_, Rcpp::seq(j + 1, n_regimes - 1));
        
        z_star = Rcpp::cbind(tmp_mat_L1, z_new, tmp_mat_L2);
      } else {
        z_star = Rcpp::cbind(tmp_mat_L1, z_new);
      }
    }
  } else {  // Merge
    int j = sample_int(n_regimes - 1) - 1;
    int n_merge = n_vec[j] + n_vec[j + 1];
    Rcpp::IntegerVector n_vec_prop = Rcpp::clone(n_vec);
    n_vec_prop.erase(j);
    n_vec_prop.erase(j);
    n_vec_prop.insert(j, n_merge);
    Rcpp::NumericMatrix alpha_star_prop = alpha_star(Rcpp::_, Rcpp::seq(0, j));
    
    if (n_regimes > j + 2) {
      Rcpp::NumericMatrix tmp_mat_N =
        alpha_star(Rcpp::_, Rcpp::seq(j + 2, n_regimes - 1));
      
      alpha_star_prop = Rcpp::cbind(alpha_star_prop, tmp_mat_N);
    }
    
    for (int r = 0; r < n_risk; r++) for (int t = j + 2; t < n_regimes; t++) {
      if (z_star(r, t)) break;
      alpha_star_prop(r, t - 1) = alpha_star(r, j);
    }
    
    double log_post_prop = log_posterior_unaugmented(
      n_vec_prop, alpha_star_prop, beta, X, time_to_event, event, cens,
      log_prior_n_vec, mu_alpha, s2_alpha
    );
    
    Rcpp::NumericVector tmp_vec1 = alpha_star(Rcpp::_, j + 1),
      tmp_vec2 = alpha_star_prop(Rcpp::_, j);

    if (runif() < std::exp(
      (n_regimes > 2L) * std::log(q) - (n_regimes < t_max) * std::log(1 - q) +
        std::log(n_regimes - 1) - std::log((double) Rcpp::sum(n_vec_prop > 1)) -
        std::log(n_merge - 1) + log_post_prop - log_post +
        Rcpp::sum(Rcpp::dnorm(
          Rcpp::as<Rcpp::NumericVector>(
            tmp_vec1[z_star(Rcpp::_, j + 1)] - tmp_vec2[z_star(Rcpp::_, j + 1)]
          ), 0.0, s_prop_alpha, true
        ))
    )) {
      log_post = log_post_prop;
      n_vec = n_vec_prop;
      alpha_star = alpha_star_prop;
      Rcpp::LogicalMatrix tmp_mat_L1 = z_star(Rcpp::_, Rcpp::seq(0, j));
      
      if (n_regimes > j + 2) {
        Rcpp::LogicalMatrix tmp_mat_L2 =
          z_star(Rcpp::_, Rcpp::seq(j + 2, n_regimes - 1));
        
        z_star = Rcpp::cbind(tmp_mat_L1, tmp_mat_L2);
      } else {
        z_star = tmp_mat_L1;
      }
    }
  }
  
  
  n_regimes = n_vec.length();  // No. of regimes is no. of change points + 1.

  if (n_regimes > 1) {  // Shuffle
    int i = sample_int(n_regimes - 1) - 1;
    int j = sample_int(n_vec[i] + n_vec[i + 1] - 1);
    Rcpp::IntegerVector n_vec_prop = Rcpp::clone(n_vec);
    n_vec_prop[i + 1] = n_vec[i] + n_vec[i + 1] - j;
    n_vec_prop[i] = j;

    double log_post_prop = log_posterior_unaugmented(
      n_vec_prop, alpha_star, beta, X, time_to_event, event, cens,
      log_prior_n_vec, mu_alpha, s2_alpha
    );
    
    if (runif() < std::exp(log_post_prop - log_post)) {
      log_post = log_post_prop;
      n_vec = n_vec_prop;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("log_post") = log_post, Rcpp::Named("n_vec") = n_vec,
    Rcpp::Named("alpha_star") = alpha_star, Rcpp::Named("z_star") = z_star
  );
}


// [[Rcpp::export]]
Rcpp::List update_z_global(
  const double log_post, Rcpp::IntegerVector n_vec,
  Rcpp::NumericMatrix alpha_star, Rcpp::LogicalMatrix z_star,
  const arma::mat beta, const Rcpp::IntegerVector time_to_event,
  const Rcpp::IntegerVector event, const Rcpp::LogicalVector cens,
  const arma::mat X, const int n_risk, const Rcpp::Function log_prior_n_vec,
  const double mu_alpha, const double s2_alpha
) {
  // Resample a randomly selected z_star.
  
  // No. of regimes is no. of change points + 1.
  int n_regimes = n_vec.length();
  
  if (n_regimes == 1) return Rcpp::List::create(
    Rcpp::Named("alpha_star") = alpha_star, Rcpp::Named("z_star") = z_star
  );
  
  int j = (n_regimes == 1) ? 2 : sample(Rcpp::seq(1, n_regimes - 1));
  Rcpp::LogicalVector z_new = sample_z(n_risk);
  
  if (Rcpp::is_true(Rcpp::all(z_new == z_star(Rcpp::_, j)))) {
    return Rcpp::List::create(
      Rcpp::Named("alpha_star") = alpha_star, Rcpp::Named("z_star") = z_star
    );
  }
  
  Rcpp::NumericMatrix alpha_star_prop = Rcpp::clone(alpha_star);
  double log_prop_alpha = 0.0;
  
  for (int r = 0; r < n_risk; r++) {
    if (z_new[r] && !z_star(r, j)) {
      
      double alpha_new = R::rnorm(alpha_star(r, j), s_prop_alpha);
      alpha_star_prop(r, j) = alpha_new;
      
      log_prop_alpha +=
        R::dnorm(alpha_new, alpha_star(r, j), s_prop_alpha, true);
    } else if (!z_new[r] && z_star(r, j)) {
      alpha_star_prop(r, j) = alpha_star(r, j - 1);
      
      log_prop_alpha -=
        R::dnorm(alpha_star(r, j), alpha_star_prop(r, j), s_prop_alpha, true);
    }
    
    for (int t = j + 1; t < n_regimes; t++) {
      if (z_star(r, t)) break;
      alpha_star_prop(r, t) = alpha_star_prop(r, j);
    }
  }
  
  double log_post_prop = log_posterior_unaugmented(
    n_vec, alpha_star_prop, beta, X, time_to_event, event, cens,
    log_prior_n_vec, mu_alpha, s2_alpha
  );
    
  if (runif() < std::exp(log_post_prop - log_post - log_prop_alpha)) {
    z_star(Rcpp::_, j) = z_new;
    alpha_star = alpha_star_prop;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("alpha_star") = alpha_star, Rcpp::Named("z_star") = z_star
  );
}
