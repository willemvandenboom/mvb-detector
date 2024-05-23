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


//' Parallel sampling of mixture indicators
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerMatrix sample_r_i_mat(
  const Rcpp::NumericMatrix U, const Rcpp::NumericMatrix z_min_eta,
  const Rcpp::NumericVector m_mix, const Rcpp::NumericVector tmp1,
  const Rcpp::NumericVector tmp2
) {
  int n_risk = z_min_eta.nrow(), n_long_k = z_min_eta.ncol();

  // Allocate the output matrix.
  Rcpp::IntegerMatrix r_i_mat(n_risk, n_long_k);
  
  // `sample_r_i` functor (Pass input and output matrices.)
  sample_r_i_worker worker(U, z_min_eta, m_mix, tmp1, tmp2, r_i_mat);
  
  // Call `parallelFor` to do the work.
  RcppParallel::parallelFor(0, U.length(), worker);
  
  return r_i_mat;  // Return the output matrix.
}
