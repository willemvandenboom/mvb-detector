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
