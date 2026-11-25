#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// Set-level denominators and trait numerators from the manuscript's
// richness-normalised incidence and relative-abundance species terms.
// [[Rcpp::export]]
NumericMatrix profile_totals(NumericMatrix x, IntegerMatrix idx, NumericMatrix traits) {
  const int n = x.nrow(), p = x.ncol(), M = idx.nrow(), k = idx.ncol(), t = traits.ncol();
  NumericMatrix site(n, t + 1);
  for (int s = 0; s < p; ++s) for (int i = 0; i < n; ++i) {
    const double v = x(i, s);
    site(i, 0) += v;
    for (int j = 0; j < t; ++j) site(i, j + 1) += v * traits(s, j);
  }
  NumericMatrix out(M, 2 * (t + 1));
  std::vector<double> mean_value(t + 1), minimum(t + 1);
  for (int r = 0; r < M; ++r) {
    std::fill(mean_value.begin(), mean_value.end(), 0);
    std::fill(minimum.begin(), minimum.end(), 0);
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j <= t; ++j) mean_value[j] += site(idx(r, i) - 1, j) / k;
    }
    for (int s = 0; s < p; ++s) {
      double v = x(idx(r, 0) - 1, s);
      for (int i = 1; i < k && v > 0; ++i) v = std::min(v, x(idx(r, i) - 1, s));
      if (v == 0) continue;
      minimum[0] += v;
      for (int j = 0; j < t; ++j) minimum[j + 1] += v * traits(s, j);
    }
    for (int j = 0; j <= t; ++j) {
      out(r, j) = minimum[j] / mean_value[0];
      out(r, t + 1 + j) = (mean_value[j] - minimum[j]) / mean_value[0];
    }
  }
  return out;
}
