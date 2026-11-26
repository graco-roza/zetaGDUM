#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Number of species shared by all sites of each set (rows of idx, 1-based site indices); PA is a sites x species 0/1 matrix.
// [[Rcpp::export]]
NumericVector shared_species_cpp(IntegerMatrix PA, IntegerMatrix idx) {
  int N = PA.nrow(), S = PA.ncol(), W = (S + 63) / 64;
  std::vector<uint64_t> bits((size_t)N * W, 0ULL);
  for (int i = 0; i < N; i++) for (int s = 0; s < S; s++) if (PA(i, s)) bits[(size_t)i * W + s / 64] |= (1ULL << (s % 64));
  int n = idx.nrow(), k = idx.ncol();
  NumericVector out(n);
  std::vector<uint64_t> m(W);
  for (int r = 0; r < n; r++) {
    int s0 = idx(r, 0) - 1;
    for (int u = 0; u < W; u++) m[u] = bits[(size_t)s0 * W + u];
    for (int a = 1; a < k; a++) { int sa = idx(r, a) - 1; for (int u = 0; u < W; u++) m[u] &= bits[(size_t)sa * W + u]; }
    int cnt = 0; for (int u = 0; u < W; u++) cnt += __builtin_popcountll(m[u]);
    out[r] = cnt;
  }
  return out;
}

// Sum over species of the minimum abundance across the sites of each set.
// [[Rcpp::export]]
NumericVector min_abundance_cpp(NumericMatrix A, IntegerMatrix idx) {
  int S = A.ncol(), n = idx.nrow(), k = idx.ncol();
  NumericVector out(n);
  for (int r = 0; r < n; r++) {
    double acc = 0;
    for (int s = 0; s < S; s++) {
      double mn = A(idx(r, 0) - 1, s);
      for (int a = 1; a < k && mn > 0; a++) { double v = A(idx(r, a) - 1, s); if (v < mn) mn = v; }
      acc += mn;
    }
    out[r] = acc;
  }
  return out;
}

// Mean pairwise distance among the sites of each set; D is a full symmetric sites x sites matrix.
// [[Rcpp::export]]
NumericVector meanpair_dist_cpp(NumericMatrix D, IntegerMatrix idx) {
  int n = idx.nrow(), k = idx.ncol();
  NumericVector out(n);
  double np = k * (k - 1) / 2.0;
  for (int i = 0; i < n; i++) {
    double acc = 0;
    for (int a = 0; a < k - 1; a++) for (int b = a + 1; b < k; b++) acc += D(idx(i, a) - 1, idx(i, b) - 1);
    out[i] = acc / np;
  }
  return out;
}

// Log elementary symmetric sum of order k, updated in descending order.
static double log_add(double a, double b) {
  double hi = std::max(a, b), lo = std::min(a, b);
  if (!std::isfinite(lo)) return hi;
  return hi + std::log1p(std::exp(lo - hi));
}

// [[Rcpp::export]]
double log_esym_cpp(NumericVector log_weights, int k) {
  const double neg_inf = -std::numeric_limits<double>::infinity();
  std::vector<double> dp(k + 1, neg_inf);
  dp[0] = 0.0;
  int seen = 0;
  for (double log_weight : log_weights) {
    if (!std::isfinite(log_weight)) continue;
    for (int j = std::min(++seen, k); j > 0; --j) {
      dp[j] = log_add(dp[j], log_weight + dp[j - 1]);
    }
  }
  return dp[k];
}

// Independent Bernoulli draws conditional on selecting exactly k sites.
// Multiplying all odds by the same constant cancels under conditioning.
// [[Rcpp::export]]
IntegerMatrix cps_sets_cpp(NumericVector omega, int k, int M) {
  std::vector<int> sites;
  std::vector<double> log_weights;
  for (R_xlen_t i = 0; i < omega.size(); ++i) {
    if (!std::isfinite(omega[i]) || omega[i] < 0.0)
      stop("`omega` must contain finite, non-negative weights.");
    if (omega[i] > 0.0) {
      sites.push_back(i + 1);
      log_weights.push_back(std::log(omega[i]));
    }
  }
  int n = sites.size();
  if (k < 1 || k > n) stop("`k` must be between 1 and the number of sites with positive weight.");
  IntegerMatrix out(M, k);
  if (n == k) {
    for (int j = 0; j < k; ++j) std::fill(out.column(j).begin(), out.column(j).end(), sites[j]);
    return out;
  }
  auto range = std::minmax_element(log_weights.begin(), log_weights.end());
  double lo = -60.0 - *range.second, hi = 60.0 - *range.first;
  for (int it = 0; it < 60; ++it) {
    double mid = 0.5 * (lo + hi), expected = 0.0;
    for (double w : log_weights) expected += R::plogis(mid + w, 0.0, 1.0, true, false);
    if (expected < k) lo = mid; else hi = mid;
  }
  double shift = 0.5 * (lo + hi);
  std::vector<double> p(n);
  for (int i = 0; i < n; ++i) p[i] = R::plogis(shift + log_weights[i], 0.0, 1.0, true, false);
  std::vector<int> selected;
  selected.reserve(k + 1);
  for (int m = 0; m < M; ++m) {
    for (int tries = 0;; ++tries) {
      selected.clear();
      for (int i = 0; i < n; ++i) {
        if (R::unif_rand() < p[i]) {
          selected.push_back(sites[i]);
          if (static_cast<int>(selected.size()) > k) break;
        }
      }
      if (static_cast<int>(selected.size()) == k) break;
      if (tries >= 5000000) stop("Conditional Poisson sampling did not accept a set; check the site weights.");
    }
    for (int j = 0; j < k; ++j) out(m, j) = selected[j];
  }
  return out;
}

// Mean pairwise absolute difference of each column of the I-spline basis B (sites x 3) of a predictor z over the sites of each set.
// With the sites of a set sorted by z, the mean over pairs equals sum_j (2j - k - 1) / choose(k, 2) * B[sorted_j].
// [[Rcpp::export]]
NumericMatrix mean_pair_isp_cpp(NumericMatrix B, NumericVector z, IntegerMatrix idx) {
  int M = idx.nrow(), k = idx.ncol(), nb = B.ncol();
  NumericMatrix out(M, nb);
  double cp = 0.5 * k * (k - 1.0);
  std::vector<std::pair<double, int> > v(k);
  for (int r = 0; r < M; r++) {
    for (int j = 0; j < k; j++) { int s = idx(r, j) - 1; v[j] = std::make_pair(z[s], s); }
    std::sort(v.begin(), v.end());
    for (int b = 0; b < nb; b++) {
      double acc = 0.0;
      for (int j = 0; j < k; j++) acc += (2.0 * (j + 1) - k - 1.0) * B(v[j].second, b);
      out(r, b) = acc / cp;
    }
  }
  return out;
}

// Partition dissimilarity on relative abundance into contributions of species
// missing from at least one site and heterogeneity among all-site species.
// [[Rcpp::export]]
NumericMatrix relative_partition_cpp(const NumericMatrix& q, const IntegerMatrix& idx) {
  const int k = idx.ncol();
  NumericMatrix result(idx.nrow(), 2);
  for (int r = 0; r < idx.nrow(); ++r) {
    for (int s = 0; s < q.ncol(); ++s) {
      double minimum = q(idx(r, 0) - 1, s);
      double sum = minimum;
      for (int j = 1; j < k; ++j) {
        const double value = q(idx(r, j) - 1, s);
        minimum = std::min(minimum, value);
        sum += value;
      }
      const double mean = sum / k;
      if (minimum == 0) result(r, 0) += mean;
      else result(r, 1) += std::max(0.0, mean - minimum);
    }
  }
  return result;
}
