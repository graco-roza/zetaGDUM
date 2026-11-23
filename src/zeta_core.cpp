#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <algorithm>
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
