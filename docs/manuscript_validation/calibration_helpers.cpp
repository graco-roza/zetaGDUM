// Validation-only accelerators for the existing conditional Poisson law and
// the unconstrained, unit-weight logit estimating equation.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

static double log_add(double a, double b) {
  const double hi = std::max(a, b), lo = std::min(a, b);
  return std::isfinite(lo) ? hi + std::log1p(std::exp(lo - hi)) : hi;
}

// [[Rcpp::export]]
IntegerMatrix calibration_cp(NumericVector omega, int k, int M) {
  const int n = omega.size(), stride = k + 1;
  if (k < 1 || k > n || M < 1) stop("Invalid set size or count.");
  double largest = 0;
  for (double value : omega) {
    if (!std::isfinite(value) || value <= 0) stop("Use positive site weights.");
    largest = std::max(largest, value);
  }
  const double neg = -std::numeric_limits<double>::infinity();
  std::vector<double> dp((n + 1) * stride, neg);
  std::vector<double> take(n * stride, 0);
  dp[n * stride] = 0;
  for (int i = n - 1; i >= 0; --i) {
    dp[i * stride] = 0;
    const double log_w = std::log(omega[i]) - std::log(largest);
    for (int r = 1; r <= std::min(k, n - i); ++r) {
      const double use = log_w + dp[(i + 1) * stride + r - 1];
      const double denom = log_add(use, dp[(i + 1) * stride + r]);
      dp[i * stride + r] = denom;
      take[i * stride + r] = std::exp(use - denom);
    }
  }
  IntegerMatrix result(M, k);
  for (int b = 0; b < M; ++b) {
    int i = 0, r = k;
    while (r > 0) {
      if (n - i == r) {
        for (; r > 0; --r, ++i) result(b, k - r) = i + 1;
      } else if (2 * r <= n / 12) {
        // Invert the right tail E_r(w_i, ..., w_n) / E_r(w_start, ..., w_n).
        const double cutoff = dp[i * stride + r] + std::log(R::unif_rand());
        int lo = i, hi = n - r + 1;
        while (hi - lo > 1) {
          const int mid = (lo + hi) / 2;
          if (dp[mid * stride + r] >= cutoff) lo = mid;
          else hi = mid;
        }
        i = lo + 1;
        result(b, k - r) = i;
        --r;
      } else {
        if (R::unif_rand() < take[i * stride + r]) {
          result(b, k - r) = i + 1;
          --r;
        }
        ++i;
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector calibration_fit(NumericMatrix x, NumericVector y, NumericVector start,
                              LogicalVector keep, NumericVector w) {
  const int n = x.nrow(), p = x.ncol();
  if (p == 0 || p > 24 || y.size() != n || start.size() != p ||
      keep.size() != n || w.size() != n) stop("Mismatched calibration inputs.");
  NumericVector beta = clone(start);
  NumericVector unavailable(p, NA_REAL);
  double r[24], g[24], h[24][24];
  for (int it = 0; it < 25; ++it) {
    std::fill(g, g + p, 0);
    for (int j = 0; j < p; ++j) std::fill(h[j], h[j] + p, 0);
    int retained = 0;
    for (int i = 0; i < n; ++i) {
      if (!keep[i] || w[i] <= 0) continue;
      ++retained;
      double eta = 0;
      for (int j = 0; j < p; ++j) {
        r[j] = x(i, j);
        eta += r[j] * beta[j];
      }
      const double mu = 1 / (1 + std::exp(-eta));
      const double v = w[i] * mu * (1 - mu);
      const double s = w[i] * (y[i] - mu);
      for (int j = 0; j < p; ++j) {
        g[j] += r[j] * s;
        const double zj = r[j] * v;
        for (int l = 0; l <= j; ++l) h[j][l] += zj * r[l];
      }
    }
    if (retained < p || h[0][0] <= 0) return unavailable;
    const double tol = 1e-14 * h[0][0];
    // Cholesky factorisation and two triangular solves.
    for (int j = 0; j < p; ++j) {
      for (int l = 0; l <= j; ++l) {
        double a = h[j][l];
        for (int k = 0; k < l; ++k) a -= h[j][k] * h[l][k];
        if (j == l && a <= tol) return unavailable;
        h[j][l] = j == l ? std::sqrt(a) : a / h[l][l];
      }
      for (int l = 0; l < j; ++l) g[j] -= h[j][l] * g[l];
      g[j] /= h[j][j];
    }
    double max_step = 0, max_abs = 0;
    for (int j = p - 1; j >= 0; --j) {
      for (int l = j + 1; l < p; ++l) g[j] -= h[l][j] * g[l];
      g[j] /= h[j][j];
      max_step = std::max(max_step, std::abs(g[j]));
      max_abs = std::max(max_abs, std::abs(beta[j]));
    }
    if (!std::isfinite(max_step)) return unavailable;
    const double step = std::min(1.0, 25 / std::max(max_step, 1e-12));
    for (int j = 0; j < p; ++j) beta[j] += step * g[j];
    if (max_step <= 1e-8 * (1 + max_abs)) return beta;
  }
  return unavailable;
}
