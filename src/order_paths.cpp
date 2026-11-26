#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// [[Rcpp::export]]
Rcpp::List path_kernels_cpp(Rcpp::List paths, Rcpp::NumericMatrix counts,
                            Rcpp::NumericMatrix relative, Rcpp::NumericMatrix spread_basis,
                            Rcpp::NumericMatrix position, Rcpp::NumericMatrix distance) {
  using namespace Rcpp;
  const int ns = counts.ncol(), nb = spread_basis.ncol(), nv = position.ncol();
  const int last = paths.size();
  std::vector<IntegerMatrix> input;
  input.reserve(last);
  for (int j = 0; j < last; ++j) input.emplace_back(as<IntegerMatrix>(paths[j]));
  const int M = input[0].nrow();
  List out(last);
  std::vector<NumericMatrix> ans;
  for (int j = 0; j < last; ++j) ans.emplace_back(M, 5 + nb + nv);
  std::vector<int> live;
  std::vector<double> minimum(ns), sums(nb + nv);
  std::vector<double> richness(counts.nrow(), 0);
  std::vector<bool> valid(counts.nrow(), false);
  for (int s = 0; s < ns; ++s) for (int i = 0; i < counts.nrow(); ++i) {
    if (counts(i,s) > 0) {
      ++richness[i];
      valid[i] = true;
    }
  }
  for (int r = 0; r < M; ++r) {
    if (r % 512 == 0) checkUserInterrupt();
    double geo = 0, rich = 0;
    bool positive = true;
    for (int target = 0; target < last; ++target) {
      const IntegerMatrix& path = input[target];
      int start = target == 0 ? 0 : input[target - 1].ncol();
      for (int q = 0; q < start; ++q) {
        if (input[target - 1](r, q) != path(r, q)) {
          start = 0;
          break;
        }
      }
      if (start == 0) {
        live.clear();
        std::fill(sums.begin(), sums.end(), 0);
        geo = rich = 0;
        positive = true;
      }
      for (int q = start; q < path.ncol(); ++q) {
        const int i = path(r, q) - 1;
        rich += richness[i];
        positive = positive && valid[i];
        if (q == 0) {
          for (int s = 0; s < ns; ++s) if (counts(i,s) > 0) {
            live.push_back(s);
            minimum[s] = relative(i,s);
          }
        } else {
          size_t kept = 0;
          for (int s : live) if (counts(i,s) > 0) {
            live[kept++] = s;
            minimum[s] = std::min(minimum[s], relative(i,s));
          }
          live.resize(kept);
        }
        for (int b = 0; b < nv; ++b) sums[nb + b] += position(i,b);
        for (int t = 0; t < q; ++t) {
          const int j = path(r,t) - 1;
          for (int b = 0; b < nb; ++b) sums[b] += std::abs(spread_basis(i,b) - spread_basis(j,b));
          if (distance.nrow()) geo += distance(i,j);
        }
      }
      NumericMatrix& a = ans[target];
      const double k = path.ncol(), pairs = k * (k - 1) / 2;
      double shared_abundance = 0;
      for (int s : live) shared_abundance += minimum[s];
      a(r,0) = live.size();
      a(r,1) = rich / k;
      a(r,2) = rich > 0 ? 1 - a(r,0) / a(r,1) : NA_REAL;
      a(r,3) = positive ? std::max(0., 1. - shared_abundance) : NA_REAL;
      a(r,4) = geo / pairs;
      for (int b = 0; b < nb; ++b) a(r,5 + b) = sums[b] / pairs;
      for (int b = 0; b < nv; ++b) a(r,5 + nb + b) = sums[nb + b];
    }
  }
  for (int j = 0; j < last; ++j) out[j] = ans[j];
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix score_inclusion_cpp(Rcpp::IntegerMatrix idx, Rcpp::NumericMatrix X,
                                        Rcpp::NumericVector weighted_residual,
                                        Rcpp::NumericVector weights, int n) {
  Rcpp::NumericMatrix out(n + 1, X.ncol() + 1);
  for (int r = 0; r < X.nrow(); ++r) {
    const double wi = weights[r], ri = weighted_residual[r];
    out(n, 0) += wi;
    for (int j = 0; j < idx.ncol(); ++j) out(idx(r,j) - 1, 0) += wi;
    for (int p = 0; p < X.ncol(); ++p) {
      const double v = X(r,p) * ri;
      out(n, p + 1) += v;
      for (int j = 0; j < idx.ncol(); ++j) out(idx(r,j) - 1, p + 1) += v;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix site_score_kernel_cpp(Rcpp::IntegerMatrix sites,
                                          Rcpp::NumericMatrix design,
                                          Rcpp::NumericVector y,
                                          Rcpp::NumericVector effective,
                                          Rcpp::NumericVector theta, int n) {
  const int k = sites.nrow(), p = design.nrow(), m = sites.ncol();
  const int stride = n + 1;
  const int* labels = INTEGER(sites);
  const double* x = REAL(design);
  const double* observed = REAL(y);
  const double* weight = REAL(effective);
  const double* coefficients = REAL(theta);
  Rcpp::NumericMatrix result(stride, p);
  double* totals = REAL(result);

  for (int r = 0; r < m; ++r) {
    double eta = 0;
    for (int j = 0; j < p; ++j) eta += x[j] * coefficients[j];
    const double exp_tail = std::exp(-std::abs(eta));
    const double mu = eta >= 0 ? 1. / (1. + exp_tail) : exp_tail / (1. + exp_tail);
    const double residual = weight[r] * (observed[r] - mu);
    for (int j = 0; j < p; ++j) {
      double* column = totals + j * stride;
      const double contribution = x[j] * residual;
      column[n] += contribution;
      for (int i = 0; i < k; ++i) column[labels[i] - 1] += contribution;
    }
    labels += k;
    x += p;
  }
  return result;
}
