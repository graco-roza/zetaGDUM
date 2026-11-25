#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Cumulative threshold values indexed by a site's positive abundance rank.
// The incidence calculation uses the same representation with one threshold.
// [[Rcpp::export]]
List conditional_zeta_cpp(NumericMatrix x, int k, bool abundance,
                          bool pairs, List include) {
  const int n = x.nrow(), p = x.ncol();
  std::vector<std::vector<int>> queries(k + 1);
  for (R_xlen_t j = 0; j < include.size(); ++j) {
    IntegerVector a = include[j];
    queries[a.size()].push_back(j);
  }
  std::vector<std::vector<double>> prob(k + 1);
  for (int a = 0; a <= k; ++a) {
    if (a > 1 && !(a == 2 && pairs) && queries[a].empty()) continue;
    prob[a].resize(n + 1, 0);
    const double den = R::lchoose(n - a, k - a);
    for (int o = k; o <= n; ++o) {
      prob[a][o] = std::exp(R::lchoose(o - a, k - a) - den);
    }
  }
  double mean = 0;
  NumericVector site(n), forced(include.size());
  NumericMatrix pair(pairs ? n : 0, pairs ? n : 0);
  std::vector<int> rank(n);
  std::vector<std::pair<double, int>> occupied;
  occupied.reserve(n);
  std::vector<double> values, cumulative;
  for (int s = 0; s < p; ++s) {
    if (s % 32 == 0) checkUserInterrupt();
    occupied.clear();
    for (int i = 0; i < n; ++i) {
      if (x(i, s) > 0) occupied.emplace_back(abundance ? x(i, s) : 1, i);
    }
    const int o = occupied.size();
    if (o < k) continue;
    if (abundance) std::sort(occupied.begin(), occupied.end());
    values.assign(1, 0);
    std::fill(rank.begin(), rank.end(), 0);
    std::vector<int> thresholds;
    // A tied group uses its first index, so all tied sites enter the threshold.
    for (int j = 0; j < o; ++j) {
      const double value = occupied[j].first;
      if (value > values.back()) {
        values.push_back(value);
        thresholds.push_back(o - j);
      }
      rank[occupied[j].second] = values.size() - 1;
    }
    for (int a = 0; a <= k; ++a) {
      if (prob[a].empty()) continue;
      cumulative.assign(values.size(), 0);
      for (size_t j = 1; j < values.size(); ++j) {
        cumulative[j] = cumulative[j - 1] + (values[j] - values[j - 1]) * prob[a][thresholds[j - 1]];
      }
      if (a == 0) {
        mean += cumulative.back();
      } else if (a == 1) {
        for (const auto& entry : occupied) site[entry.second] += cumulative[rank[entry.second]];
      } else if (a == 2 && pairs) {
        // For incidence this is one weighted sparse outer product per species.
        for (int i = 0; i < o; ++i) for (int j = 0; j < i; ++j) {
          const int u = occupied[i].second, v = occupied[j].second;
          pair(std::max(u, v), std::min(u, v)) += cumulative[std::min(rank[u], rank[v])];
        }
      }
      for (int j : queries[a]) {
        IntegerVector members = include[j];
        int r = cumulative.size() - 1;
        for (int m : members) {
          r = std::min(r, rank[m - 1]);
          if (r == 0) break;
        }
        forced[j] += cumulative[r];
      }
    }
  }
  if (pairs) {
    for (int i = 0; i < n; ++i) {
      pair(i, i) = NA_REAL;
      for (int j = 0; j < i; ++j) pair(j, i) = pair(i, j);
    }
  }
  return List::create(_["mean"] = mean, _["site"] = site, _["pair"] = pair, _["forced"] = forced);
}
