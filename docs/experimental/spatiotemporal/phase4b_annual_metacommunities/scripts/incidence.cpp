#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <limits>
#include <vector>
using namespace Rcpp;

using Word = std::uint64_t;
struct Rectangle {
  std::vector<Word> shared;
  int richness = 0;
  int minimum = std::numeric_limits<int>::max();
  int count = 0;
  Rectangle(int words) : shared(words, ~Word(0)) {}
};

struct Community {
  int words, rows;
  std::vector<Word> bits;
  std::vector<int> richness;
  Community(int nr, int nc) : words((nc + 63) / 64), rows(nr),
    bits(nr * words, 0), richness(nr, 0) {}
  void add(int row, int species) {
    bits[row * words + species / 64] |= Word(1) << (species % 64);
    ++richness[row];
  }
  const Word* row(int i) const { return &bits[i * words]; }
  Rectangle rectangle(const IntegerMatrix& sets, int set, int sites) const {
    Rectangle result(words);
    for (int t = 0; t < rows / sites; ++t) {
      for (int col = 0; col < sets.ncol(); ++col) {
        const int i = sets(set, col) - 1 + t * sites;
        result.richness += richness[i];
        result.minimum = std::min(result.minimum, richness[i]);
        for (int w = 0; w < words; ++w) result.shared[w] &= row(i)[w];
      }
    }
    for (Word w : result.shared) result.count += __builtin_popcountll(w);
    return result;
  }
};

double sorensen(const Word* a, const Word* b, int words, int total) {
  if (!total) return NA_REAL;
  int shared = 0;
  for (int i = 0; i < words; ++i) shared += __builtin_popcountll(a[i] & b[i]);
  return 1.0 - 2.0 * shared / total;
}

using Metrics = std::array<double, 6>;
Metrics metrics(const Community& a, const Community& b, const Rectangle& base,
                const IntegerMatrix& sets, int set, int sites) {
  const auto now = b.rectangle(sets, set, sites);
  const int obs = sets.ncol() * b.rows / sites;
  const double total = now.richness ? 1.0 - double(now.count) * obs / now.richness : NA_REAL;
  const double turnover = now.minimum ? 1.0 - double(now.count) / now.minimum : NA_REAL;
  double local = 0;
  for (int t = 0; t < b.rows / sites; ++t) {
    for (int col = 0; col < sets.ncol(); ++col) {
      const int i = sets(set, col) - 1 + t * sites;
      local += sorensen(a.row(i), b.row(i), b.words, a.richness[i] + b.richness[i]);
    }
  }
  return {double(now.count), total, turnover, total - turnover,
          sorensen(base.shared.data(), now.shared.data(), b.words, base.count + now.count),
          local / obs};
}

Community encode(const IntegerMatrix& a) {
  Community result(a.nrow(), a.ncol());
  for (int s = 0; s < a.ncol(); ++s)
    for (int i = 0; i < a.nrow(); ++i)
      if (a(i,s)) result.add(i,s);
  return result;
}

// [[Rcpp::export]]
List annual_incidence(const IntegerMatrix& a, const IntegerMatrix& b,
                      const IntegerMatrix& sets, int sites, const List& masks) {
  const int n = sets.nrow(), draws = masks.size(), responses = 6;
  Community first = encode(a), last = encode(b);
  NumericMatrix before(n, responses), after(n, responses);
  std::vector<Rectangle> base;
  base.reserve(n);
  for (int i = 0; i < n; ++i) {
    base.push_back(first.rectangle(sets, i, sites));
    const auto x = metrics(first, first, base.back(), sets, i, sites);
    const auto y = metrics(first, last, base.back(), sets, i, sites);
    for (int j = 0; j < responses; ++j) {
      before(i,j) = x[j];
      after(i,j) = y[j];
    }
  }
  NumericVector contributions(R_xlen_t(n) * responses * 3 * draws, NA_REAL);
  contributions.attr("dim") = IntegerVector::create(n, responses, 3, draws);
  double max_error = 0;
  for (int d = 0; d < draws; ++d) {
    const IntegerMatrix mask = as<IntegerMatrix>(masks[d]);
    std::vector<Community> states(6, Community(a.nrow(), a.ncol()));
    for (int s = 0; s < a.ncol(); ++s) {
      for (int i = 0; i < a.nrow(); ++i) {
        for (int coalition = 1; coalition < 7; ++coalition)
          if ((coalition & mask(i,s)) ? b(i,s) : a(i,s))
            states[coalition - 1].add(i,s);
      }
    }
    for (int i = 0; i < n; ++i) {
      std::array<Metrics, 8> values;
      for (int coalition = 1; coalition < 7; ++coalition)
        values[coalition] = metrics(first, states[coalition - 1], base[i], sets, i, sites);
      for (int j = 0; j < responses; ++j) {
        values[0][j] = before(i,j);
        values[7][j] = after(i,j);
        bool valid = true;
        for (const auto& value : values) valid &= std::isfinite(value[j]);
        if (!valid) continue;
        double total = 0;
        for (int p = 0; p < 3; ++p) {
          double phi = 0;
          const int bit = 1 << p;
          for (int c = 0; c < 8; ++c) {
            if (c & bit) continue;
            const double weight = c == 0 || (c | bit) == 7 ? 1.0/3 : 1.0/6;
            phi += weight * (values[c | bit][j] - values[c][j]);
          }
          contributions[i + R_xlen_t(n) * (j + responses * (p + 3 * d))] = phi;
          total += phi;
        }
        max_error = std::max(max_error, std::abs(total - (after(i,j) - before(i,j))));
      }
    }
  }
  return List::create(_["before"] = before, _["after"] = after,
                      _["draws"] = contributions, _["error"] = max_error);
}
