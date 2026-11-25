suppressPackageStartupMessages(library(zetaGDUM))

iid_sets <- function(n, k, m) {
  t(vapply(seq_len(m), function(i) sample.int(n, k), integer(k)))
}

logistic_loss <- function(y, mu, w) {
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  -weighted.mean(y * log(mu) + (1 - y) * log1p(-mu), w)
}

write_result <- function(x, file) {
  write.csv(x, file.path("docs/phase3/results", file), row.names = FALSE, na = "")
}
