.st3_temporal_distance <- function(a, b) .st_ratio(rowSums(abs(a - b)), rowSums(a + b))

.st3_temporal_metrics <- function(a, b, idx, responses) {
  index_mean <- function(d) rowMeans(matrix(d[idx], nrow(idx), ncol(idx)))
  out <- list()
  if ("temporal_presence" %in% responses) {
    out$temporal_presence <- .st3_temporal_distance(a$shared, b$shared)
  }
  if ("temporal_abundance" %in% responses) {
    out$temporal_abundance <- .st3_temporal_distance(a$minima, b$minima)
  }
  if ("site_temporal_presence" %in% responses) {
    out$site_temporal_presence <- index_mean(.st3_temporal_distance(a$incidence, b$incidence))
  }
  if ("site_temporal_abundance" %in% responses) {
    out$site_temporal_abundance <- index_mean(.st3_temporal_distance(a$relative, b$relative))
  }
  out
}

st3_site_reorganisation <- function(plan) {
  a <- plan$before
  b <- plan$after
  n <- length(plan$sites)
  m <- nrow(plan$alignment)
  q0 <- .st3_normalize(a)$q
  q1 <- .st3_normalize(b)$q
  data.frame(site = rep(plan$sites, m), position = rep(seq_len(m), each = n),
             incidence = .st3_temporal_distance(a > 0, b > 0),
             abundance = .st3_temporal_distance(q0, q1))
}
