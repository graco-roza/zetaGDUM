.st3_responses <- c("zeta_presence", "zeta_abundance", "zeta_dominance", "incidence",
                    "turnover", "nestedness", "abundance", "occupancy", "numerical",
                    "temporal_presence", "temporal_abundance", "site_temporal_presence",
                    "site_temporal_abundance")

.st3_community_matrix <- function(x) {
  if (!is.matrix(x) || !(is.numeric(x) || is.logical(x)) || !length(x) ||
      any(!is.finite(x) | x < 0)) stop("Use a finite nonnegative assemblage-by-species matrix.")
  x
}

.st3_normalize <- function(a) {
  total <- rowSums(a)
  list(q = a / ifelse(total > 0, total, NA_real_), total = total)
}

st3_abundance_partition <- function(a) {
  a <- .st3_community_matrix(a)
  norm <- .st3_normalize(a)
  x <- a > 0
  r <- colSums(x)
  q <- norm$q
  positive_mean <- colSums(q) / ifelse(r > 0, r, 1)
  minimum <- apply(ifelse(x, q, Inf), 2L, min)
  minimum[r == 0] <- 0
  p <- r / nrow(a)
  h <- r == nrow(a)
  occ <- (p - h) * minimum
  num <- p * pmax(positive_mean - minimum, 0)
  if (any(norm$total == 0)) {
    positive_mean[] <- minimum[] <- occ[] <- num[] <- NA_real_
  }
  species <- data.frame(species = if (is.null(colnames(a))) seq_len(ncol(a)) else colnames(a),
                        occurrences = r, p, shared = h, positive_mean, minimum,
                        occupancy = occ, numerical = num, total = occ + num,
                        row.names = NULL)
  list(species = species, diversity = c(occupancy = sum(occ), numerical = sum(num),
                                      abundance = sum(occ + num)))
}

# One row per selected rectangle; the columns of idx are the matched site-time positions.
.st3_rectangle_state <- function(a, idx, responses) {
  need_q <- any(responses %in% c("zeta_abundance", "abundance", "occupancy", "numerical",
                                 "temporal_abundance", "site_temporal_abundance"))
  need_d <- "zeta_dominance" %in% responses
  norm <- if (need_q) .st3_normalize(a) else list(total = rowSums(a))
  x <- a > 0
  r <- matrix(0L, nrow(idx), ncol(a))
  k <- ncol(idx)
  qsum <- qmin <- dominance <- NULL
  if (need_q) {
    qsum <- r * 0.0
    qmin <- qsum + Inf
  }
  if (need_d) {
    d <- zetaGDUM::dominance_matrix(a, .8) > 0
    dominance <- r == 0L
  }
  rich <- rowSums(x)
  richness_min <- rep(Inf, nrow(idx))
  empty <- rep(FALSE, nrow(idx))
  for (j in seq_len(k)) {
    ii <- idx[, j]
    incidence <- x[ii, , drop = FALSE]
    r <- r + incidence
    richness_min <- pmin(richness_min, rich[ii])
    if (need_d) dominance <- dominance & d[ii, , drop = FALSE]
    if (need_q) {
      q <- norm$q[ii, , drop = FALSE]
      qsum <- qsum + q
      qmin[] <- pmin(qmin, ifelse(incidence, q, Inf))
      empty <- empty | norm$total[ii] == 0
    }
  }
  shared <- r == k
  h <- rowSums(shared)
  result <- list(zeta_presence = h,
                 incidence = .st_ratio(rowSums(r) / k - h, rowSums(r) / k),
                 turnover = .st_ratio(richness_min - h, richness_min))
  result$nestedness <- result$incidence - result$turnover
  if (need_d) result$zeta_dominance <- rowSums(dominance)
  minima <- NULL
  if (need_q) {
    qmin[r == 0] <- 0
    minima <- shared * qmin
    minima[empty, ] <- NA_real_
    result$zeta_abundance <- rowSums(minima)
    result$abundance <- 1 - result$zeta_abundance
    result$occupancy <- rowSums((r / k - shared) * qmin)
    result$occupancy[empty] <- NA_real_
    result$numerical <- rowSums(pmax(qsum - r * qmin, 0) / k)
    result$numerical[empty] <- NA_real_
  }
  list(metrics = result, shared = shared, minima = minima, incidence = x, relative = norm$q)
}

.st3_site_indices <- function(sets, n, m) {
  do.call(cbind, lapply(seq_len(m), function(j) sets + (j - 1L) * n))
}
