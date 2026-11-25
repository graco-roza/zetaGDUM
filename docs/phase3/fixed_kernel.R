source("docs/phase3/common.R")
set.seed(93331)
specs <- read.csv("docs/manuscript_validation/communities.csv")
reps <- 250L
rows <- list()
for (i in seq_len(nrow(specs))) {
  n <- specs$sites[i]
  k <- specs$high[i]
  m <- n %/% k
  B <- 800L %/% m
  for (r in seq_len(reps)) {
    site <- rnorm(n)
    h <- function(idx) rowMeans(matrix(site[idx], nrow(idx)))
    sample <- sample_sets_matched(n, k, B, m)
    values <- h(sample$idx)
    matched <- u_matched_covariance(values, sample)
    overlap <- u_overlap_covariance(n, k, h, pairs_per_overlap = 100)
    rows[[length(rows) + 1L]] <- data.frame(n, k, r, observed_complete = mean(site),
      observed_incomplete = mean(values), truth_complete = 1 / n,
      truth_incomplete = (1 - 1 / B) / n + 1 / (k * B * m),
      complete_matched = matched$complete[1], incomplete_matched = matched$incomplete[1],
      complete_overlap = overlap$complete[1],
      incomplete_overlap = (1 - 1 / B) * overlap$complete[1] + overlap$contrast[[k + 1L]][1] / (B * m))
  }
}
write_result(do.call(rbind, rows), "fixed_kernel.csv")
