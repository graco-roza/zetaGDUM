# bci_lib.R: local species contributions from the package's exact conditional zeta (per-species calls of the compiled routine used by conditional_zeta()).
# For site i and order k: contribution of species s to expected local dissimilarity  c_is(k) = q_is - zeta_{s|i}(k)  (relative abundance), or E[set-mean occurrence of s | i in set] - zeta_{s|i}(k) (incidence, unnormalised: mean set richness minus shared richness);
# these sum exactly to the expected dissimilarity returned by zeta_lczd(). Split into occupancy part  O_is = q_is (1 - P_s)  (P_s = prob. the other k-1 sites all contain s) and
# heterogeneity part  H_is = q_is P_s - zeta_{s|i}  (the shared-species abundance deficit).
.libPaths(c(Sys.getenv("ZG_LIB6"), .libPaths())); suppressMessages({library(zetaGDUM); library(dplyr); library(tidyr)})
species_local <- function(x, k, type = c("relative", "incidence")) { type <- match.arg(type); n <- nrow(x); rel <- type == "relative"; q <- if (rel) x / rowSums(x) else 1 * (x > 0); o <- colSums(x > 0)
  Z <- vapply(seq_len(ncol(x)), function(s) zetaGDUM:::conditional_zeta_cpp(q[, s, drop = FALSE], as.integer(k), rel, FALSE, list())$site, numeric(n)); P <- exp(lchoose(o - 1, k - 1) - lchoose(n - 1, k - 1)); P[o < k] <- 0
  if (rel) { C <- q - Z; Oc <- sweep(q, 2, 1 - P, "*"); Oc[q == 0] <- 0; H <- C - Oc } else { marg <- (q + (k - 1) * (matrix(o, n, ncol(x), byrow = TRUE) - q) / max(1, n - 1)) / k; C <- marg - Z; Oc <- C; H <- C * 0 }; dimnames(C) <- dimnames(Oc) <- dimnames(H) <- dimnames(x); list(C = C, O = Oc, H = H, P = P, q = q) }
