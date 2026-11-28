# 70_splits.R: 30 independent random site-level 70/30 splits (Latombe et al. 2017 design) for each EPA community and pool.
# Training and test site pools are disjoint within a split. The same split is used by every model, order and response representation
# (splits depend only on the community pool's site ids and the split index). In the shared pool the macroinvertebrate and zooplankton
# tables contain the same lakes, so they receive identical splits. Run from final/ with ZG_LIB6 set.
source("common.R"); library(dplyr); NS <- 30L; PROP <- .70; BASE_SEED <- 910000L
rows <- list()
for (pool in c("shared", "full")) for (cc in c("BM", "ZP")) { ids <- sort(POOLS[[pool]][[cc]]$env$id); n <- length(ids); ntr <- round(PROP * n)
  for (s in seq_len(NS)) { set.seed(BASE_SEED + s); tr <- sample(ids, ntr)
    rows[[length(rows) + 1]] <- data.frame(pool = pool, community = COMM[[cc]], split = s, site_id = ids, role = ifelse(ids %in% tr, "train", "test")) } }
S <- bind_rows(rows); stopifnot(all(S |> group_by(pool, community, split) |> summarise(ok = sum(role == "train") == round(PROP * n()), .groups = "drop") |> pull(ok)))
write.csv(S, "revision_random/tables/random_splits.csv", row.names = FALSE)
chk <- S |> group_by(pool, community, split) |> summarise(n_train = sum(role == "train"), n_test = sum(role == "test"), .groups = "drop"); write.csv(chk, "revision_random/tables/random_splits_sizes.csv", row.names = FALSE)
# shared pool: identical splits across taxa
a <- S |> filter(pool == "shared", community == COMM[["BM"]]) |> arrange(split, site_id); b <- S |> filter(pool == "shared", community == COMM[["ZP"]]) |> arrange(split, site_id); stopifnot(identical(a$role, b$role))
# overlap between the 30 splits (test-set reuse is expected in random resampling; reported for transparency)
print(chk |> group_by(pool, community) |> summarise(n_train = unique(n_train), n_test = unique(n_test), .groups = "drop"))
