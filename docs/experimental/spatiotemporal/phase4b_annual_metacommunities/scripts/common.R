phase <- "docs/experimental/spatiotemporal"
directory <- file.path(phase, "phase4b_annual_metacommunities")
source(file.path(phase, "phase3_load.R"))
st3_load(phase)
raw <- file.path(directory, "data", "raw")
processed <- file.path(directory, "data", "processed")
results <- file.path(directory, "results")
figures <- file.path(directory, "figures")
for (path in c(raw, processed, results, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
responses <- c("zeta_presence", "incidence", "turnover", "nestedness",
               "temporal_presence", "site_temporal_presence")
processes <- st3_processes[1:3]
write_result <- function(x, name) {
  data.table::fwrite(x, file.path(results, paste0(name, ".csv")), na = "NA")
}
read_result <- function(name) data.table::fread(file.path(results, paste0(name, ".csv")))
load_kernel <- function() {
  Rcpp::sourceCpp(file.path(directory, "scripts", "incidence.cpp"),
                  cacheDir = file.path(processed, "cpp-cache"))
}
sets_for <- function(n, k, budget, seed) {
  if (choose(n,k) <= budget) return(t(utils::combn(n,k)))
  set.seed(seed)
  t(replicate(budget, sort(sample.int(n,k))))
}
exact_zeta <- function(x, k) {
  site_sharing <- apply(x, c(1,3), all)
  occupancy <- colSums(matrix(site_sharing, dim(x)[1], dim(x)[3]))
  weight <- numeric(length(occupancy))
  valid <- occupancy >= k
  weight[valid] <- exp(lchoose(occupancy[valid],k) - lchoose(dim(x)[1],k))
  weight
}
evaluate <- function(plan, sets, allocations = NULL) {
  if (any(plan$before != (plan$before > 0)) ||
      any(plan$after != (plan$after > 0))) stop("Supply binary incidence matrices.")
  masks <- if (is.null(allocations)) list() else allocations$masks
  x <- annual_incidence(plan$before * 1L, plan$after * 1L, sets, length(plan$sites), masks)
  colnames(x$before) <- colnames(x$after) <- responses
  dimnames(x$draws) <- list(NULL, responses, processes, NULL)
  if (length(masks)) {
    x$phi <- apply(x$draws, c(1,2,3), mean)
    x$exact <- allocations$exact
  }
  x
}
