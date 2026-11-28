source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
load_kernel()

checks <- list()
check <- function(name, error, tolerance = 1e-12) {
  stopifnot(is.finite(error), error <= tolerance)
  checks[[name]] <<- data.frame(test = name, error, tolerance)
}
set.seed(701)
for (m in c(1L, 3L)) {
  a <- array(as.integer(runif(7 * m * 2 * 137) < 0.34), c(7, m * 2, 137))
  a[1, , ] <- 0L
  panel <- st_panel(a)
  plan <- st3_plan(panel, seq_len(m), m + seq_len(m))
  for (k in 1:4) {
    sets <- t(combn(7,k))
    allocations <- st3_allocations(plan, draws = 6, seed = 99)
    candidate <- evaluate(plan, sets, allocations)
    reference <- st3_shapley(plan, sets, responses, draws = 6, seed = 99)
    error <- function(x,y) max(c(0, abs(x-y)), na.rm = TRUE)
    stopifnot(identical(is.na(candidate$phi), is.na(reference$contributions[,,1:3,drop=FALSE])))
    check(sprintf("four-block dummy equivalence m%d k%d",m,k),
          error(candidate$phi, reference$contributions[,,1:3,drop=FALSE]))
    check(sprintf("efficiency m%d k%d",m,k), candidate$error)
    check(sprintf("dummy is zero m%d k%d",m,k), error(reference$contributions[,,4L],0))
    check(sprintf("exhaustive vs exact zeta m%d k%d",m,k),
          abs(mean(candidate$before[,1]) - sum(exact_zeta(a[, seq_len(m), , drop=FALSE],k))))
  }
}

# A moving taxon and a ubiquitous anchor. Mean sharing and spatial dissimilarity
# are invariant, but the four sites touched by the swap change in composition.
a <- array(0L,c(6,2,2))
a[,,1] <- 1L
a[1:2,1,2] <- a[3:4,2,2] <- 1L
plan <- st3_plan(st_panel(a),1,2)
for (k in 1:6) {
  z <- evaluate(plan, t(combn(6,k)), st3_allocations(plan))
  check(paste("redistribution preserves mean raw zeta k",k), abs(mean(z$after[,1] - z$before[,1])))
  check(paste("pure redistribution has no expansion/contraction k",k), max(abs(z$phi[,,1:2]),na.rm=TRUE))
}
reorg <- evaluate(plan, t(combn(6,1)))
check("local temporal distance retains redistribution", abs(mean(reorg$after[,6])-2/9))
# Match the package's native Sørensen–Simpson special case and higher-order response.
set.seed(911)
a <- array(runif(12*2*30) > .61,c(12,2,30))
a[,,1] <- TRUE
plan <- st3_plan(st_panel(a),1,2)
native_curve <- zetaGDUM::zeta_curve(plan$before*1L,kmax=6,type="presence")
for (k in 2:6) {
  sets <- t(combn(12,k))
  z <- evaluate(plan,sets)
  for (component in c("total","turnover","nestedness")) {
    native <- zetaGDUM::ugdm_response(plan$before*1L,sets,component)$y
    name <- if (component == "total") "incidence" else component
    check(paste("native",component,"k",k),max(abs(native-z$before[,name])))
  }
  check(paste("native zeta k",k),abs(mean(z$before[,1])-native_curve$zeta[k]))
}
write_result(data.table::rbindlist(checks), "mathematical_validation")
