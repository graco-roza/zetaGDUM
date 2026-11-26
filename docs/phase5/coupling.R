source("docs/phase5/common.R")
args <- commandArgs(TRUE)
case <- as.integer(args[1])
n <- c(14L,70L,150L)[case]
s <- simulate_metacommunity(n, c(30,58,200)[case], breadth = c(.9,.36,.64)[case],
  spatial = .55, abundance_mean = 1.4, seed = 1501 + case)
s$comm <- cbind(s$comm, diag(1L * (rowSums(s$comm) == 0)))
orders <- if (case == 1) c(3,5) else c(round(.05*n),round(.12*n))
dm <- coords_dist(s$coords)
nomen <- effect_names[c(1:7,11:14)]
full_summary <- function(a, component) {
  d <- model_data(a, component, "B")
  f <- fit5(d, "B")
  u <- if (f$converged && is.finite(f$D2)) ugdm_local_distinctiveness(s$comm,f,d$design,
    if (component == "relative") "relative" else "normalised",d$sampling_weights)
  setNames(effect5(d,f,s$env,u)[c(1:7,11:14)],nomen)
}
Mref <- if (case == 1L) choose(n,max(orders)) else 40000L
base <- ugdm_path_data(s$comm,s$env,sample_sets_coupled(s$comm,orders,Mref,
  exact_limit = if(case==1L) Mref else 0,seed=531),dm)
run <- function(M, method, coupled) {
  paths <- sample_sets_coupled(s$comm, orders, M, sampler = method, coupled = coupled, exact_limit = 0)
  a <- ugdm_path_data(s$comm, s$env, paths, dm, ref = base[[1L]]$design)
  list(paths = paths, a = a)
}
set.seed(1752 + case)
output <- costs <- list()
adaptive_only <- length(args) > 1L && args[2] == "adaptive"
if (!adaptive_only) {
for (component in c("incidence","relative")) {
  target <- lapply(base,full_summary,component=component)
  for (budget in c(600L,2400L)) for (method in c("uniform","species","site")) for (coupled in c(FALSE,TRUE)) {
    for (r in seq_len(24L)) {
      tick <- sum(proc.time()[1:2])
      dat <- run(budget, method, coupled)
      fit <- lapply(dat$a, full_summary, component=component)
      cpu <- sum(proc.time()[1:2]) - tick
      for (j in 1:3) {
        estimate <- if (j == 3L) fit[[2]] - fit[[1]] else fit[[j]]
        truth <- if (j == 3L) target[[2]] - target[[1]] else target[[j]]
        output[[length(output)+1L]] <- data.frame(case,n,order_low=orders[1],order_high=orders[2],
          component, budget, sampler = method, coupled, r, comparison = c("low","high","change")[j],
          effect=nomen,estimate=unname(estimate),reference=unname(truth),cpu)
      }
    }
  }
}
path <- sample_sets_coupled(s$comm,orders,8000,seed=13,exact_limit=0)
for (r in seq_len(16L)) {
  tick <- sum(proc.time()[1:2])
  a <- ugdm_path_data(s$comm,s$env,path,dm,base[[1]]$design)
  bundled <- sum(proc.time()[1:2]) - tick
  tick <- sum(proc.time()[1:2])
  standalone <- lapply(path, function(p) {
    d <- ugdm_design(s$env,p$idx,dm,ref=base[[1]]$design)
    r <- ugdm_response(s$comm,p$idx,"relative")
    i <- ugdm_response(s$comm,p$idx,"total")
    ugdm_model_design(d,"B",s$env)
  })
  conventional <- sum(proc.time()[1:2]) - tick
  costs[[r]] <- data.frame(case,n,order_low=orders[1],order_high=orders[2],budget=8000,r,
                             path_cpu=bundled,standalone_cpu=conventional)
}
write5(do.call(rbind,output),paste0("coupling_",case))
write5(do.call(rbind,costs),paste0("path_cost_",case))
}
if (case %in% 2:3) {
  callback <- function(M, method) {
    coupled <- startsWith(method,"coupled")
    family <- substr(method,regexpr("_",method)+1,nchar(method))
    q <- lapply(run(M, family, coupled)$a,full_summary,component="relative")
    ans <- c(q[[2]][1:2],setNames((q[[2]]-q[[1]])[1:2], c("spread_change","position_change")))
    ans
  }
  effort <- ugdm_contrast_budget(callback,c("coupled_uniform","independent_uniform","coupled_species","coupled_site"),
    c(500,1500,3000),c(spread=.010,position=.014,spread_change=.010,position_change=.014),
    nrep=12,max_budget=15000,seed=571)
  write5(cbind(sampler=effort$sampler,effort$contrasts), paste0("adaptive_confirmation_",case))
  write5(effort$pilots,paste0("adaptive_pilots_",case))
}
