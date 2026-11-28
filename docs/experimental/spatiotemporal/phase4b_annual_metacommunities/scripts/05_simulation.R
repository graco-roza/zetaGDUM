source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
source(file.path(directory,"scripts","trajectory.R"))
source(file.path(directory,"scripts","simulate.R"))
load_kernel()
scenarios <- c("expansion","contraction","opposing","stable_redistribution","reversal","shocks","detection","null")
replicates <- as.integer(Sys.getenv("ST4B_SIM_REPLICATES","24"))
stopifnot(replicates > 0)
jobs <- CJ(scenario=scenarios,replicate=seq_len(replicates))
# Workers inherit the compiled kernel on Unix.
run_job <- function(i) {
  scenario <- jobs$scenario[i]
  j <- jobs$replicate[i]
  is_small <- j > replicates
  sim <- simulate_annual(scenario,3201+j,n=if (is_small) 30L else 60L,
                         species=if (is_small) 120L else 200L,years=if (is_small) 15L else 30L)
  simulation_rep(sim,scenario,j)
}
jobs <- rbind(jobs,CJ(scenario=scenarios,replicate=replicates+1L))
runs <- suppressWarnings(parallel::mclapply(seq_len(nrow(jobs)),run_job,mc.cores=4L))
stopifnot(all(vapply(runs,is.list,logical(1))))
for (name in names(runs[[1]])) write_result(rbindlist(lapply(runs,`[[`,name)),paste0("simulation_",name))
s <- read_result("simulation_summary")[replicate<=replicates & scenario!="detection" & abs(truth)>1e-10]
stopifnot(s[k==1,max(abs(accumulated-truth))] < 1e-10)
s[, `:=`(error=accumulated-truth,
         cover=abs(accumulated-truth)<=1.96*combination_se+1e-12)]
cal <- s[, .(n_rep=.N, bias=mean(error), rmse=sqrt(mean(error^2)),
             empirical_sd=sd(error), mean_mc_se=mean(combination_se),
             coverage=mean(cover), coverage_mcse=sqrt(mean(cover)*(1-mean(cover))/.N),
             coverage_low=binom.test(sum(cover),.N)$conf.int[1L],
             coverage_high=binom.test(sum(cover),.N)$conf.int[2L]),
         by=.(scenario,k,process)]
write_result(cal,"simulation_calibration")
# Representative complete 30-year trajectories evaluated across the full order grid.
support <- lapply(c("stable_redistribution","opposing","shocks","detection"),function(scenario) {
  sim <- simulate_annual(scenario,7771L)
  z <- analyse_trajectory(sim$panel,sim$blocks,budget=192L,draws=8L, exhaustive_limit=2000)
  z$diagnostics[,scenario:=scenario]
})
write_result(rbindlist(support),"simulation_order_support")
