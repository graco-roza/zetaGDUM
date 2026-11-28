root <- "docs/experimental/spatiotemporal"
source(file.path(root,"phase3_load.R"))
st3_load(root)
out <- Sys.getenv("ST3_OUTPUT",file.path(root,"phase3/results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reps <- as.integer(Sys.getenv("ST3_REPS","200"))
truth_reps <- as.integer(Sys.getenv("ST3_TRUTH_REPS","240"))
bootstrap <- as.integer(Sys.getenv("ST3_BOOT","32"))
workers <- as.integer(Sys.getenv("ST3_WORKERS","4"))
if (anyNA(c(reps,truth_reps,bootstrap,workers)) || min(reps,truth_reps,bootstrap)<2 || workers<1) stop("Invalid replication settings.")
cases <- st3_scenarios
if (Sys.getenv("ST3_PILOT")=="1") cases <- c("opposing","null","shared_shock","filtering","persistence")
save <- function(d,name) write.csv(d,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
jobs <- function(args,fun) {
  if (.Platform$OS.type=="windows" || workers==1L) return(lapply(args,fun))
  z <- parallel::mclapply(args,fun,mc.cores=workers,mc.set.seed=FALSE)
  if (any(vapply(z,inherits,logical(1),"try-error"))) stop(paste(unlist(z[vapply(z,inherits,logical(1),"try-error")]),collapse="\n"))
  z
}
template <- st3_simulate()
set.seed(417)
sets <- st3_sample_sets(length(template$blocks),2,36)
w <- st3_pair_windows(template$panel$times)
d <- st3_environment_design(template$panel,template$env,sets,w)
climate <- d$X[,"position_change.climate"]/2
group <- ifelse(climate>=quantile(climate,.75),1,ifelse(climate<=quantile(climate,.25),-1,0))
strata <- function(y) colMeans(y[group==1,,drop=FALSE])-colMeans(y[group==-1,,drop=FALSE])
contrast <- c(intercept=0,strata(d$X))
full <- function(y) cbind(y,total=rowSums(y))
truth_job <- function(scenario) {
  z <- t(vapply(seq_len(truth_reps),function(i) {
    panel <- st3_simulate(scenario,seed=870000L+i)$panel
    a <- st3_collect(panel,sets,w,draws=8,seed=i)
    strata(full(a$y))
  },setNames(numeric(5),c(st3_processes,"total"))))
  data.frame(scenario,process=colnames(z),mc_truth=colMeans(z),mc_se=apply(z,2,sd)/sqrt(truth_reps),
             truth_replicates=truth_reps, row.names=NULL)
}
t <- do.call(rbind,jobs(cases,truth_job))
# A truly null design has zero population contrast by exchangeability.
known_zero <- t$scenario %in% c("null","shared_shock") | (t$scenario=="opposing" & t$process=="total")
t$truth <- ifelse(known_zero,0,t$mc_truth)
save(t,"independent_truth")

one <- function(scenario,r) {
  start <- proc.time()[["elapsed"]]
  sim <- st3_simulate(scenario,seed=460000L+r)
  a <- st3_collect(sim$panel,sets,w,draws=8,seed=r)
  y <- full(a$y)
  f <- st3_fit_process_gdum(d,y)
  effect <- contrast[colnames(f$X)]
  u <- st3_contrast_uncertainty(f,effect,sim$blocks,draws=bootstrap,width=2L)
  p <- st3_predict_process_gdum(f)
  estimated <- unname(drop(effect %*% f$coefficients))
  linearity <- max(abs(p$contributions[,"total"]-rowSums(p$contributions[,st3_processes,drop=FALSE])))
  composition <- st3_signed_summary(a$y)$by_combination
  estimate <- data.frame(scenario,replicate=r,process=colnames(y),observed_mean=colMeans(y),
                         mean_absolute=colMeans(abs(y)),difference=strata(y),estimate=estimated,
                         support=mean(a$valid),allocation_mc=a$allocation_mc,exact=a$exact,
                         design_rows=nrow(d$X),rank=f$rank,effective_df=f$effective_df,
                         info_position=f$profile_position,info_spread=f$profile_spread,
                         identity_error=max(a$error,linearity),elapsed_seconds=proc.time()[["elapsed"]]-start,
                         gross=mean(composition$gross),net=mean(composition$net),row.names=NULL)
  structural <- apply(abs(y),2,max) <= 1e-12
  uncertainty <- do.call(rbind,lapply(colnames(u$se),function(method) data.frame(
    scenario,replicate=r,process=colnames(y),method,estimate=estimated,
    se=ifelse(structural,NA_real_,u$se[,method]),structural_zero=structural,
    finite_temporal=u$finite_temporal,finite_crossed=u$finite_crossed,row.names=NULL)))
  list(estimate=estimate,uncertainty=uncertainty)
}
tasks <- expand.grid(scenario=cases, batch=seq_len(min(8L,reps)),stringsAsFactors=FALSE)
runs <- jobs(seq_len(nrow(tasks)),function(i) {
  r <- tasks[i,]
  lapply(seq.int(r$batch,reps,by=min(8L,reps)),function(j) one(r$scenario,j))
})
runs <- unlist(runs,recursive=FALSE)
raw <- do.call(rbind,lapply(runs,`[[`,"estimate"))
u <- merge(do.call(rbind,lapply(runs,`[[`,"uncertainty")),t,by=c("scenario","process"),sort=FALSE)
u$error <- u$estimate-u$truth
# Pointwise interval evaluation also reports the Monte Carlo error on truth.
u$covered <- abs(u$error)<=1.96*u$se
u$reject_zero <- abs(u$estimate)>1.96*u$se
save(raw,"main_panels")
save(u,"interval_panels")

keys <- c("scenario","process","method")
summary <- do.call(rbind,lapply(split(u,u[keys],drop=TRUE),function(z) {
  q <- z[is.finite(z$se),]
  coverage <- mean(q$covered)
  rejection <- mean(q$reject_zero)
  n <- nrow(q)
  data.frame(z[1L,keys],panels=n,structural_zero_panels=sum(z$structural_zero),truth=z$truth[1],truth_mc_se=z$mc_se[1],
             bias=mean(q$error),rmse=sqrt(mean(q$error^2)),sampling_sd=sd(q$estimate),
             mean_se=mean(q$se),coverage,coverage_mc_se=sqrt(coverage*(1-coverage)/n),
             reject_zero=rejection,rejection_mc_se=sqrt(rejection*(1-rejection)/n),
             wrong_sign=if(z$truth[1]==0)NA_real_ else mean(sign(q$estimate)!=sign(z$truth[1])),
             row.names=NULL)
}))
save(summary,"interval_summary")
