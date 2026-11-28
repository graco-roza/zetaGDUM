root <- "docs/experimental/spatiotemporal"
source(file.path(root,"phase3_load.R"))
st3_load(root)
out <- Sys.getenv("ST3_OUTPUT",file.path(root,"phase3/results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reps <- as.integer(Sys.getenv("ST3_HOLD_REPS","24"))
workers <- as.integer(Sys.getenv("ST3_WORKERS","4"))
run <- function(scenario,r) {
  sim <- st3_simulate(scenario,seed=73000L+r)
  x <- sim$panel
  n <- length(sim$blocks)
  t <- length(x$times)
  take <- sort(sample.int(n,8L))
  results <- list()
  for (split in c("sites","forward")) {
    sites0 <- if (split=="sites") setdiff(seq_len(n),take) else seq_len(n)
    sites1 <- if (split=="sites") take else seq_len(n)
    time0 <- if (split=="forward") seq_len(t-4L) else seq_len(t)
    time1 <- if (split=="forward") (t-3L):t else seq_len(t)
    train <- st3_subset(x,sites0,time0)
    test <- st3_subset(x,sites1,time1)
    e0 <- sim$env[,time0,,drop=FALSE]
    e1 <- sim$env[,time1,,drop=FALSE]
    basis <- st3_basis(train,e0)
    for (k in c(2L,4L)) for (m in 1:2) {
      s0 <- st3_sample_sets(length(sites0),k,30)
      s1 <- st3_sample_sets(length(sites1),k,30)
      w0 <- st3_pair_windows(train$times,m)
      w1 <- st3_pair_windows(test$times,m)
      d0 <- st3_environment_design(train,e0,s0,w0,basis)
      d1 <- st3_environment_design(test,e1,s1,w1,basis)
      y0 <- st3_collect(train,s0,w0,draws=6,seed=r)
      y1 <- st3_collect(test,s1,w1,draws=6,seed=r)
      append_total <- function(z) cbind(z,total=rowSums(z))
      Y <- append_total(y0$y)
      y <- append_total(y1$y)
      f <- st3_fit_process_gdum(d0,Y)
      i <- st3_fit_process_gdum(d0,Y,interaction=TRUE)
      pp <- st3_predict_process_gdum(f,d1)
      ipp <- st3_predict_process_gdum(i,d1)
      base <- matrix(colMeans(Y),nrow(y),ncol(y),byrow=TRUE)
      flags <- function(p) mean(y1$before+p$contributions[,"total"] < -1e-8 |
                                 y1$before+p$contributions[,"total"]>1+1e-8 |
                                 rowSums(abs(p$contributions)>1+1e-8)>0)
      results[[length(results)+1L]] <- data.frame(scenario,replicate=r,split,spatial_order=k,temporal_order=m,
        process=colnames(y),rmse=sqrt(colMeans((pp$contributions-y)^2)),
        interaction_rmse=sqrt(colMeans((ipp$contributions-y)^2)),baseline_rmse=sqrt(colMeans((base-y)^2)),
        bias=colMeans(pp$contributions-y),extrapolation=mean(pp$extrapolation),
        physically_incompatible=flags(pp),support=sum(y1$valid),row.names=NULL)
    }
  }
  do.call(rbind,results)
}
cases <- if (Sys.getenv("ST3_PILOT")=="1") c("opposing","persistence","filtering") else st3_scenarios
tasks <- expand.grid(scenario=cases,r=seq_len(reps),stringsAsFactors=FALSE)
job <- function(i) run(tasks$scenario[i],tasks$r[i])
ans <- if (.Platform$OS.type=="windows" || workers==1L) lapply(seq_len(nrow(tasks)),job) else
  parallel::mclapply(seq_len(nrow(tasks)),job,mc.cores=workers,mc.set.seed=FALSE)
if (any(vapply(ans,inherits,logical(1),"try-error"))) stop(paste(unlist(ans[vapply(ans,inherits,logical(1),"try-error")]),collapse="\n"))
write.csv(do.call(rbind,ans),file.path(out,"holdout_panels.csv"),row.names=FALSE,na="")
