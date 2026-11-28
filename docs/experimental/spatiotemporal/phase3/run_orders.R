root <- "docs/experimental/spatiotemporal"
source(file.path(root,"phase3_load.R"))
st3_load(root)
out <- Sys.getenv("ST3_OUTPUT",file.path(root,"phase3/results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
reps <- as.integer(Sys.getenv("ST3_GRID_REPS","12"))
workers <- as.integer(Sys.getenv("ST3_WORKERS","4"))
save <- function(x,name) write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
cases <- c(st3_scenarios,"support_stress")
if (Sys.getenv("ST3_PILOT")=="1") cases <- c("opposing","persistence","support_stress")
stress_only <- Sys.getenv("ST3_SUPPORT_ONLY")=="1"
if (stress_only) cases <- "support_stress"
responses <- c("zeta_presence","zeta_abundance","zeta_dominance","incidence","turnover","nestedness",
               "abundance","occupancy","numerical","temporal_presence","temporal_abundance",
               "site_temporal_presence","site_temporal_abundance")
work <- function(task) {
  stress <- task$scenario=="support_stress"
  sim <- st3_simulate(if(stress)"filtering" else task$scenario,seed=64000L+task$replicate)
  p <- sim$panel
  if (stress) {
    # Replace the ubiquitous anchors with persistent site-specific taxa. This
    # preserves nonempty normalized assemblages but permits zero shared minima.
    n <- dim(p$x)[1]
    T <- dim(p$x)[2]
    S <- dim(p$x)[3]
    labels <- dimnames(p$x)
    labels[[3L]] <- c(labels[[3L]],paste0("anchor_",seq_len(n)))
    a <- array(0,c(n,T,S+n),labels)
    a[,,seq_len(S)] <- p$x
    a[,,1:2] <- 0
    for (j in seq_len(n)) a[j,,S+j] <- 2
    p$x <- a
  }
  basis <- st3_basis(p,sim$env)
  grid <- expand.grid(k=2:6,m=1:4)
  if (Sys.getenv("ST3_PILOT")=="1") grid <- subset(grid,k %in% c(2,6)&m %in% c(1,4))
  output <- support <- list()
  for (g in seq_len(nrow(grid))) {
    start <- proc.time()[["elapsed"]]
    k <- grid$k[g]
    m <- grid$m[g]
    sets <- st3_sample_sets(length(sim$blocks),k,26)
    windows <- st3_pair_windows(p$times,m)
    result <- lapply(seq_along(windows),function(j) {
      w <- windows[[j]]
      st3_shapley(st3_plan(p,w$before,w$after),sets,responses,exact_limit=6,draws=6,seed=task$replicate*17L+j)
    })
    a <- lapply(seq_along(responses),function(i) {
      phi <- do.call(rbind,lapply(result,function(z) z$contributions[,i,]))
      valid <- complete.cases(phi)
      before <- unlist(lapply(result,function(z) z$before[,i]))
      after <- unlist(lapply(result,function(z) z$after[,i]))
      change <- after-before
      data.frame(task,spatial_order=k,temporal_order=m,response=responses[i],process=st3_processes,
                 before=mean(before[valid]),after=mean(after[valid]),net=mean(change[valid]),
                 mean=colMeans(phi[valid,,drop=FALSE]),mean_absolute=colMeans(abs(phi[valid,,drop=FALSE])),
                 mean_positive=colMeans(pmax(phi[valid,,drop=FALSE],0)),
                 mean_homogenising=colMeans(pmax(-phi[valid,,drop=FALSE],0)),
                 absolute_mean=abs(colMeans(phi[valid,,drop=FALSE])),
                 mean_importance=colMeans(.st_ratio(abs(phi[valid,,drop=FALSE]),rowSums(abs(phi[valid,,drop=FALSE]))),na.rm=TRUE),
                 common_support=mean(valid),nonzero_support=colMeans(abs(phi[valid,,drop=FALSE])>1e-12),
                 defined_endpoints=mean(is.finite(change)),
                 allocation_mc=max(c(0,unlist(lapply(result,function(z) z$allocation_mc_se[,i,]))),na.rm=TRUE),
                 identity_error=max(vapply(result,`[[`,numeric(1),"max_efficiency_error")),row.names=NULL)
    })
    output[[g]] <- do.call(rbind,a)
    d <- st3_environment_design(p,sim$env,sets,windows,basis)
    Y <- do.call(rbind,lapply(result,st3_model_response,crossed=TRUE))
    model <- st3_fit_process_gdum(d,Y)
    pr <- st3_predict_process_gdum(model)
    available <- if(length(windows)>=4L) {
      climate <- d$X[,"position_change.climate"]/k
      high <- climate>=quantile(climate,.75)
      low <- climate<=quantile(climate,.25)
      X <- cbind(intercept=1,d$X[,model$columns,drop=FALSE])
      cvec <- colMeans(X[high,,drop=FALSE])-colMeans(X[low,,drop=FALSE])
      fit <- st3_fit_process_gdum(d,cbind(total=rowSums(Y)))
      v <- st3_contrast_uncertainty(fit,cvec,sim$blocks,draws=16,width=2L)
      unname(v$se[1,"crossed"])
    } else NA_real_
    support[[g]] <- data.frame(task,spatial_order=k,temporal_order=m,matches=length(windows),sets=nrow(sets),
      common_support=mean(complete.cases(Y)),rank=model$rank,columns=ncol(model$X),effective_df=model$effective_df,
      profiled_position=model$profile_position,profiled_spread=model$profile_spread,
      crossed_contrast_se=available,zero_share=mean(unlist(lapply(result,function(x) x$before[,"zeta_presence"]==0))),
      fit_bytes=as.numeric(object.size(model)),elapsed_seconds=proc.time()[["elapsed"]]-start,
      incompatible=mean(pr$total+unlist(lapply(result,function(z) z$before[,"abundance"]))< -1e-8 |
                          pr$total+unlist(lapply(result,function(z) z$before[,"abundance"]))>1+1e-8))
  }
  list(estimates=do.call(rbind,output),support=do.call(rbind,support))
}
tasks <- expand.grid(scenario=cases,replicate=seq_len(reps),stringsAsFactors=FALSE)
job <- function(i) work(tasks[i,])
runs <- if (.Platform$OS.type=="windows" || workers==1L) lapply(seq_len(nrow(tasks)),job) else
  parallel::mclapply(seq_len(nrow(tasks)),job,mc.cores=workers,mc.set.seed=FALSE)
if (any(vapply(runs,inherits,logical(1),"try-error"))) stop(paste(unlist(runs[vapply(runs,inherits,logical(1),"try-error")]),collapse="\n"))
for (item in c(estimates="order_processes",support="order_support")) {
  key <- if (item=="order_processes") "estimates" else "support"
  z <- do.call(rbind,lapply(runs,`[[`,key))
  path <- file.path(out,paste0(item,".csv"))
  if (stress_only && file.exists(path)) {
    old <- read.csv(path,stringsAsFactors=FALSE)
    z <- rbind(old[old$scenario!="support_stress",],z)
  }
  save(z,item)
}
