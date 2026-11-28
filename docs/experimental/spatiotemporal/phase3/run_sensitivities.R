root <- "docs/experimental/spatiotemporal"
source(file.path(root,"phase3_load.R"))
st3_load(root)
out <- Sys.getenv("ST3_OUTPUT",file.path(root,"phase3/results"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
save <- function(x,name) write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")

allocation <- list()
for (r in 1:24) {
  set.seed(51000+r)
  n <- 8L
  a <- cbind(exp(rnorm(n)),exp(rnorm(n)),c(runif(3,1,4),rep(0,5)))
  b <- a
  b[1:2,3] <- 0
  b[4:7,3] <- runif(4,1,4)
  b[,1] <- a[,1]*exp(rnorm(n,.1,.25))
  p <- st3_plan(aperm(array(cbind(a,b),c(n,3,2)),c(1,3,2)),1,2)
  idx <- t(combn(n,3))
  truth <- st3_shapley(p,idx,c("abundance","occupancy","numerical"),retain_draws=TRUE)
  for (draws in c(4,16,64)) {
    z <- st3_shapley(p,idx,c("abundance","occupancy","numerical"),exact_limit=1L,draws=draws,
                     seed=r+1500L,retain_draws=TRUE)
    for (j in seq_along(st3_processes)) {
      reference <- vapply(truth$draws,function(q) mean(q[,"abundance",j]),numeric(1))
      drawn <- vapply(z$draws,function(q) mean(q[,"abundance",j]),numeric(1))
      allocation[[length(allocation)+1L]] <- data.frame(replicate=r,draws,process=st3_processes[j],
        truth=mean(reference),estimate=mean(drawn),mc_se=sd(drawn)/sqrt(draws),
        allocation_range=diff(range(reference)),
        largest_local_range=max(apply(vapply(truth$draws,function(q) q[,"abundance",j],numeric(nrow(idx))),1,function(x) diff(range(x)))),
        common_support=mean(z$valid[,"abundance"]))
    }
  }
}
allocation <- do.call(rbind,allocation)
allocation$error <- allocation$estimate-allocation$truth
save(allocation,"allocation_sensitivity")

spacing <- list()
for (r in 1:24) {
  for (irregular in c(FALSE,TRUE)) {
    z <- st3_simulate("persistence",seed=58000+r,irregular=irregular)
    sets <- st3_sample_sets(length(z$blocks),3,34)
    for (window in c("aligned","initial_shift","gap")) {
      before <- if (window=="initial_shift") 2:4 else 1:3
      after <- if (window=="gap") 12:14 else 9:11
      p <- st3_plan(z$panel,before,after)
      a <- st3_shapley(p,sets,c("zeta_presence","abundance","temporal_presence","site_temporal_presence"),
                        draws=8,seed=r)
      for (response in colnames(a$delta)) {
        x <- a$contributions[,response,]
        spacing[[length(spacing)+1L]] <- data.frame(replicate=r,irregular,window,response,
          process=colnames(x),attributed=colMeans(x,na.rm=TRUE),support=mean(a$valid[,response]),
          before_duration=diff(range(p$alignment$initial)),after_duration=diff(range(p$alignment$final)),
          mean_lag=mean(p$alignment$lag),max_spacing_difference=max(abs(diff(p$alignment$initial)-diff(p$alignment$final))),
          row.names=NULL)
      }
    }
  }
}
save(do.call(rbind,spacing),"window_sensitivity")

if (requireNamespace("ecopart",quietly=TRUE)) {
  sim <- st3_simulate("filtering",seed=19101)
  sets <- st3_sample_sets(length(sim$blocks),2,36)
  cases <- list()
  for (j in c(1L,3L,5L,7L)) {
    p <- st3_plan(sim$panel,j,j+1L)
    z <- st3_shapley(p,sets,c("incidence","abundance"),draws=12,seed=j)
    for (i in seq_len(nrow(sets))) for (measure in c("incidence","abundance")) {
      a <- p$before[sets[i,],,drop=FALSE]
      b <- p$after[sets[i,],,drop=FALSE]
      if (measure=="abundance") {
        a <- .st3_normalize(a)$q
        b <- .st3_normalize(b)$q
      }
      parts <- ecopart::ecopart.pair(a,b,index=if(measure=="incidence")"sorensen" else "bray-curtis",components="four")
      phi <- z$contributions[i,measure,]
      cases[[length(cases)+1L]] <- data.frame(contrast=j,set=i,measure,
        observed=z$delta[i,measure],shapley_net=sum(phi),tatsumi_net=sum(parts),
        shapley_differentiating=sum(pmax(phi,0)),shapley_homogenising=sum(pmax(-phi,0)),
        tatsumi_loss_homogenising=unname(-parts[1]),tatsumi_loss_differentiating=unname(parts[2]),
        tatsumi_gain_homogenising=unname(-parts[3]),tatsumi_gain_differentiating=unname(parts[4]))
    }
  }
  save(do.call(rbind,cases),"pairwise_comparator")
}

curves <- list()
for (r in 1:40) {
  sim <- st3_simulate("filtering",seed=15200+r)
  p <- sim$panel
  basis <- st3_basis(p,sim$env)
  k <- 3L
  sets <- st3_sample_sets(length(sim$blocks),k,46)
  windows <- st3_pair_windows(p$times)
  d <- st3_environment_design(p,sim$env,sets,windows,basis)
  y <- st3_collect(p,sets,windows,draws=8,seed=r)
  model <- st3_fit_process_gdum(d,y$y)
  contrast <- seq(-2,2,length.out=35)
  tested <- vector("list",length(contrast))
  start <- seq_len(nrow(sets))
  # Start at an actual training set near mean climate position. Hold both soil
  # endpoints fixed and recompute climate spread after a common site-level shift.
  reference_sites <- sets[which.min(abs(d$X[start,"position.climate"])),]
  world <- st3_subset(p,reference_sites,1:2)
  for (j in seq_along(contrast)) {
    e <- sim$env[reference_sites,1:2,,drop=FALSE]
    e[,2,"climate"] <- e[,1,"climate"] + contrast[j]*basis$native$scale["climate"]
    g <- st3_environment_design(world,e,matrix(seq_len(k),1),list(list(before=1L,after=2L)),basis)
    pr <- st3_predict_process_gdum(model,g)
    tested[[j]] <- data.frame(replicate=r,delta=contrast[j],process=colnames(pr$contributions),
      estimate=drop(pr$contributions),extrapolation=pr$extrapolation,row.names=NULL)
  }
  curves[[r]] <- do.call(rbind,tested)
}
save(do.call(rbind,curves),"environmental_curves")

sites <- list()
for (scenario in c("redistribution","opposing")) for (r in 1:12) {
  sim <- st3_simulate(scenario,seed=78000L+r)
  w <- st3_pair_windows(sim$panel$times)
  for (j in seq_along(w)) {
    p <- st3_plan(sim$panel,w[[j]]$before,w[[j]]$after)
    sites[[length(sites)+1L]] <- cbind(scenario,replicate=r,contrast=j,
                                      initial=p$alignment$initial,final=p$alignment$final,
                                      st3_site_reorganisation(p))
  }
}
save(do.call(rbind,sites),"site_reorganisation")
