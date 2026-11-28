source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
reps <- as.integer(Sys.getenv("ST5_REPS","100"))
workers <- as.integer(Sys.getenv("ST5_WORKERS","4"))
assertions <- c("expansion","contraction","redistribution","numerical","opposing","stable","reversal","irregular")
tasks <- bind_rows(list(expand.grid(scenario=c("multiple","null"),rep=seq_len(reps)),
                        expand.grid(scenario=assertions,rep=seq_len(max(8,floor(reps/10)))),
                        data.frame(scenario="reference",rep=seq_len(reps))))
extract <- function(sim,records=sim$records,budget=32L,k=2,m=1,abundance=FALSE,draws=4L) {
  p <- if(abundance) st5_input(records,abundance=sim$abundance,env=sim$env) else
    st5_input(records,incidence=sim$abundance>0,env=sim$env)
  s <- st5_support(p,"change",k,m,budget=budget)
  r <- st5_responses(p,s,TRUE,draws=draws)
  z <- st5_extract(r,if(abundance) "abundance" else "incidence")
  list(p=p,z=z,s=s,valid=mean(complete.cases(z$phi)),
       error=r$efficiency_error,seed=sim$times[1])
}
stat <- function(obj) {
  z <- obj$z
  z$phi <- z$phi[complete.cases(z$phi),,drop=FALSE]
  c(net=mean(rowSums(z$phi)),colMeans(z$phi))
}
fit_one <- function(obj,method="temporal",grid=10^seq(-6,2)) {
  z <- obj$z
  p <- obj$p
  st5_validate(p,z,method,grid=grid,blocks=3L)
}
scenario <- function(task) {
  n <- as.character(task$scenario)
  j <- task$rep
  seed <- 8000+j
  key <- data.frame(scenario=n,replicate=j)
  if(n=="reference") {
    # Independent replicate seeds and four times the principal combination budget.
    seed <- 223000+j
    intact <- st5_simulate(seed)
    a <- extract(intact,budget=128L)
    fit <- fit_one(a,grid=c(.0001,.01,1))
    rank <- st5_importance(fit)
    rank <- rank[rank$fold=="pooled",]
    drivers <- lapply(1:2,function(v) {
      effect <- rep(1,2); effect[v] <- 0
      off <- extract(st5_simulate(seed,effect=effect),budget=128L)
      difference <- apply((a$z$phi-off$z$phi)^2,2,mean,na.rm=TRUE)^.5
      cbind(key,driver=paste0("E",v),response=names(difference),contrast_rms=unname(difference))
    })
    return(list(importance=cbind(key,rank),reference=bind_rows(drivers)))
  }
  sim <- st5_simulate(seed,scenario=n)
  use_abundance <- n%in%c("numerical","opposing")
  obj <- extract(sim,abundance=use_abundance)
  sx <- stat(obj)
  temporal_order <- match(sim$records$time,sim$times)
  # Analytical spatial mean incidence zeta, by survey, at order two.
  zeta <- vapply(seq_along(sim$times),function(t) {
    count <- colSums(sim$abundance[temporal_order==t,,drop=FALSE]>0)
    sum(choose(count,2))/choose(sum(temporal_order==t),2)
  },numeric(1))
  local <- st5_responses(obj$p,st5_support(obj$p,"local",budget=Inf),FALSE)
  metric <- if(use_abundance) "site_temporal_abundance" else "site_temporal_presence"
  mechanism <- data.frame(key,net=sx["net"],mean_local=mean(local$after[,metric],na.rm=TRUE),
    valid=obj$valid,shapley_error=obj$error,zeta_endpoint_change=tail(zeta,1)-zeta[1],
    max_zeta_increase=max(diff(zeta)),min_zeta_change=min(diff(zeta)),check.names=FALSE)
  for(v in setdiff(names(sx),"net")) mechanism[[v]] <- sx[v]
  if(!n%in%c("multiple","null")) {
    # Test process activity directly for isolated generating mechanisms.
    expected <- intersect(n,colnames(obj$z$phi))
    off <- setdiff(colnames(obj$z$phi),expected)
    inactive <- if(length(expected)) max(abs(obj$z$phi[,off,drop=FALSE]),na.rm=TRUE) else NA_real_
    mechanism$max_inactive <- if(n=="numerical") inactive else if(length(expected)) inactive else NA_real_
    mechanism$nonzero_expected <- if(length(expected)) mean(abs(obj$z$phi[,expected])>1e-10,na.rm=TRUE) else NA_real_
    return(list(mechanisms=mechanism))
  }
  validation <- fit_one(obj)
  rank <- st5_importance(validation)
  rank <- rank[rank$fold=="pooled" | rank$response=="joint",]
  perf <- validation$scores
  if(n=="multiple" && j<=25L) {
    gap <- lapply(c("random","network","period","environment","cluster"),function(scheme) {
      lapply(c(.10,.25,.50),function(rate) {
        records <- st5_missing(sim,rate,scheme,seed+j+114)
        # Include one all-observation estimate on each recorded date, even when matching fails.
        h <- extract(sim,records)
        sp <- st5_responses(h$p,st5_support(h$p,"spatial",budget=96L))
        full <- st5_responses(obj$p,st5_support(obj$p,"spatial",budget=96L))
        ts <- cbind(head(sim$times,-1),tail(sim$times,-1))
        all <- st5_support_change(sp,list(rows=h$z$rows,before=cbind(incidence=h$z$baseline),
                    after=cbind(incidence=h$z$y),phi=NULL),"incidence",ts)
        full_all <- st5_support_change(full,list(rows=obj$z$rows,before=cbind(incidence=obj$z$baseline),
                    after=cbind(incidence=obj$z$y),phi=NULL),"incidence",ts)
        valid <- is.finite(all$matched)&is.finite(full_all$matched)
        availability <- st5_eligibility(h$p,h$s)
        v <- fit_one(h,grid=c(.0001,.01,1))
        rank <- st5_importance(v)
        pool <- rank[rank$fold=="pooled" & rank$response=="joint",]
        metrics <- v$scores[v$scores$model=="full",]
        data.frame(key,scheme,rate,observed=mean(records$available),
          used_fraction=availability$eligible_observation_fraction,
          relative_eligible_pairs=availability$eligible_rows/st5_eligibility(obj$p,obj$s)$eligible_rows,
          missing_periods=sum(!sim$times%in%h$p$records$time),
          matched_contrasts=sum(valid),mean_all_minus_matched=mean(all$discrepancy[valid]),
          rmse_all_minus_matched=sqrt(mean(all$discrepancy[valid]^2)),
          bias_matched=mean(all$matched[valid]-full_all$matched[valid]),
          bias_all=mean(all$all[valid]-full_all$all[valid]),
          endpoint_rmse=mean(metrics$rmse),joint_rmse=sqrt(mean(metrics$joint)),
          best_predictor=if(nrow(pool)&&any(is.finite(pool$importance)))
            pool$feature[which.max(pool$importance)] else NA_character_)
      })
    })
    gap <- bind_rows(unlist(gap,recursive=FALSE))
  } else gap <- NULL
  order <- if(n=="multiple" && j<=12) bind_rows(lapply(2:4,function(k)
    bind_rows(lapply(1:3,function(m) {
      h <- extract(sim,budget=32L,k=k,m=m)
      data.frame(key,k,m,st5_eligibility(h$p,h$s)[c("combination_rows","eligible_rows")],
                 responses=length(h$z$y),valid=h$valid,zero_sharing=mean(h$z$y==1,na.rm=TRUE),
                 rms_process=sqrt(mean(h$z$phi^2,na.rm=TRUE)))
    })))) else NULL
  list(importance=cbind(key,rank),performance=cbind(key,perf),mechanisms=mechanism,missing=gap,orders=order)
}
res <- if(workers>1L) parallel::mclapply(seq_len(nrow(tasks)),function(i) scenario(tasks[i,]),
                  mc.cores=workers,mc.preschedule=FALSE) else lapply(seq_len(nrow(tasks)),function(i) scenario(tasks[i,]))
if(any(vapply(res,inherits,logical(1),"try-error"))) stop(res[[which(vapply(res,inherits,logical(1),"try-error"))[1]]])
for(n in unique(unlist(lapply(res,names)))) write_result(bind_rows(lapply(res,`[[`,n)),paste0("simulation_",n))
