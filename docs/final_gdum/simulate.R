source("docs/final_gdum/design.R")
args <- commandArgs(TRUE)
s <- args[1L]
n <- as.integer(args[2L])
reps <- as.integer(args[3L])
orders <- as.integer(strsplit(args[4L],",")[[1L]])
M <- as.integer(args[5L])
phase <- args[6L]
id <- scenario_id(s)
workers <- as.integer(Sys.getenv("MS_GDUM_WORKERS","1"))
targets <- lapply(orders, function(k) read_target(s,k))
one <- function(i) {
  seed <- 380881L+id*10000L+31L*i+if(phase=="pilot") 400000L else 0L
  a <- generate(n,s,seed)
  do.call(rbind,lapply(seq_along(orders),function(j) {
    k <- orders[j]
    t <- targets[[j]]
    idx <- sampling(n,k,M,seed+1187L+k)
    y <- ugdm_component_response(a$x,idx)
    idx2 <- sampling(n,k,M,seed+5071L+k)
    repeat_y <- ugdm_component_response(a$x,idx2)
    do.call(rbind,lapply(models,function(m) {
      d <- ugdm_component_design(a$e,idx,m)
      z <- ugdm_component_design(t$e,t$idx,reference=d)
      repeat_d <- ugdm_component_design(a$e,idx2,reference=d)
      do.call(rbind,lapply(architectures,function(arch) {
        start <- proc.time()[["elapsed"]]
        f <- ugdm_fit_components(d,y,arch)
        ft <- proc.time()[["elapsed"]]-start
        p <- ugdm_predict_components(f,z)
        scores <- ugdm_assess_components(f,z,t$observed,y)
        quality <- ugdm_component_diagnostics(f,d,y)
        g <- importance(f,d,y,z,t$expected)
        repeated <- ugdm_fit_components(repeat_d,repeat_y,arch)
        p2 <- ugdm_predict_components(repeated,z)
        numerics <- vapply(responses,function(v) sqrt(mean((p[[v]]-p2[[v]])^2)/2),numeric(1))
        err <- accuracy(p,t$expected)
        b <- merge(scores$metrics,g$point,by.x="component",by.y="response",suffixes=c("_observed","_oracle"))
        v <- b$component
        res <- quality$interaction_residual$residual_fraction
        cbind(data.frame(scenario=s,n,replicate=i,seed,order=k,model=m,architecture=arch,pool_sets=M),
          b,oracle_RMSE=unname(err[v]),numerical_prediction_RMSE=unname(numerics[v]),closure_RMSE=scores$closure_RMSE,
          converged=f$converged,repeat_converged=repeated$converged,
          min_rank=min(quality$stages$information_rank),parameters=max(quality$stages$parameters),
          max_condition=max(quality$stages$condition),active_total=quality$stages$active[1L],
          min_interaction_residual=if(length(res)) min(res) else NA_real_,
          failed_coalitions=g$failures,fit_seconds=ft,total_seconds=proc.time()[["elapsed"]]-start)
      }))
    }))
  }))
}
run <- if (workers>1L && .Platform$OS.type!="windows") {
  parallel::mclapply(seq_len(reps), one, mc.cores=workers, mc.set.seed=FALSE)
} else lapply(seq_len(reps),one)
stopifnot(all(vapply(run,is.data.frame,logical(1))))
atomic_csv(do.call(rbind,run),paste0(s,"_",phase,".csv.gz"))
