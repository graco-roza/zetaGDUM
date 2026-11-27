source("docs/adaptive_refinement/ecology_common.R")
source("docs/final_gdum/design.R")
args <- commandArgs(TRUE)
meta <- read_manifest()[as.integer(args[2L]),,drop=FALSE]
a <- load_dd2(meta,args[1L])
fold <- study_blocks(a$dist,partitions=1L)[[1L]]
orders <- unique(c(2L,min(meta$mid,min(table(fold)))))
M <- 900L
profiles <- zeta_order_information(a$x,a$env,orders=sort(unique(c(2:min(12L,nrow(a$x)),orders))),
  response="all",combinations=400,pools=2,seed=9213)
atomic_csv(profiles,paste0("emp_profile_",meta$id,".csv"))
metrics <- do.call(rbind,lapply(orders,function(k) {
  summary <- ugdm_validate_components(a$x,a$env,order=k,fold=fold,combinations=M,seed=5501+k)
  cbind(meta[c("id","community")],summary)
}))
atomic_csv(metrics,paste0("emp_validation_",meta$id,".csv"))
out <- lapply(orders,function(k) {
  n <- nrow(a$x)
  idx <- sampling(n,k,M,8921+k)
  idx2 <- sampling(n,k,M,17321+k)
  fixed <- sampling(n,k,1800L,37122+k)
  y <- ugdm_component_response(a$x,idx)
  repeat_y <- ugdm_component_response(a$x,idx2)
  v <- ugdm_component_response(a$x,fixed)
  do.call(rbind,lapply(c("additive","interaction"),function(model) {
    d <- ugdm_component_design(a$env,idx,model)
    r <- ugdm_component_design(a$env,idx2,reference=d)
    z <- ugdm_component_design(a$env,fixed,reference=d)
    f <- ugdm_fit_components(d,y)
    alt <- ugdm_fit_components(r,repeat_y)
    q <- ugdm_component_diagnostics(f,d,y)
    p1 <- ugdm_predict_components(f,z)
    p2 <- ugdm_predict_components(alt,z)
    support <- ugdm_component_support(f,z)
    g <- importance(f,d,y,z,v,permutations=24L)
    err <- sqrt(colMeans((p1[responses]-p2[responses])^2)/2)
    detail <- q$interaction_residual$residual_fraction
    cbind(meta[c("id","community")],order=k,model,
      g$point, converged=f$converged, failed_coalitions=g$failures,
      condition=max(q$stages$condition),rank=min(q$stages$information_rank),
      parameters=max(q$stages$parameters),min_interaction_residual=if(length(detail))min(detail)else NA_real_,
      support_fraction=mean(support$supported), numerical_prediction_RMSE=unname(err[g$point$response]))
  }))
})
atomic_csv(do.call(rbind,out),paste0("emp_importance_",meta$id,".csv"))
