source("docs/final_gdum/design.R")
orders <- c(2L,4L,8L)
families <- c("contrasting","occurrence")
rows <- lapply(families,function(s) {
 do.call(rbind,lapply(orders,function(k) {
  t <- read_target(s,k)
  do.call(rbind,lapply(1:32,function(j) {
    a <- generate(110L,s,380881L+scenario_id(s)*10000L+31L*j)
    idx <- sampling(nrow(a$e),k,900L,3911L+k+j)
    y <- ugdm_component_response(a$x,idx)
    d <- ugdm_component_design(a$e,idx,"additive")
    z <- ugdm_component_design(t$e,t$idx,reference=d)
    f <- ugdm_fit_components(d,y)
    g <- ugdm_component_design(a$e,idx,"interaction",reference=d)
    v <- ugdm_component_design(t$e,t$idx,reference=g)
    base <- accuracy(ugdm_predict_components(f,z),t$expected)
    do.call(rbind,lapply(c(.01,.001,.0001,0),function(lambda) {
      h <- ugdm_fit_components(g,y,interaction_shrink=lambda)
      accuracy_result <- accuracy(ugdm_predict_components(h,v),t$expected)
      info <- ugdm_component_diagnostics(h,g,y)
      data.frame(scenario=s,order=k,replicate=j,penalty=lambda,response=responses,
        additive_RMSE=unname(base),interaction_RMSE=unname(accuracy_result),
        converged=h$converged,condition=max(info$stages$condition))
    }))
  }))
 }))
})
atomic_csv(do.call(rbind,rows),"interaction_penalties.csv")
# Noise-free positive control on the fitted link-scale basis. This is a
# numerical identification experiment, not a species-assembly simulation.
positive <- do.call(rbind,lapply(orders,function(k) {
  n <- 150L
  env <- environment(n,"contrasting",52013L)[1L]
  idx <- sampling(n,k,3200L,3031L+k)
  add <- ugdm_component_design(env,idx)
  inter <- ugdm_component_design(env,idx,"interaction",reference=add)
  eta <- -1.3+drop(inter$dir %*% c(.25,.20,.12,1.8,1.4,.9))
  mu <- plogis(eta)
  pi <- plogis(.6-drop(inter$dir%*%c(.12,.2,.3,.8,1.1,1.3)))
  y <- structure(cbind(total=mu,occupancy=mu*(1-pi),heterogeneity=mu*pi),
    sets=inter$idx,class=c("ugdm_component_response","matrix","array"))
  do.call(rbind,lapply(list(add,inter),function(d) {
    f <- ugdm_fit_components(d,y,ridge=1e-7,interaction_shrink=1e-7)
    p <- ugdm_predict_components(f,d)
    data.frame(order=k,model=d$model,response=responses,
      RMSE=unname(accuracy(p,y)),converged=f$converged)
  }))
}))
atomic_csv(positive,"interaction_positive_control.csv")
