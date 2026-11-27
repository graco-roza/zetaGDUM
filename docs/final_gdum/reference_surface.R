source("docs/final_gdum/design.R")
rows <- lapply(mechanisms$name,function(s) {
  i <- scenario_id(s)
  a <- generate(2100L,s,6817L+10L*i+1L)
  do.call(rbind,lapply(c(2L,4L,8L),function(k) {
    t <- read_target(s,k)
    idx <- sampling(nrow(a$e),k,16000L,9030L+19L*k)
    y <- ugdm_component_response(a$x,idx)
    do.call(rbind,lapply(c("additive","interaction"),function(m) {
      d <- ugdm_component_design(a$e,idx,m)
      z <- ugdm_component_design(t$e,t$idx,reference=d)
      f <- ugdm_fit_components(d,y)
      p <- ugdm_predict_components(f,z)
      support <- ugdm_component_support(f,z)$supported
      oracle <- as.data.frame(unclass(t$expected))
      names(oracle) <- paste0("oracle_",names(oracle))
      data.frame(scenario=s,order=k,model=m,row=seq_len(nrow(p)),
        position=z$raw_position[,"phosphorus"],spread=z$raw_spread[,"phosphorus"],
        supported=as.logical(support),p,oracle,converged=f$converged)
    }))
  }))
})
atomic_csv(do.call(rbind,rows),"reference_surfaces.csv.gz")
