source("docs/final_gdum/design.R")
cells <- expand.grid(n=c(70L,110L,250L),order=c(2L,4L,8L),replicate=1:4)
a <- do.call(rbind,lapply(seq_len(nrow(cells)),function(i) {
  cell <- cells[i,]
  draw <- generate(cell$n,"contrasting",71238L+i)
  t <- read_target("contrasting",cell$order)
  do.call(rbind,lapply(models,function(m) {
    prep_start <- proc.time()[["elapsed"]]
    idx <- sampling(cell$n,cell$order,900L,7314L+i)
    train <- ugdm_component_design(draw$e,idx,m)
    y <- ugdm_component_response(draw$x,idx)
    test <- ugdm_component_design(t$e,t$idx,reference=train)
    prep <- proc.time()[["elapsed"]]-prep_start
    do.call(rbind,lapply(architectures,function(arch) {
      start <- proc.time()[["elapsed"]]
      fit <- ugdm_fit_components(train,y,arch)
      fitting <- proc.time()[["elapsed"]]-start
      start <- proc.time()[["elapsed"]]
      g <- importance(fit,train,y,test,t$expected,resolution="hierarchy")
      inference <- proc.time()[["elapsed"]]-start
      data.frame(cell,model=m,architecture=arch,preparation_seconds=prep,
        fit_seconds=fitting,allocation_seconds=inference,elapsed=prep+fitting+inference,
        converged=fit$converged,failed_coalitions=g$failures)
    }))
  }))
}))
atomic_csv(a,"benchmark.csv")
