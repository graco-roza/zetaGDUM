source("docs/final_gdum/design.R")
args <- commandArgs(TRUE)
s <- args[1L]
orders <- as.integer(strsplit(args[2L],",")[[1L]])
i <- scenario_id(s)
dir <- reference_dir
dir.create(dir,recursive=TRUE,showWarnings=FALSE)
rows <- lapply(orders,function(k) {
  target <- truth_data(s,k)
  rows <- lapply(1:2,function(r) {
    a <- generate(2100L,s,6817L+10L*i+r)
    idx <- sampling(nrow(a$e),k,16000L,9029L+r+19L*k)
    y <- ugdm_component_response(a$x,idx)
    do.call(rbind,lapply(models,function(m) {
      d <- ugdm_component_design(a$e,idx,m)
      z <- ugdm_component_design(target$e,target$idx,reference=d)
      do.call(rbind,lapply(architectures,function(arch) {
        f <- ugdm_fit_components(d,y,arch)
        g <- importance(f,d,y,z,target$expected,"hierarchy")
        p <- ugdm_predict_components(f,z)
        err <- accuracy(p,target$expected)
        data.frame(scenario=s,order=k,reference=r,model=m,architecture=arch,
          g$hierarchy,oracle_RMSE=unname(err[g$hierarchy$response]),converged=f$converged,
          failed_coalitions=g$failures)
      }))
    }))
  })
  saveRDS(target,file.path(dir,paste0(s,"_",k,".rds")))
  do.call(rbind,rows)
})
atomic_csv(do.call(rbind,rows),paste0("reference_",s,".csv"))
