source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
source(file.path(phase,"scripts/archives.R"))
archive <- st5_lake_archive(file.path(base,"phase4_swedish_lakes"))
out <- list()
for(dataset in names(archive$inputs)) {
  q <- archive$inputs[[dataset]]$p
  p <- st5_input(q$records,incidence=q$x,env=q$env)
  block <- archive$blocks[p$sites]
  support <- function(sites,seed) {
    shape <- st5_support(p,"change",lag=1,budget=2L)
    set.seed(seed)
    sampled <- lapply(shape,function(g) {
      g$focal <- intersect(g$region,sites)
      if(length(g$focal)<2L) return(NULL)
      idx <- st5_sets(length(g$focal),2L,64L)
      g$sets <- matrix(g$focal[idx],ncol=2)
      g$eligible <- choose(length(g$focal),2)
      g
    })
    structure(Filter(Negate(is.null),sampled),type="change",k=2L,m=1L)
  }
  relabel <- function(g) {
    ans <- lapply(g,function(x) {x$region <- x$focal; x})
    attributes(ans) <- attributes(g)
    ans
  }
  response <- function(x) st5_extract(st5_responses(p,x,TRUE,draws=4L),"incidence")
  for(b in sort(unique(block))) {
    tr <- support(which(block!=b),100+b)
    te <- support(which(block==b),200+b)
    a <- response(tr)
    b0 <- response(te)
    ar <- response(relabel(tr))
    br <- response(relabel(te))
    key <- data.frame(dataset,block=b)
    same_train <- complete.cases(a$phi)&complete.cases(ar$phi)
    same_test <- complete.cases(b0$phi)&complete.cases(br$phi)
    a <- st5_subset(a,same_train); ar <- st5_subset(ar,same_train)
    b0 <- st5_subset(b0,same_test); br <- st5_subset(br,same_test)
    i <- seq_along(a$y); j <- seq_along(b0$y)
    if(length(i)<20L || length(j)<10L) next
    outputs <- lapply(c("conditional","independent_regions"),function(target) {
      z <- if(target=="conditional") a else ar
      x <- if(target=="conditional") b0 else br
      f <- .st5_fold_data(p,z,list(train=i,test=i),colnames(p$env),"coherent")
      design <- st5_design(p,x,f$ref)
      gg <- .st5_valid(x,design,"coherent")
      tune <- st5_tune(p,z,st5_folds(z,"temporal",3),colnames(p$env),"coherent")
      f$mode <- "coherent"; f$penalty <- tune; f$equal_rows <- f$monotonic <- FALSE
      pred <- st5_predict(.st5_refit(f),x,design)
      data.frame(key,target,train=length(i),test=length(gg),as.list(.st5_loss(x,pred,"coherent",gg)),
        mean_absolute_label_change=mean(abs(b0$phi[gg,]-br$phi[gg,])))
    })
    out[[paste(dataset,b)]] <- bind_rows(outputs)
  }
}
write_result(bind_rows(out),"empirical_spatial_estimands")
