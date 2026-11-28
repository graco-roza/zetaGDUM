Rcpp::sourceCpp(file.path(root,"scripts/rectangle_kernel.cpp"),cacheDir=file.path(tempdir(),"st4_cpp"))
st4_prepare <- function(a,ids,initial,final,taxonomic="harmonised") {
  p <- st4_panel(a,ids,c(initial,final),taxonomic)
  active <- colSums(p$x>0,dims=2L)>0
  p$x <- p$x[,,active,drop=FALSE]
  st3_plan(p,seq_along(initial),length(initial)+seq_along(final))
}
st4_indices <- function(sets,n,m) {
  I <- .st3_site_indices(sets,n,m)
  lapply(seq_len(nrow(I)),function(i) as.integer(I[i,]))
}
st4_eval <- function(plan,sets,draws=32L,seed=109,endpoint=FALSE) {
  idx <- st4_indices(sets,length(plan$sites),nrow(plan$alignment))
  masks <- if (endpoint) list(masks=list(plan$mask),exact=TRUE) else st3_allocations(plan,32L,draws,seed)
  z <- st4_evaluate_cpp(plan$before,plan$after,idx,masks$masks,
                        if(endpoint)matrix(0,16,0) else .st3_coalitions())
  colnames(z$before) <- colnames(z$after) <- metric_labels
  if (!endpoint) {
    dimnames(z$phi) <- list(NULL,metric_labels,st3_processes,NULL)
    colnames(z$support) <- metric_labels
    z$valid <- z$support==1
    z$exact <- masks$exact
    z$draws <- length(masks$masks)
    z$mean <- apply(z$phi,c(1,2,3),mean)
    dim(z$mean) <- dim(z$phi)[1:3]
    dimnames(z$mean) <- dimnames(z$phi)[1:3]
    stopifnot(z$efficiency_error<1e-10)
  }
  z
}
st4_inputs <- function(group,design,region_drop=NULL) {
  p <- panels()
  m <- p$membership
  ids <- m$id[m$assemblage==group&m$design==design&m$included]
  if (!is.null(region_drop)) ids <- setdiff(ids,region_drop)
  list(ids=ids,years=p$years[[design]],data=archive()[[group]])
}
