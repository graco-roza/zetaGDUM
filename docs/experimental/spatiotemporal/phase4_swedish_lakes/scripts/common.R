phase <- "docs/experimental/spatiotemporal"
root <- file.path(phase,"phase4_swedish_lakes")
source(file.path(phase,"phase3_load.R"))
st3_load(phase)
assemblages <- c("phytoplankton","littoral","profundal")
raw <- function(n) file.path(root,"data/raw",n)
processed <- function(n) file.path(root,"data/processed",n)
result <- function(n) file.path(root,"results",n)
fread <- function(...) data.table::fread(...,data.table=FALSE)
save_csv <- function(x,n) data.table::fwrite(x,result(paste0(n,".csv")),na="")
archive <- function() readRDS(processed("assemblages.rds"))
panels <- function() readRDS(processed("panels.rds"))
years <- 1992:2022
close <- function(x) x / rowSums(x)
sgn <- function(z) ifelse(z %in% c("zeta_presence","zeta_abundance","zeta_dominance"),-1,1)
st4_panel <- function(a, ids, years, taxonomic="harmonised") {
  z <- a[[taxonomic]]
  key <- as.vector(outer(ids,years,paste))
  j <- match(key,paste(a$samples$id,a$samples$year))
  if (anyNA(j) || any(a$samples$status[j]!="observed_nonempty")) stop("A matched panel contains unobserved or unsupported records.")
  x <- array(z[j,,drop=FALSE],c(length(ids),length(years),ncol(z)),list(ids,as.character(years),colnames(z)))
  st_panel(x,years)
}
metric_labels <- c("zeta_presence","zeta_abundance","zeta_dominance","incidence","turnover","nestedness",
"abundance","occupancy","numerical","temporal_presence","temporal_abundance","site_temporal_presence",
"site_temporal_abundance","site_temporal_dominance","temporal_dominance","dominance_dissimilarity")
proc_labels <- c("Expansion","Contraction","Redistribution","Relative shares")
jobs <- function(x,fn,cores=as.integer(Sys.getenv("ST4_WORKERS","5"))) {
  z <- parallel::mclapply(x,fn,mc.cores=if(.Platform$OS.type=="windows")1L else cores,mc.set.seed=FALSE)
  bad <- vapply(z,inherits,logical(1),"try-error")
  if (any(bad)) stop(paste(unlist(z[bad]),collapse="\n"))
  z
}
draw_sets <- function(n,k,budget) {
  if (choose(n,k)<=budget) return(.st_indices(n,k,budget))
  st3_sample_sets(n,k,budget)
}
input_panel <- function(group,design,data=archive(),definition=panels()) {
  m <- definition$membership
  list(data=data[[group]],ids=m$id[m$assemblage==group & m$design==design & m$included],
       years=definition$years[[design]])
}
