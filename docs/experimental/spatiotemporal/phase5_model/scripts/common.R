base <- "docs/experimental/spatiotemporal"
phase <- file.path(base,"phase5_model")
source(file.path(phase,"load.R"))
cache <- file.path(phase,".cpp-cache")
dir.create(cache,showWarnings=FALSE)
st5_load(base,cache=cache)
write_result <- function(x,name) write.csv(x,file.path(phase,"results",paste0(name,".csv")),row.names=FALSE)
read_result <- function(name) read.csv(file.path(phase,"results",paste0(name,".csv")),check.names=FALSE)
bind_rows <- function(x) data.table::rbindlist(x,fill=TRUE,use.names=TRUE)
spatial_blocks <- function(site,x,y,n=3L) {
  z <- ord <- order(x,y)
  z[ord] <- cut(seq_along(ord),n,labels=FALSE)
  setNames(z,site)
}
