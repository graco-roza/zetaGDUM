source("docs/order_information/design.R")
name <- commandArgs(TRUE)[1L]
workers <- as.integer(Sys.getenv("ZETA_ORDER_WORKERS", "1"))
sources <- list.files(output, paste0("^",name,"_.*_(train|test|principal|validation)_(total|relative)_calibration[.]csv[.]gz$"),
                      full.names=TRUE)
cells <- do.call(rbind,lapply(sources,function(f) {
  a <- read.csv(gzfile(f))
  unique(a[c("name","n","split","replicate","seed","order")])
}))
cells <- unique(cells)
groups <- split(cells, interaction(cells$n,cells$split,cells$replicate,drop=TRUE))
compute <- function(a) {
  c <- community(a$n[1],a$name[1],a$seed[1])
  d <- suppressWarnings(zeta_order_information(c$x,c$env,orders=sort(a$order),response="all",
    method="sample", combinations=320L, pools=2L, threshold_limit=128L,
    seed=18L+a$seed[1]))
  p <- d$mean_zeta[d$order==2]
  names(p) <- d$response[d$order==2]
  d$zeta_remaining <- unname(d$mean_zeta/p[d$response])
  data.frame(a[1,c("name","n","split","replicate","seed")],d)
}
rows <- if (workers>1L && .Platform$OS.type!="windows") {
  parallel::mclapply(groups,compute,mc.cores=workers,mc.set.seed=FALSE)
} else lapply(groups,compute)
stopifnot(all(vapply(rows,is.data.frame,logical(1))))
write_result(do.call(rbind,rows),file.path(output,paste0("profile_",name,".csv.gz")))
