source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
source(file.path(phase,"R/metrics_abundance.R"))
data <- archive(); definition <- panels()
surface <- species <- list()
for(g in assemblages) for(period in c("early","late")) for(tax in c("harmonised","genus")) {
  x <- input_panel(g,"A",data,definition)
  J <- if(period=="early")x$years[1:5] else x$years[6:10]
  p <- st4_panel(x$data,x$ids,J,tax)
  keep <- colSums(p$x>0,dims=2L)>0
  p$x <- p$x[,,keep,drop=FALSE]
  for(mode in c("all","consecutive")) {
    z <- st2_zeta(p,6,5,mode=mode)
    main <- z$surface
    main$threshold <- ifelse(main$response=="dominance",.8,NA_real_)
    surface[[length(surface)+1L]] <- cbind(assemblage=g,period,taxonomy=tax,main)
    index <- main[main$response=="presence",c("spatial_order","temporal_order")]
    if(tax=="harmonised") for(type in names(z$weights)) {
      i <- rep(seq_len(nrow(index)),each=ncol(z$weights[[type]]))
      v <- as.vector(t(z$weights[[type]]))
      species[[length(species)+1L]] <- data.frame(assemblage=g,period,mode,response=type,
          index[i,,drop=FALSE],taxon=rep(colnames(z$weights[[type]]),nrow(index)),
          weight=v,contribution=as.vector(t(z$contributions[[type]])))
    }
    # Exact threshold sensitivity on both explicitly defined temporal schedules.
    if(tax=="harmonised") for(threshold in c(.5,.9)) {
      d <- dim(p$x)
      dom <- zetaGDUM::dominance_matrix(matrix(p$x,nrow=d[1]*d[2]),threshold)
      pp <- p
      pp$x[] <- dom
      zz <- st2_zeta(pp,6,5,mode=mode,empty="undefined")
      alt <- zz$surface[zz$surface$response=="presence",]
      alt$response <- "dominance"
      alt$threshold <- threshold
      surface[[length(surface)+1L]] <- cbind(assemblage=g,period,taxonomy=tax,alt)
      ww <- zz$weights$presence
      ii <- rep(seq_len(nrow(index)),each=ncol(ww))
      species[[length(species)+1L]] <- data.frame(assemblage=g,period,mode,
        response=paste0("dominance_",threshold),index[ii,,drop=FALSE],
        taxon=rep(colnames(ww),nrow(index)),weight=as.vector(t(ww)),
        contribution=as.vector(t(zz$contributions$presence)))
    }
  }
}
save_csv(do.call(rbind,surface),"zeta_surfaces")
data.table::fwrite(do.call(rbind,species),result("zeta_species_weights.csv.gz"))

# Native, analytically exact snapshot trajectories on the independent annual
# panel. These coexist with the common-combination dissimilarity accounting.
annual <- list()
for(g in assemblages) {
  x <- input_panel(g,"B",data,definition)
  p <- st4_panel(x$data,x$ids,x$years)
  for(j in seq_along(p$times)) for(type in c("presence","abundance","dominance")) {
    a <- .st_slice(p,j)
    q <- zetaGDUM::zeta_curve(a,6,type,p=.8)
    annual[[length(annual)+1L]] <- cbind(assemblage=g,year=p$times[j],response=type,q)
  }
}
save_csv(do.call(rbind,annual),"annual_exact_zeta")
