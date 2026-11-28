source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
source(file.path(root,"scripts/engine.R"))
data <- archive(); def <- panels(); geo <- fread(processed("lake_geography.csv"))
p <- data.table::fread(result("process_replicates.csv"))
rl <- c("metric","process")
safe_mean <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
s <- p[,.(before=safe_mean(before),after=safe_mean(after),delta=safe_mean(delta),
          signed_mean=safe_mean(mean_raw),negative=safe_mean(negative),positive=safe_mean(positive),
          valid_fraction=mean(defined/sets),endpoint_fraction=mean(defined_endpoints/sets),
          undefined_intermediates=sum(undefined_intermediate),zero_sharing=mean(zero_before)),
       by=.(assemblage,spatial_order,temporal_order,metric,process,replicate)]
main <- s[,.(before=safe_mean(before),after=safe_mean(after),delta=safe_mean(delta),signed_mean=safe_mean(signed_mean),
             negative=safe_mean(negative),positive=safe_mean(positive),
             combination_mc_se=if(.N>1L)sd(signed_mean)/sqrt(.N) else 0,
             combination_replicates=.N,valid_fraction=mean(valid_fraction),
             endpoint_fraction=mean(endpoint_fraction),undefined_intermediates=sum(undefined_intermediates),
             zero_sharing=mean(zero_sharing)),by=.(assemblage,spatial_order,temporal_order,metric,process)]
main[,gross:=negative+positive]
main[,negative_share:=100*.st_ratio(negative,sum(negative)),by=.(assemblage,spatial_order,temporal_order,metric)]
main[,positive_share:=100*.st_ratio(positive,sum(positive)),by=.(assemblage,spatial_order,temporal_order,metric)]
save_csv(main,"process_summary")
z <- data.table::fread(result("zeta_surfaces.csv"))
z <- z[taxonomy=="harmonised" & window=="consecutive" & (is.na(threshold)|threshold==.8)]
z <- z[,.(early_zeta=zeta[period=="early"],late_zeta=zeta[period=="late"],
          early_effective=effective_species[period=="early"],late_effective=effective_species[period=="late"],
          spatial_retention=spatial_retention[period=="late"],temporal_retention=temporal_retention[period=="late"]),
       by=.(assemblage,spatial_order,temporal_order,response)]
# Response-specific validity and zero-sharing diagnostics. The latter uses
# the same evaluation windows and a separate uniform endpoint sample.
d <- main[grepl("^zeta_",metric),
          .(valid_fraction=valid_fraction[1],combination_mc_se=max(combination_mc_se,na.rm=TRUE)),
          by=.(assemblage,spatial_order,temporal_order,response=sub("^zeta_","",metric))]
d$combination_mc_se[!is.finite(d$combination_mc_se)] <- NA_real_
tmp <- main[process=="expansion" & metric%in%c("temporal_presence","temporal_abundance","temporal_dominance"),
            .(assemblage,spatial_order,temporal_order,response=sub("^temporal_","",metric),
              temporal_support=valid_fraction,undefined_intermediates)]
zero <- do.call(rbind,jobs(seq_along(assemblages),function(i) {
    g <- assemblages[i]; a <- input_panel(g,"A",data,def); n <- length(a$ids)
    ans <- list()
    for(m in 1:5) {
      set.seed(8200L+i*100L+m)
      sets <- lapply(1:6,function(k)draw_sets(n,k,240L))
      sizes <- vapply(sets,nrow,integer(1))
      ids <- do.call(c,lapply(sets,st4_indices,n=n,m=m))
      membership <- rep(1:6,sizes)
      total <- matrix(0,6,3)
      for(j in seq_len(6-m)) {
        ix <- j+seq_len(m)-1L
        q <- st4_prepare(a$data,a$ids,a$years[ix],a$years[ix+5L])
        z <- st4_evaluate_cpp(q$before,q$after,ids,list(q$mask),matrix(0,16,0))
        total <- total+rowsum((z$before[,1:3]==0)*1.0,membership)/(sizes*(6-m))
      }
      ans[[m]] <- data.frame(assemblage=g,spatial_order=rep(1:6,3),temporal_order=m,
                 response=rep(c("presence","abundance","dominance"),each=6),zero_sharing=as.vector(total))
    }
    do.call(rbind,ans)
},cores=3L))
cols <- c("assemblage","spatial_order","temporal_order","response")
z <- Reduce(function(x,y)merge(x,y,by=cols),list(z,d,tmp,zero))
z[,informative_shared:=temporal_support>=.8 & pmin(early_effective,late_effective)>=2 & zero_sharing<=.5]
z[,qualified_spatial:=spatial_order>=2 & valid_fraction>=.8]
save_csv(z,"informative_orders")
comb <- data.table::fread(result("process_combinations.csv.gz"))
spatial <- setNames(geo$spatial_block,geo$EU_ID)
comb[,in_block:=lapply(strsplit(sites,"|",fixed=TRUE),function(ids)unique(spatial[ids]))]
ld <- comb[k==2 & metric%in%c("abundance","occupancy","numerical","incidence","dominance_dissimilarity")]
ld[,total:=expansion+contraction+redistribution+numerical]
components <- c(st3_processes,"total")
jack <- lapply(1:6,function(b) {
  x <- ld[!vapply(in_block,function(w)b%in%w,logical(1))]
  y <- x[,lapply(.SD,safe_mean),by=.(assemblage,k,m,metric,offset),.SDcols=components]
  cbind(y[,lapply(.SD,safe_mean),by=.(assemblage,k,m,metric),.SDcols=components],block=b)
})
jack <- data.table::rbindlist(jack)
u <- data.table::melt(jack,id.vars=c("assemblage","k","m","metric","block"),
                    measure.vars=components,variable.name="process",value.name="delete_estimate")
bs <- u[,.(conditional_block_se=if(all(is.finite(delete_estimate)))sqrt(5*mean((delete_estimate-mean(delete_estimate))^2)) else NA_real_,
           spatial_blocks=.N),by=.(assemblage,spatial_order=k,temporal_order=m,metric,process)]
save_csv(bs,"long_spatial_uncertainty")
cross <- merge(main[metric%in%c("occupancy","numerical")],bs,
                  by=c("assemblage","spatial_order","temporal_order","metric","process"),all.x=TRUE)
cross[,negative_share_across_eight:=100*.st_ratio(negative,sum(negative)),by=.(assemblage,spatial_order,temporal_order)]
cross[,positive_share_across_eight:=100*.st_ratio(positive,sum(positive)),by=.(assemblage,spatial_order,temporal_order)]
save_csv(cross,"crossed_summary")

# Nested budgets with one shared allocation stream across spatial orders.
sensitivity <- function(i) {
  g <- assemblages[i]; a <- input_panel(g,"A",data,def); n <- length(a$ids)
  p <- st4_prepare(a$data,a$ids,a$years[1],a$years[6])
  set.seed(1800+i)
  all <- lapply(c(2L,4L,6L),function(k)draw_sets(n,k,640))
  ids <- do.call(c,lapply(all,st4_indices,n=n,m=1L))
  alloc <- st3_allocations(p,32,32,4900+i)
  z <- st4_evaluate_cpp(p$before,p$after,ids,alloc$masks,.st3_coalitions())
  metric <- match(c("abundance","occupancy","numerical"),metric_labels)
  cursor <- 0L; budgets <- list()
  for(k in seq_along(all)) {
    M <- nrow(all[[k]])
    block <- cursor+seq_len(M);cursor<-max(block)
    for(b in c(40L,80L,160L,320L,640L)) {
      v <- apply(z$phi[head(block,b),metric,,,drop=FALSE],c(2,3),mean)
      budgets[[length(budgets)+1L]] <- data.frame(assemblage=g,k=c(2,4,6)[k],budget=min(M,b),
        metric=rep(metric_labels[metric],times=4),process=rep(st3_processes,each=3),mean=as.vector(v))
    }
  }
  # 512 draws from a single stream, compared by identical prefixes and support.
  allocation <- list()
  for(m in c(1L,3L)) {
    p <- st4_prepare(a$data,a$ids,a$years[seq_len(m)],a$years[5+seq_len(m)])
    set.seed(5800+i+m)
    s <- draw_sets(n,if(m==1)2 else 4,40)
    z <- st4_eval(p,s,512L,6400+i+m)
    ids <- match(c("abundance","occupancy","numerical","incidence","dominance_dissimilarity"),metric_labels)
    for(j in ids) for(b in c(32L,128L,512L)) {
      supported <- z$valid[,j]
      x <- z$phi[supported,j,,seq_len(min(b,z$draws)),drop=FALSE]
      draw_means <- apply(x,c(3,4),mean)
      allocation[[length(allocation)+1L]] <- data.frame(assemblage=g,k=if(m==1)2 else 4,m,
        metric=metric_labels[j],process=st3_processes,budget=b,used=ncol(draw_means),
        signed_mean=rowMeans(draw_means),allocation_sd=apply(draw_means,1,sd),
        allocation_mc_se=if(z$exact)0 else apply(draw_means,1,sd)/sqrt(ncol(draw_means)),
        supported_combinations=sum(supported),sets=nrow(s),row.names=NULL)
    }
  }
  list(budget=do.call(rbind,budgets),allocation=do.call(rbind,allocation))
}
diag <- jobs(seq_along(assemblages),sensitivity,cores=3)
save_csv(do.call(rbind,lapply(diag,`[[`,"budget")),"combination_budget")
save_csv(do.call(rbind,lapply(diag,`[[`,"allocation")),"allocation_sensitivity")

# Reclassify only in separately labelled panel and taxonomy sensitivities.
scenario <- expand.grid(g=assemblages,type=c("genus","drop_block_3","original_dates"),stringsAsFactors=FALSE)
reclass <- function(i) {
  g<-scenario$g[i]; type<-scenario$type[i]
  a<-input_panel(g,if(type=="original_dates")"A_original" else "A",data,def)
  if(type=="drop_block_3")a$ids<-a$ids[spatial[a$ids]!=3]
  if(length(a$ids)<2L)return(NULL)
  set.seed(1800+match(g,assemblages))
  sets<-draw_sets(length(a$ids),2,160)
  do.call(rbind,lapply(1:5,function(j){
    p<-st4_prepare(a$data,a$ids,a$years[j],a$years[5+j],if(type=="genus")"genus" else "harmonised")
    z<-st4_eval(p,sets,32L,7000+j)
    do.call(rbind,lapply(c("abundance","incidence","occupancy","numerical","dominance_dissimilarity"),function(k)
      data.frame(assemblage=g,type,initial=a$years[j],final=a$years[5+j],lakes=length(a$ids),metric=k,
                 process=st3_processes,mean=colMeans(z$mean[,k,],na.rm=TRUE))))
  }))
}
save_csv(do.call(rbind,jobs(seq_len(nrow(scenario)),reclass)),"panel_taxonomy_sensitivity")
annual <- data.table::fread(result("B_annual.csv"))
deletes <- data.table::fread(result("B_deletion.csv"))
moving <- function(n,width=4L) {
  start<-sample.int(n-width+1L,ceiling(n/width),replace=TRUE)
  head(as.vector(outer(0:(width-1L),start,"+")),n)
}
set.seed(81401)
annual_u <- lapply(split(deletes,deletes[,.(assemblage,k,metric,process)],drop=TRUE),function(q) {
  actual<-annual[assemblage==q$assemblage[1]&k==q$k[1]&metric==q$metric[1]&process==q$process[1]]
  actual<-actual[order(year)]
  del<-data.table::dcast(q,year~block,value.var="mean")
  del<-as.matrix(del[match(actual$year,del$year),-1])
  if(nrow(actual)<12L || any(!is.finite(c(actual$mean,del)))) return(NULL)
  jack<-function(x) {v<-colMeans(x);sqrt((length(v)-1)*mean((v-mean(v))^2))}
  b<-replicate(400L,{
    j<-moving(nrow(actual))
    c(mean(actual$mean[j]),jack(del[j,,drop=FALSE]))
  })
  data.frame(actual[1,.(assemblage,k,metric,process)],comparisons=nrow(actual),mean=mean(actual$mean),
    spatial_se=jack(del),temporal_se=sd(b[1,]),crossed_se=sqrt(var(b[1,])+mean(b[2,]^2)),
    temporal_block_width=4L,interpretation="conditional annual mean; crossed combines temporal variance and mean conditional spatial variance")
})
save_csv(do.call(rbind,annual_u),"annual_block_uncertainty")

# Recompute local dominance and all its counterfactuals at three thresholds,
# retaining the same spatial sample and allocation stream.
dominance <- lapply(seq_along(assemblages),function(i){
  g <- assemblages[i]; a <- input_panel(g,"A",data,def); n <- length(a$ids)
  set.seed(6660L+i)
  sets <- draw_sets(n,2,120)
  idx <- st4_indices(sets,n,1L)
  do.call(rbind,lapply(1:5,function(j){
    p <- st4_prepare(a$data,a$ids,a$years[j],a$years[5+j])
    mask <- st3_allocations(p,32,32,7120L+i*10L+j)
    do.call(rbind,lapply(c(.5,.8,.9),function(threshold){
      z <- st4_evaluate_cpp(p$before,p$after,idx,mask$masks,.st3_coalitions(),threshold)
      metrics <- match(c("dominance_dissimilarity","site_temporal_dominance"),metric_labels)
      do.call(rbind,lapply(metrics,function(metric) {
        x <- apply(z$phi[,metric,,,drop=FALSE],c(1,3),mean)
        k <- complete.cases(x)
        data.frame(assemblage=g,offset=j,threshold,metric=metric_labels[metric],
                   before=mean(z$before[k,metric]),after=mean(z$after[k,metric]),
                   process=st3_processes,signed_mean=colMeans(x[k,,drop=FALSE]),support=mean(k))
      }))
    }))
  }))
})
save_csv(do.call(rbind,dominance),"dominance_threshold_sensitivity")
