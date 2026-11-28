source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
source(file.path(root,"scripts/engine.R"))
data <- archive(); definition <- panels()
save_rows <- function(z,group,design,k,m,offset,initial,final,replicate,budget) {
  selected <- if(k==1) which(grepl("^zeta_|temporal",metric_labels)) else seq_along(metric_labels)
  do.call(rbind,lapply(selected,function(j) {
    q <- z$mean[,j,,drop=FALSE]
    phi <- matrix(q,nrow=length(replicate))
    take <- z$valid[,j]
    a <- if(any(take)) mean(z$before[take,j]) else NA_real_
    b <- if(any(take)) mean(z$after[take,j]) else NA_real_
    raw <- phi[take,,drop=FALSE]
    y <- raw*sgn(metric_labels[j])
    H <- colMeans(pmax(-y,0))
    D <- colMeans(pmax(y,0))
    if(any(take)) stopifnot(abs(sum(colMeans(raw))-(b-a))<1e-9)
    is_distance <- startsWith(metric_labels[j],"temporal_") || startsWith(metric_labels[j],"site_temporal")
    data.frame(assemblage=group,design,spatial_order=k,temporal_order=m,offset,
      initial=paste(initial,collapse=":"),final=paste(final,collapse=":"),
      duration=diff(range(initial)),lag=mean(final-initial),replicate=replicate[1],
      metric=metric_labels[j],process=st3_processes,sets=length(take),budget,
      defined=sum(take),defined_endpoints=sum(is.finite(z$before[,j])&is.finite(z$after[,j])),
      mean_allocation_support=mean(z$support[,j]),undefined_intermediate=sum(!take & is.finite(z$before[,j]) & is.finite(z$after[,j])),
      before=a,after=b,delta=b-a,mean_raw=colMeans(raw),negative=H,positive=D,gross=H+D,
      negative_share=100*.st_ratio(H,sum(H)),positive_share=100*.st_ratio(D,sum(D)),
      sign_convention=if(is_distance)"temporal distance increases/decreases" else "differentiation positive",
      zero_before=mean(z$before[,"zeta_presence"]==0),allocations=z$draws,exact_allocations=z$exact,
      identity_error=z$efficiency_error)
  }))
}
slice <- function(z,q) {
  r <- z
  for(n in c("before","after","support","valid"))r[[n]]<-z[[n]][q,,drop=FALSE]
  r$mean <- z$mean[q,,,drop=FALSE]
  r
}
study <- expand.grid(assemblage=assemblages,temporal_order=1:5,stringsAsFactors=FALSE)
job <- function(i) {
  task <- study[i,]
  g <- task$assemblage; m <- task$temporal_order
  a <- input_panel(g,"A",data,definition); n <- length(a$ids)
  R <- 3L
  budget <- as.integer(Sys.getenv("ST4_BUDGET","320"))
  draws <- as.integer(Sys.getenv("ST4_DRAWS","32"))
  summary <- list()
  combo <- list()
  for(offset in seq_len(6-m)) {
    positions <- offset+seq_len(m)-1L
    from <- a$years[positions]; to <- a$years[positions+5L]
    plan <- st4_prepare(a$data,a$ids,from,to)
    set.seed(1000+100*i+offset)
    sets <- lapply(1:6,function(k)lapply(seq_len(if(k==1)1L else R),function(r)draw_sets(n,k,budget)))
    all <- do.call(c,unlist(lapply(sets,function(s)lapply(s,st4_indices,n=n,m=m)),recursive=FALSE))
    masks <- st3_allocations(plan,32L,draws,7000+100*i+offset)
    full <- st4_evaluate_cpp(plan$before,plan$after,all,masks$masks,.st3_coalitions())
    colnames(full$before)<-colnames(full$after)<-colnames(full$support)<-metric_labels
    dimnames(full$phi)<-list(NULL,metric_labels,st3_processes,NULL)
    full$mean<-apply(full$phi,c(1,2,3),mean)
    full$valid<-full$support==1
    full$draws<-length(masks$masks);full$exact<-masks$exact
    cursor <- 0L
    for(k in 1:6)for(r in seq_along(sets[[k]])) {
      block <- sets[[k]][[r]]
      q <- cursor+seq_len(nrow(block));cursor<-max(q)
      z <- slice(full,q)
      sum_cross <- rowSums(matrix(z$mean[,c("occupancy","numerical"),],nrow(block),8))
      sum_total <- rowSums(matrix(z$mean[,"abundance",],nrow(block),4))
      stopifnot(all(abs(sum_cross-sum_total)<1e-10,na.rm=TRUE),
                all(z$mean[,c("zeta_presence","incidence","turnover","nestedness"),"numerical"]==0,na.rm=TRUE))
      summary[[length(summary)+1L]]<-save_rows(z,g,"A_C",k,m,offset,from,to,rep(r,nrow(block)),budget)
      if((k==2||k==1)&&m%in%c(1,5)) {
        ids <- apply(block,1,function(idx)paste(a$ids[idx],collapse="|"))
        combo[[length(combo)+1L]] <- do.call(rbind,lapply(c("abundance","occupancy","numerical","incidence",
                   "dominance_dissimilarity","site_temporal_presence","site_temporal_abundance",
                   "site_temporal_dominance","temporal_presence","temporal_abundance","temporal_dominance"),function(use) {
          if(k==1&&!grepl("temporal",use))return(NULL)
          data.frame(assemblage=g,k,m,offset,replicate=r,sites=ids,metric=use,
            before=z$before[,use],after=z$after[,use],support=z$support[,use],
            z$mean[,use,,drop=FALSE][,1,],row.names=NULL)
        }))
      }
    }
  }
  list(summary=do.call(rbind,summary),combo=do.call(rbind,combo))
}
all <- jobs(seq_len(nrow(study)),job)
save_csv(do.call(rbind,lapply(all,`[[`,"summary")),"process_replicates")
data.table::fwrite(do.call(rbind,lapply(all,`[[`,"combo")),result("process_combinations.csv.gz"),na="")
