source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
source(file.path(root,"scripts/engine.R"))
data <- archive(); def <- panels()
geo <- fread(processed("lake_geography.csv"))
job <- function(i) {
 g <- assemblages[i]
 a <- input_panel(g,"B",data,def)
 n <- length(a$ids)
 blocks <- geo$spatial_block[match(a$ids,geo$EU_ID)]
 set.seed(990+i)
 sets <- lapply(1:6,function(k)draw_sets(n,k,150))
 local <- matrix(seq_len(n),ncol=1L)
 take <- c(1L,2L,4L)
 all <- do.call(c,lapply(sets,st4_indices,n=n,m=1L))
 span <- c(0L,cumsum(lengths(lapply(sets,function(x)seq_len(nrow(x))))))
 trajectory <- annual <- site <- deletion <- list()
 for(t in seq_len(length(a$years)-1L)) {
   y <- a$years[t]
   plan <- st4_prepare(a$data,a$ids,y,y+1L)
   end <- st4_evaluate_cpp(plan$before,plan$after,all,list(plan$mask),matrix(0,16,0))
   colnames(end$before)<-colnames(end$after)<-metric_labels
   for(k in 2:6) {
     ids <- seq.int(span[k]+1L,span[k+1L])
     types <- c("zeta_presence","zeta_abundance","zeta_dominance","incidence","turnover","nestedness",
                "abundance","occupancy","numerical","dominance_dissimilarity")
     if(t==1L) trajectory[[length(trajectory)+1L]] <- data.frame(assemblage=g,year=y,k,metric=types,
                                  value=colMeans(end$before[ids,types,drop=FALSE]),sets=length(ids))
     trajectory[[length(trajectory)+1L]] <- data.frame(assemblage=g,year=y+1L,k,metric=types,
                                  value=colMeans(end$after[ids,types,drop=FALSE]),sets=length(ids))
   }
   # Evaluate a common allocation across both selected spatial orders and every local result.
   active <- do.call(c,lapply(sets[take],st4_indices,n=n,m=1L))
   al <- st3_allocations(plan,32,16,55000L+i*100+t)
   q <- st4_evaluate_cpp(plan$before,plan$after,active,al$masks,.st3_coalitions())
   colnames(q$before)<-colnames(q$after)<-colnames(q$support)<-metric_labels
   dimnames(q$phi)<-list(NULL,metric_labels,st3_processes,NULL)
   mean_phi<-apply(q$phi,c(1,2,3),mean)
   stopifnot(q$efficiency_error<1e-10)
   cursor <- 0L
   for(k in take) {
     ids<-cursor+seq_len(nrow(sets[[k]]));cursor<-max(ids)
     types <- c("site_temporal_presence","site_temporal_abundance","site_temporal_dominance",
                 "temporal_presence","temporal_abundance","temporal_dominance")
     if(k!=1)types<-c("abundance","occupancy","numerical","incidence","dominance_dissimilarity",types)
     for(metric in types) {
       f<-matrix(mean_phi[ids,metric,],length(ids),4L)
       valid<-complete.cases(f)
       mean_or_na<-function(a) if(sum(valid))mean(a[valid]) else NA_real_
       h<-colMeans(pmax(-f[valid,,drop=FALSE],0));d<-colMeans(pmax(f[valid,,drop=FALSE],0))
       annual[[length(annual)+1L]]<-data.frame(assemblage=g,year=y,end=y+1L,k,metric,process=st3_processes,
         defined=sum(valid),sets=length(ids),before=mean_or_na(q$before[ids,metric]),
         after=mean_or_na(q$after[ids,metric]),mean=colMeans(f[valid,,drop=FALSE]),negative=h,positive=d,
         allocations=length(al$masks),row.names=NULL)
       if(k==1L && startsWith(metric,"site_")) site[[length(site)+1L]]<-data.frame(
         assemblage=g,year=y,end=y+1L,id=a$ids,metric,distance=q$after[ids,metric],
         support=q$support[ids,metric],f,row.names=NULL)
       if(k%in%c(2L,4L)&&metric%in%c("abundance","occupancy","numerical","incidence","dominance_dissimilarity",
                 "site_temporal_abundance","site_temporal_presence")) {
         b <- matrix(blocks[sets[[k]]],nrow(sets[[k]]))
         for(block in seq_len(6L)) {
           retain <- valid & rowSums((b==block) & !is.na(b))==0L
           deletion[[length(deletion)+1L]]<-data.frame(assemblage=g,year=y,k,metric,block,
             process=st3_processes,mean=colMeans(f[retain,,drop=FALSE]),sets=sum(retain),row.names=NULL)
         }
       }
     }
   }
 }
 list(trajectory=do.call(rbind,trajectory),annual=do.call(rbind,annual),site=do.call(rbind,site),
      deletion=do.call(rbind,deletion))
}
z <- jobs(seq_along(assemblages),job)
for(n in c("trajectory","annual","site","deletion"))
 save_csv(do.call(rbind,lapply(z,`[[`,n)),paste0("B_",n))
