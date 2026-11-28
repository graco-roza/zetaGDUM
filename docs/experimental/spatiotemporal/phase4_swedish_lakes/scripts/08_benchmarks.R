source("docs/experimental/spatiotemporal/phase4_swedish_lakes/scripts/common.R")
source(file.path(root,"scripts/engine.R"))
if(!requireNamespace("ecopart",quietly=TRUE))stop("Install the Tatsumi et al. ecopart implementation.")
d <- archive(); definition <- panels()
x <- fread(result("process_combinations.csv.gz"))
bench <- validations <- local <- list()
for(g in assemblages) {
 a <- input_panel(g,"A",d,definition)
 for(j in 1:5) {
   plan <- st4_prepare(a$data,a$ids,a$years[j],a$years[j+5L])
   subset <- x[x$assemblage==g&x$m==1&x$k==2&x$offset==j&x$replicate==1&x$metric=="abundance",]
   subset <- head(subset,40)
   idx <- t(vapply(strsplit(subset$sites,"|",fixed=TRUE),function(s)match(s,a$ids),integer(2)))
   z <- st4_eval(plan,idx,4,741)
   if(j==1) {
     r <- st3_shapley(plan,idx,.st3_responses,exact_limit=32,draws=4,seed=741)
     fast <- z$mean[,.st3_responses,]
     both <- is.finite(fast)&is.finite(r$contributions)
     stopifnot(identical(is.finite(fast),is.finite(r$contributions)),
               max(abs(fast[both]-r$contributions[both]))<1e-10)
     val <- zetaGDUM::pair_partition(plan$before)
     lookup <- paste(val$i,val$j)
     key <- apply(idx,1,function(y)paste(sort(y),collapse=" "))
     e <- val[match(key,lookup),]
     base <- max(abs(z$before[,c("incidence","turnover","nestedness")]-as.matrix(e[,c("sorensen","simpson","nestedness")])))
     stopifnot(base<1e-12)
     validations[[length(validations)+1L]] <- data.frame(assemblage=g,n_sets=nrow(idx),
       native_phase3_max_error=max(abs(fast[both]-r$contributions[both])),baselga_max_error=base,
       phase3_kernel_efficiency=z$efficiency_error)
   }
   for(k in seq_len(nrow(idx)))for(metric in c("incidence","abundance")) {
     z0<-plan$before[idx[k,],];z1<-plan$after[idx[k,],]
     index<-if(metric=="abundance")"bray-curtis" else "sorensen"
     p<-ecopart::ecopart.pair(z0,z1,index=index,components="four")
     q<-x[x$assemblage==g &x$m==1&x$k==2&x$offset==j&x$replicate==1 &
           x$metric==metric & x$sites==subset$sites[k],][1,]
     if(all(is.finite(p)))stopifnot(abs(sum(p)-(q$after-q$before))<1e-9)
     phi<-as.numeric(unlist(q[1,st3_processes],use.names=FALSE))
     bench[[length(bench)+1L]]<-data.frame(assemblage=g,offset=j,lake_pair=q$sites,metric,
       change=q$after-q$before,tatsumi_net=sum(p),tatsumi_defined=all(is.finite(p)),
       tatsumi_loss_homogenisation=-p[1],tatsumi_loss_differentiation=p[2],
       tatsumi_gain_homogenisation=-p[3],tatsumi_gain_differentiation=p[4],
       expansion=phi[1],contraction=phi[2],redistribution=phi[3],relative_shares=phi[4],
       shapley_homogenisation=sum(pmax(-phi,0)),shapley_differentiation=sum(pmax(phi,0)),row.names=NULL)
   }
   # Legendre's incidence loss/gain decomposition on the complete local panel.
   P<-rowSums(plan$before>0 & plan$after>0)
   G<-rowSums(plan$before==0 & plan$after>0)
   L<-rowSums(plan$before>0 & plan$after==0)
   denom<-2*P+G+L
   z<-st4_eval(plan,matrix(seq_along(a$ids),ncol=1),endpoint=TRUE)
   D<-(G+L)/denom
   error<-max(abs(z$after[,"site_temporal_presence"]-D))
   stopifnot(error<1e-12)
   local[[length(local)+1L]]<-data.frame(assemblage=g,offset=j,id=a$ids,persistence=P,gain=G,loss=L,
                      gain_sorensen=G/denom,loss_sorensen=L/denom,distance=D,identity_error=error)
 }
}
save_csv(do.call(rbind,bench),"pairwise_tatsumi")
save_csv(do.call(rbind,validations),"mathematical_validation")
save_csv(do.call(rbind,local),"local_legendre")

window_validation <- lapply(assemblages,function(g){
  a<-input_panel(g,"A",d,definition)
  plan<-st4_prepare(a$data,a$ids,a$years[1:3],a$years[6:8])
  set.seed(8123)
  sets<-draw_sets(length(a$ids),4,10)
  fast<-st4_eval(plan,sets,4,9921)
  base<-st3_shapley(plan,sets,.st3_responses,exact_limit=32,draws=4,seed=9921)
  value<-fast$mean[,.st3_responses,]
  finite<-is.finite(value)
  stopifnot(identical(finite,is.finite(base$contributions)))
  error<-max(abs(value[finite]-base$contributions[finite]))
  stopifnot(error<1e-10)
  data.frame(assemblage=g,k=4,m=3,sets=nrow(sets),phase3_error=error,valid_cells=sum(finite))
})
save_csv(do.call(rbind,window_validation),"window_mathematical_validation")
