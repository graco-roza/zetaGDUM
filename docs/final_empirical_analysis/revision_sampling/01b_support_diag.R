# 01b_support_diag.R: copy of prepare_sim that also returns the per-criterion pass arrays (diagnostic only; nothing is relaxed in the simulation results).
prepare_diag<-function(E,xy,idx,k){
 stopifnot(all(is.finite(as.matrix(E))),all(is.finite(xy)))
 dd<-ugdm_component_design(E,idx,'additive',coordinates=xy,lonlat=TRUE)
 ii<-dd$idx;stats<-raw_stats(dd,xy);uniq<-!duplicated(t(apply(ii,1,sort)));tr<-unique(stats[uniq,,drop=FALSE]);mu<-colMeans(tr);sd0<-apply(tr,2,sd);sd0[!is.finite(sd0)|sd0==0]<-1
 standard<-function(x)sweep(sweep(x,2,mu,'-'),2,sd0,'/');trz<-standard(tr);stopifnot(nrow(trz)>1)
 threshold<-as.numeric(quantile(RANN::nn2(trz,trz,k=2)$nn.dists[,2],.95))
 trainref<-list(position=dd$raw_position,spread=dd$raw_spread,idx=ii)
 cells<-setNames(lapply(dd$vars,function(v)traincell(trainref,v,BINS)),dd$vars)
 # Reference profiles are observed distinct training configurations. Select a
 # systematic sample of the central half of non-TP environmental configurations.
 # Reference selection never uses fitted responses or coefficients.
 tp<-matrix(E$TP[ii],nrow(ii));spr<-spread_rows(tp);nuis<-stats[,c(2,3,5,6,7),drop=FALSE];med<-apply(nuis,2,median);ds<-apply(nuis,2,sd);ds[ds==0]<-1
 rank<-order(rowSums(sweep(sweep(nuis,2,med,'-'),2,ds,'/')^2));rank<-rank[uniq[rank]&spr[rank]>1e-10];central<-head(rank,max(min(NREF,length(rank)),floor(length(rank)/2)));if(length(central)<MINREF)return(list(k=k,E=E,xy=xy,valid=list(),ok=matrix(FALSE,1,length(GRID)),threshold=NA_real_,dd=dd,ref=NULL,shape=NULL))
 ri<-central[unique(round(seq(1,length(central),length.out=min(NREF,length(central)))))];ref<-ii[ri,,drop=FALSE];refstats<-stats[ri,,drop=FALSE]
 shape<-sweep(sweep(tp[ri,,drop=FALSE],1,rowMeans(tp[ri,,drop=FALSE]),'-'),1,spr[ri],'/')
 stopifnot(max(abs(rowMeans(shape)))<1e-9,max(abs(spread_rows(shape)-1))<1e-8)
 other_ok<-rep(TRUE,nrow(ref));for(v in setdiff(dd$vars,'TP')){j<-match(v,dd$vars);ct<-cells[[v]];cc<-cellfn(refstats[,j],refstats[,j+length(dd$vars)],BINS,ct$position_range,ct$spread_range)$cell;other_ok<-other_ok & !is.na(cc)&ct$sets[cc]>=MINSET&ct$sites[cc]>=MINSITE}
 other_ok[is.na(other_ok)]<-FALSE
 tpct<-cells[['TP']]; CR<-list(inside=array(FALSE,c(nrow(ref),length(GRID),length(POS))),cell=array(FALSE,c(nrow(ref),length(GRID),length(POS))),nn=array(FALSE,c(nrow(ref),length(GRID),length(POS)))); ok<-matrix(FALSE,nrow(ref),length(GRID)); rawok<-array(FALSE,c(nrow(ref),length(GRID),length(POS)))
 # Simultaneous per-predictor support + joint nearest-neighbour diagnostic in
 # raw position/spread/mean geographic distance. This does not prove support in
 # the full spline-feature space; it is a conservative extrapolation diagnostic.
 for(j in seq_along(GRID))for(p in seq_along(POS)){
   m<-log10(TARGET[p]);s<-GRID[j];newtp<-m+s*shape
   cc<-cellfn(rep(m,nrow(ref)),rep(s,nrow(ref)),BINS,tpct$position_range,tpct$spread_range)$cell
   inside<-apply(newtp,1,function(v)all(v>=min(E$TP)&v<=max(E$TP)))
   r<-refstats;r[,1]<-m;r[,1+length(dd$vars)]<-s
   nd<-RANN::nn2(trz,standard(r),k=1)$nn.dists[,1]
   cellok<-!is.na(cc) & tpct$sets[cc]>=MINSET & tpct$sites[cc]>=MINSITE; cellok[is.na(cellok)]<-FALSE; nnok<-nd<=threshold+1e-12
   CR$inside[,j,p]<-inside; CR$cell[,j,p]<-cellok; CR$nn[,j,p]<-nnok
   z<-other_ok & inside & cellok & nnok
   z[is.na(z)]<-FALSE;rawok[,j,p]<-z
 }
 ok<-apply(rawok,c(1,2),all)
 # Count profiles shared at ALL positions across every cell in candidate intervals.
 # Require >=3 grid points, >=.20 log-spread width and >=2 distinct TP-spread bins.
 valid<-list();for(i in seq_along(GRID))for(j in seq.int(i,length(GRID))){
   if(j-i+1L<3L||GRID[j]-GRID[i]+1e-10<MINWIDTH)next
   ids<-which(rowSums(ok[,i:j,drop=FALSE])==j-i+1L);if(length(ids)<MINREF)next
   cb<-cellfn(rep(log10(35),j-i+1L),GRID[i:j],BINS,tpct$position_range,tpct$spread_range)$spread_bin
   if(length(unique(cb[!is.na(cb)]))<2L)next
   valid[[paste(i,j,sep=':')]]<-ids
 }
  list(k=k,other_ok=other_ok,CR=CR,ok=ok,nref=nrow(ref))
}
