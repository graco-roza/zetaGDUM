# 95_order_diagnostics.R (part 2 helper): support machinery copied unchanged from figure6_gdm_synthesis.R (prepare(), cell support, matched-TP-profile screening) so that every evaluated order can be screened with exactly the Figure 6 rules.
# Not used to select orders by effect size; reports why orders are (or are not) usable. Sourced by 95_order_diagnostics.R.
POS<-c('Low','Intermediate','High');TARGET<-c(Low=10,Intermediate=35,High=100);COMP<-c('occupancy','heterogeneity','total')
MINREF<-10L;NREF<-1000L;MINWIDTH<-.20;GRID<-seq(0,1.50,by=.05); BINS<-12L; MINSET<-8L; MINSITE<-12L
OUT<-Sys.getenv('SUPPORT_OUT','revision_random/tables/order_support'); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
ns<-asNamespace('zetaGDUM');intern<-function(n)get(n,ns,inherits=FALSE)
cellfn<-intern('.component_cells'); traincell<-intern('.component_support_cells'); layoutfn<-intern('.interaction_layout')
tp_only<-function(d){for(v in setdiff(d$vars,'TP')){j<-match(v,d$vars);ii<-6L*(j-1L)+1:6;blk<-d$dir[,ii[1:3],drop=FALSE]+d$dir[,ii[4:6],drop=FALSE];d$dir[,ii[1:3]]<-blk/2;d$dir[,ii[4:6]]<-blk/2};d}
hav<-function(x1,y1,x2,y2){r<-pi/180;a<-sin((y2-y1)*r/2)^2+cos(y1*r)*cos(y2*r)*sin((x2-x1)*r/2)^2;2*6371*asin(sqrt(pmin(1,pmax(0,a))))}
set_geo<-function(xy,idx){v<-numeric(nrow(idx));k<-ncol(idx);for(i in seq_len(k-1L))for(j in seq.int(i+1L,k))v<-v+hav(xy[idx[,i],1],xy[idx[,i],2],xy[idx[,j],1],xy[idx[,j],2]);v/choose(k,2)}
raw_stats<-function(d,xy)cbind(d$raw_position,d$raw_spread,geographic_distance=set_geo(xy,d$idx))
spread_rows<-function(m){k<-ncol(m);v<-numeric(nrow(m));for(i in seq_len(k-1L))for(j in seq.int(i+1L,k))v<-v+abs(m[,i]-m[,j]);v/choose(k,2)}
prepare<-function(pool,tx,k,d){
 E<-data.frame(TP=log10(d$env$TP),TN=log10(d$env$TN),Chla=log10(d$env$Chla));xy<-as.matrix(d$env[,c('lon','lat')]);stopifnot(all(is.finite(as.matrix(E))),all(is.finite(xy)))
 idx<-sample_sets(nrow(E),k,6000,seed=7000+k)
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
 rank<-order(rowSums(sweep(sweep(nuis,2,med,'-'),2,ds,'/')^2));rank<-rank[uniq[rank]&spr[rank]>1e-10];central<-head(rank,max(min(NREF,length(rank)),floor(length(rank)/2)));stopifnot(length(central)>=MINREF)
 ri<-central[unique(round(seq(1,length(central),length.out=min(NREF,length(central)))))];ref<-ii[ri,,drop=FALSE];refstats<-stats[ri,,drop=FALSE]
 shape<-sweep(sweep(tp[ri,,drop=FALSE],1,rowMeans(tp[ri,,drop=FALSE]),'-'),1,spr[ri],'/')
 stopifnot(max(abs(rowMeans(shape)))<1e-9,max(abs(spread_rows(shape)-1))<1e-8)
 other_ok<-rep(TRUE,nrow(ref));for(v in setdiff(dd$vars,'TP')){j<-match(v,dd$vars);ct<-cells[[v]];cc<-cellfn(refstats[,j],refstats[,j+length(dd$vars)],BINS,ct$position_range,ct$spread_range)$cell;other_ok<-other_ok & !is.na(cc)&ct$sets[cc]>=MINSET&ct$sites[cc]>=MINSITE}
 other_ok[is.na(other_ok)]<-FALSE
 tpct<-cells[['TP']]; ok<-matrix(FALSE,nrow(ref),length(GRID)); rawok<-array(FALSE,c(nrow(ref),length(GRID),length(POS)))
 # Simultaneous per-predictor support + joint nearest-neighbour diagnostic in
 # raw position/spread/mean geographic distance. This does not prove support in
 # the full spline-feature space; it is a conservative extrapolation diagnostic.
 for(j in seq_along(GRID))for(p in seq_along(POS)){
   m<-log10(TARGET[p]);s<-GRID[j];newtp<-m+s*shape
   cc<-cellfn(rep(m,nrow(ref)),rep(s,nrow(ref)),BINS,tpct$position_range,tpct$spread_range)$cell
   inside<-apply(newtp,1,function(v)all(v>=min(E$TP)&v<=max(E$TP)))
   r<-refstats;r[,1]<-m;r[,1+length(dd$vars)]<-s
   nd<-RANN::nn2(trz,standard(r),k=1)$nn.dists[,1]
   z<-other_ok & inside & !is.na(cc) & tpct$sets[cc]>=MINSET & tpct$sites[cc]>=MINSITE & nd<=threshold+1e-12
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
 write.csv(data.frame(pool=pool,group=tx,order=k,spread=GRID,common_profiles=colSums(ok),low_profiles=colSums(rawok[,,1]),intermediate_profiles=colSums(rawok[,,2]),high_profiles=colSums(rawok[,,3]),nn_threshold=threshold),file.path(OUT,sprintf('support_%s_%s_k%d.csv',pool,tx,k)),row.names=FALSE)
 list(pool=pool,tx=tx,k=k,E=E,xy=xy,ref=ref,shape=shape,dd=dd,valid=valid,ok=ok,threshold=threshold)
}
