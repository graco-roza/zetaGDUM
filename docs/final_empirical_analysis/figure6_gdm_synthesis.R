# Figure 6: eutrophication synthesis, GDM-style transformation panels. Run from docs/final_empirical_analysis.
# ZG_LIB6=<f7e115f library> Rscript figure6_gdm_synthesis.R fig6_interaction/tables figures/figure6_gdm_synthesis
# Panel A reuses existing coherent fits. Panel B fits independent marginal components.
# A: changes in predicted total, occupancy-associated, heterogeneity (percentage points).
# B: fitted environmental transformations f(TP | position), not raw spread contrasts.
# Predeclared orders: BM 2,5; ZP 2,4. No effect-based order selection.
# All choices are made from environment/support and sharing, not predicted effects.
.libPaths(c(Sys.getenv('ZG_LIB6'),.libPaths()))
suppressPackageStartupMessages({library(zetaGDUM);library(dplyr);library(tidyr);library(ggplot2);library(scales);library(grid);library(gtable)})
stopifnot(requireNamespace('RANN',quietly=TRUE),requireNamespace('svglite',quietly=TRUE),requireNamespace('ggh4x',quietly=TRUE))
a<-commandArgs(TRUE); IN<-if(length(a)>=1)a[1]else'fig6_interaction/tables'; OUT<-if(length(a)>=2)a[2]else'figures/figure6_gdm_synthesis';dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
PEN<-as.numeric(Sys.getenv('FIG6_PENALTY','1e-6'));stopifnot(PEN%in%c(.01,1e-4,1e-6))
TX<-c('BM','ZP');TXL<-c(BM='EPA macroinvertebrates',ZP='EPA zooplankton');POOLS<-c('shared','full')
POS<-c('Low','Intermediate','High');TARGET<-c(Low=10,Intermediate=35,High=100);COMP<-c('occupancy','heterogeneity','total')
MINREF<-10L;NREF<-1000L;MINWIDTH<-.20;GRID<-seq(0,1.50,by=.05); BINS<-12L; MINSET<-8L; MINSITE<-12L
# Fixed ecological contrast; keep original >=.10 informative-order criterion.
ORD<-list(BM=c(2L,5L),ZP=c(2L,4L))
# Optional explicitly labelled alternative, never selected by effect sizes.
if(nzchar(Sys.getenv('FIG6_ZP_HIGH'))) ORD$ZP[2]<-as.integer(Sys.getenv('FIG6_ZP_HIGH'))
stopifnot(ORD$ZP[2]%in%c(3L,4L,5L))
ns<-asNamespace('zetaGDUM');intern<-function(n)get(n,ns,inherits=FALSE)
cellfn<-intern('.component_cells'); traincell<-intern('.component_support_cells'); layoutfn<-intern('.interaction_layout')
orders<-read.csv('tables/evaluated_orders.csv',stringsAsFactors=FALSE) |> filter(community%in%unname(TXL),pool%in%POOLS) |> mutate(group=names(TXL)[match(community,TXL)],eligible=informative)
write.csv(orders,file.path(OUT,'Figure6_order_information.csv'),row.names=FALSE)
# Fixed orders, checked against the existing informative-order table.
candidates<-ORD
for(tx in TX)for(pool in POOLS)stopifnot(all(ORD[[tx]] %in% orders$order[orders$group==tx & orders$pool==pool & orders$eligible]))
# Quantify taxonomic representation independently of model coefficients.
# These exact taxon contribution tables were produced by Figure 3.
rep_rows<-list()
for(pool in POOLS){
 fn<-file.path('tables',paste0('fig3_taxon_contributions_',pool,'.csv'))
 if(!file.exists(fn))stop('Missing Figure 3 taxon contribution table: ',fn)
 v<-read.csv(fn,stringsAsFactors=FALSE) |> filter(community %in% unname(TXL)) |>
   mutate(group=names(TXL)[match(community,TXL)])
 for(cutoff in c(.25,.50,.75)){
   rr<-v |> group_by(community,group,order) |> summarise(
    n_taxa=n(), n_widespread=sum(occupancy>=cutoff),
    weighted_occupancy_similarity=sum(occupancy*sim_rel)/sum(sim_rel),
    weighted_occupancy_dissimilarity=sum(occupancy*dis_rel)/sum(dis_rel),
    widespread_similarity=sum(sim_rel[occupancy>=cutoff])/sum(sim_rel),
    widespread_dissimilarity=sum(dis_rel[occupancy>=cutoff])/sum(dis_rel),.groups='drop') |>
    mutate(pool=pool,occupancy_threshold=cutoff)
   rep_rows[[length(rep_rows)+1L]]<-rr
 }
}
REP<-bind_rows(rep_rows)
REP$displayed<-mapply(function(g,k) k %in% ORD[[g]],REP$group,REP$order)
write.csv(REP,file.path(OUT,'Figure6_taxonomic_representation.csv'),row.names=FALSE)
coef_tab<-bind_rows(lapply(POOLS,function(p)read.csv(file.path(IN,paste0('interaction_coefficients_allsites_',p,'.csv')),stringsAsFactors=FALSE))) |> filter(model=='additive' | (variant=='interaction_TP' & abs(penalty-PEN)<1e-12))
link_tab<-bind_rows(lapply(POOLS,function(p)read.csv(file.path(IN,paste0('interaction_link_decomposition_',p,'.csv')),stringsAsFactors=FALSE))) |> filter(model=='additive' | (variant=='interaction_TP' & abs(penalty-PEN)<1e-12))
# Package TP-only implementation, identical to 61_spline_shapes.R.
tp_only<-function(d){for(v in setdiff(d$vars,'TP')){j<-match(v,d$vars);ii<-6L*(j-1L)+1:6;blk<-d$dir[,ii[1:3],drop=FALSE]+d$dir[,ii[4:6],drop=FALSE];d$dir[,ii[1:3]]<-blk/2;d$dir[,ii[4:6]]<-blk/2};d}
hav<-function(x1,y1,x2,y2){r<-pi/180;a<-sin((y2-y1)*r/2)^2+cos(y1*r)*cos(y2*r)*sin((x2-x1)*r/2)^2;2*6371*asin(sqrt(pmin(1,pmax(0,a))))}
set_geo<-function(xy,idx){v<-numeric(nrow(idx));k<-ncol(idx);for(i in seq_len(k-1L))for(j in seq.int(i+1L,k))v<-v+hav(xy[idx[,i],1],xy[idx[,i],2],xy[idx[,j],1],xy[idx[,j],2]);v/choose(k,2)}
raw_stats<-function(d,xy)cbind(d$raw_position,d$raw_spread,geographic_distance=set_geo(xy,d$idx))
spread_rows<-function(m){k<-ncol(m);v<-numeric(nrow(m));for(i in seq_len(k-1L))for(j in seq.int(i+1L,k))v<-v+abs(m[,i]-m[,j]);v/choose(k,2)}
# Rebuild same training basis and sampled sets; no fit or response-dependent choices.
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
cache<-list();key<-function(pool,tx,k)paste(pool,tx,k,sep='_')
for(pool in POOLS){d<-readRDS(file.path(paste0('fig6_',pool),'out','prepared.rds'));for(tx in TX)for(k in candidates[[tx]]){
 message('Screen support: ',pool,' ',tx,' k=',k);cache[[key(pool,tx,k)]]<-prepare(pool,tx,k,d[[tx]])
}}
getp<-function(pool,tx,k)cache[[key(pool,tx,k)]]
parse_interval<-function(s)as.integer(strsplit(s,':',fixed=TRUE)[[1]])
best_interval<-function(keys){if(!length(keys))return(NA_character_);z<-t(vapply(keys,parse_interval,integer(2)));keys[order(-(z[,2]-z[,1]),z[,1],z[,2])[1]]}
# Use a common spread interval at both FIXED orders. Prefer an interval also
# supported in both full pools within each taxon; if none, use the widest
# shared-pool interval and explicitly mark missing full-pool panels.
# If no shared-pool common interval exists, do not manufacture a comparison.
# In that case display only a supported k=2 contrast and mark high k unavailable.
SEL<-list();screen<-list()
for(tx in TX){
 k<-ORD[[tx]][2]; low<-getp('shared',tx,2); high<-getp('shared',tx,k)
 kk<-intersect(names(low$valid),names(high$valid))
 both<-Reduce(intersect,list(kk,names(getp('full',tx,2)$valid),names(getp('full',tx,k)$valid)))
 iv<-best_interval(if(length(both))both else kk)
 matched<-!is.na(iv)
 rule<-if(length(both))'Common interval: both orders and pools' else if(matched)'Common interval: both shared-pool orders' else 'No matched higher-order contrast; k=2 only'
 if(!matched){warning(tx,': no common supported interval for fixed orders; higher-order bars unavailable.');iv<-best_interval(names(low$valid))}
 SEL[[tx]]<-list(high=k,interval=iv,matched=matched,rule=rule)
 screen[[length(screen)+1L]]<-data.frame(group=tx,order=k,n_common_shared=length(kk),n_common_both_pools=length(both),selected_interval=iv,matched=matched,selection_rule=rule)
}
write.csv(bind_rows(screen),file.path(OUT,'Figure6_candidate_order_support.csv'),row.names=FALSE)
selection<-bind_rows(lapply(TX,function(tx){s<-SEL[[tx]];ij<-if(is.na(s$interval))c(NA_integer_,NA_integer_)else parse_interval(s$interval);data.frame(group=tx,k_selected=s$high,matched=s$matched,spread_from=GRID[ij[1]],spread_to=GRID[ij[2]],sharing_at_selected=orders$sharing_probability[match(paste(tx,'shared',s$high),paste(orders$group,orders$pool,orders$order))],selection_rule=s$rule)}))
write.csv(selection,file.path(OUT,'Figure6_selected_orders.csv'),row.names=FALSE)
print(selection)
# Evaluate exact synthetic configurations only after support/order selection.
# Each profile keeps TN, Chla, location and TP profile shape fixed. TP geometric
# mean and mean absolute pairwise log difference are set by construction.
# Rebuild synthetic evaluation design with original reference; multiply by the
# exported full-layout coefficients, not by refitted coefficients.
coef_for<-function(p,kind,stage,X){
 vv<-if(kind=='Additive')'additive'else'interaction_TP'
 z<-coef_tab |> filter(pool==p$pool,group==p$tx,order==p$k,.data$stage==.env$stage,variant==vv)
 stopifnot(nrow(z)>0,all(z$converged),!any(z$stage_constant),!anyDuplicated(z$term),setequal(z$term,colnames(X)))
 unname(setNames(z$coef,z$term)[colnames(X)])
}
# Verify reconstruction on all-site training fit, against exported site spline
# coefficient vectors and the package's full model layout. Batch synthetic sets
# to avoid very large pairwise geographic matrices.
eval_cell<-function(p,ids,position,s,kind){
 ref<-p$ref[ids,,drop=FALSE];sh<-p$shape[ids,,drop=FALSE];k<-p$k;ans<-list()
 for(b in split(seq_along(ids),ceiling(seq_along(ids)/80L))){
  ii<-as.vector(t(ref[b,,drop=FALSE]));en<-p$E[ii,,drop=FALSE];xy<-p$xy[ii,,drop=FALSE];en$TP<-as.vector(t(log10(TARGET[position])+s*sh[b,,drop=FALSE]));idx<-matrix(seq_len(length(ii)),ncol=k,byrow=TRUE)
  zz<-ugdm_component_design(en,idx,if(kind=='Additive')'additive'else'interaction',coordinates=xy,lonlat=TRUE,reference=p$reference_design)
  if(kind=='Interaction')zz<-tp_only(zz)
  stopifnot(max(abs(zz$raw_position[,'TP']-log10(TARGET[position])))<1e-7,max(abs(zz$raw_spread[,'TP']-s))<1e-7)
  # Verify unchanged nuisance statistics against the corresponding original sets.
  rr<-p$ref[ids[b],,drop=FALSE]; nmat<-function(v)matrix(v[rr],nrow(rr)); for(nm in c('TN','Chla')){stopifnot(max(abs(zz$raw_position[,nm]-rowMeans(nmat(p$E[[nm]]))))<1e-7,max(abs(zz$raw_spread[,nm]-spread_rows(nmat(p$E[[nm]]))))<1e-7)}
  X<-layoutfn(zz,!is.null(zz$geo))$X
  et<-drop(X%*%coef_for(p,kind,'total',X));eh<-drop(X%*%coef_for(p,kind,'share',X));tt<-plogis(et);hh<-tt*plogis(eh)
  ans[[length(ans)+1L]]<-data.frame(profile=ids[b],eta_total=et,eta_share=eh,total=tt,heterogeneity=hh,occupancy=tt-hh)
 };bind_rows(ans)
}
# Interaction evaluation must reference an interaction training design. Additive
# and interaction share the original environmental bases, but model labels must
# be rebuilt explicitly. Keep separately reconstructed training references.
curve_rows<-list();profile_rows<-list();status<-list()
for(pool in POOLS)for(tx in TX){ss<-SEL[[tx]];ks<-c(2L,ss$high)
 for(slot in seq_along(ks)){
  k<-ks[slot];iv<-ss$interval;available<-!is.na(k)&&!is.na(iv)&&(slot==1L||ss$matched)&&iv%in%names(getp(pool,tx,k)$valid)
  status[[length(status)+1L]]<-data.frame(pool=pool,group=tx,slot=slot,order=k,available=available,n_profiles=if(available)length(getp(pool,tx,k)$valid[[iv]])else 0L,interval=iv)
  if(!available){warning(pool,' ',tx,' slot ',slot,': no common supported contrast.');next}
  p<-getp(pool,tx,k);ids<-p$valid[[iv]];ij<-parse_interval(iv);sp<-GRID[ij[1]:ij[2]]
  write.csv(data.frame(profile=ids,p$ref[ids,,drop=FALSE]),file.path(OUT,sprintf('profiles_%s_%s_k%d.csv',pool,tx,k)),row.names=FALSE)
  # Reconstruct interaction training design using exactly the original indices.
  di<-ugdm_component_design(p$E,p$dd$idx,'interaction',coordinates=p$xy,lonlat=TRUE)
  for(kind in c('Additive','Interaction')){
   # The original training reference is either additive or interaction, and
   # the original 61 script fitted TP-only by collapsing TN and Chla columns.
   p$reference_design<-if(kind=='Additive')p$dd else di
   train_design<-if(kind=='Additive')p$dd else tp_only(di)
   Xtrain<-layoutfn(train_design,!is.null(train_design$geo))$X
   # Verify reconstructed coefficients/design against the original exported
   # median linear predictor at the 200 nearest training sets.
   vv<-if(kind=='Additive')'additive'else'interaction_TP'
   lv<-link_tab |> filter(.data$pool==.env$pool,group==tx,order==k,variant==vv)
   stopifnot(nrow(lv)==3L,all(lv$converged))
   meanTP<-rowMeans(matrix(p$E$TP[p$dd$idx],nrow(p$dd$idx)))
   for(st in c('total','share')){eta<-drop(Xtrain%*%coef_for(p,kind,st,Xtrain));for(pos in POS){near<-order(abs(meanTP-log10(TARGET[pos])))[1:200];original<-lv[[paste0('eta0_',st)]][lv$position==pos];stopifnot(length(original)==1L,abs(median(eta[near])-original)<1e-7)}}
   for(pos in POS)for(s in sp){
    v<-eval_cell(p,ids,pos,s,kind) |> mutate(pool=pool,group=tx,slot=slot,order=k,model=kind,position=pos,TP_position_ugL=unname(TARGET[pos]),spread=s)
    profile_rows[[length(profile_rows)+1L]]<-v
   }
  }
 }
}
STATUS<-bind_rows(status);write.csv(STATUS,file.path(OUT,'Figure6_panel_support.csv'),row.names=FALSE)
if(!length(profile_rows))stop('No supported profiles; see support exports. No figure created.')
P<-bind_rows(profile_rows);stopifnot(all(is.finite(P$eta_total)),all(abs(P$total-P$occupancy-P$heterogeneity)<1e-10))
P<-P |> group_by(pool,group,slot,order,model,position,profile) |> arrange(spread,.by_group=TRUE) |> mutate(delta_eta=eta_total-first(eta_total),across(all_of(COMP),~100*(.x-first(.x)),.names='delta_{.col}')) |> ungroup()
write.csv(P,file.path(OUT,'Figure6_profile_predictions.csv'),row.names=FALSE)
C<-P |> group_by(pool,group,slot,order,model,position,TP_position_ugL,spread) |> summarise(n_profiles=n_distinct(profile),across(c(delta_eta,delta_total,delta_occupancy,delta_heterogeneity),mean),.groups='drop')
stopifnot(max(abs(C$delta_total-C$delta_occupancy-C$delta_heterogeneity))<1e-8)
B<-C |> filter(model=='Interaction') |> group_by(pool,group,slot,order,position,TP_position_ugL) |> arrange(spread,.by_group=TRUE) |> summarise(spread_from=first(spread),spread_to=last(spread),n_profiles=first(n_profiles),total=last(delta_total),occupancy=last(delta_occupancy),heterogeneity=last(delta_heterogeneity),.groups='drop') |> pivot_longer(all_of(COMP),names_to='component',values_to='change_pp')
write.csv(C,file.path(OUT,'FigureS6_matched_spread_curves.csv'),row.names=FALSE);write.csv(B,file.path(OUT,'Figure6_A_bars.csv'),row.names=FALSE)
# ------------------ Panel B: genuine fitted environmental transformations ------------------
# Independent fits are needed for marginal occupancy and heterogeneity f(TP).
# Their predictions need not close. Coherent Panel A above remains unchanged.
# Total is independently refitted and checked against the existing coherent fit.
# All models reuse the original 6000 combinations, seeds, spline bases, ridge,
# predictor set and TP-only interaction construction. No sign constraints on
# occupancy/heterogeneity; their signed f is not a nonnegative ecological distance.
STAGES<-c('total','occupancy','heterogeneity')
STLAB<-c(total='Total',occupancy='Occupancy-associated',heterogeneity='Abundance heterogeneity')
# Identical central site-concentration interval for both orders and pools per taxon.
DOM<-setNames(lapply(TX,function(tx){
 qq<-lapply(POOLS,function(pool)quantile(getp(pool,tx,2)$E$TP,c(.02,.98),names=FALSE))
 c(max(vapply(qq,`[`,numeric(1),1L)),min(vapply(qq,`[`,numeric(1),2L)))
}),TX)
stopifnot(all(vapply(DOM,function(z)all(is.finite(z))&&diff(z)>0,logical(1))))
write.csv(bind_rows(lapply(TX,function(tx)data.frame(group=tx,TP_min=10^DOM[[tx]][1],TP_max=10^DOM[[tx]][2]))),file.path(OUT,'Figure6_B_site_domain.csv'),row.names=FALSE)
SROWS<-list();FCOEF<-list();FITCHECK<-list();SUPP<-list()
for(pool in POOLS){
 D<-readRDS(file.path(paste0('fig6_',pool),'out','prepared.rds'))
 for(tx in TX)for(k in ORD[[tx]]){
  p<-getp(pool,tx,k); dd0<-p$dd; y<-ugdm_component_response(D[[tx]]$x,dd0$idx)
  # Same basis as the fitted additive design; interaction uses same training sites.
  zg<-seq(DOM[[tx]][1],DOM[[tx]][2],length.out=300)
  for(kind in c('Additive','Interaction')){
   dd<-if(kind=='Additive')dd0 else tp_only(ugdm_component_design(p$E,dd0$idx,'interaction',coordinates=p$xy,lonlat=TRUE))
   X<-layoutfn(dd,!is.null(dd$geo))$X
   fit<-ugdm_fit_components(dd,y,architecture='independent',ridge=.001,
       interaction_shrink=if(kind=='Additive').01 else PEN,monotone_total=TRUE)
   stopifnot(identical(fit$architecture,'independent'),setequal(names(fit$fits),STAGES))
   tb<-dd$basis; ctr<-tb$center[['TP']]; scl<-tb$scale[['TP']]
   spl<-tb$splines[[match('TP',dd$vars)]]; bnd<-attr(spl,'Boundary.knots')
   Bz<-function(z)unclass(predict(spl,pmin(pmax((z-ctr)/scl,bnd[1]),bnd[2])))[,1:3,drop=FALSE]
   uof<-function(z)pmin(1,pmax(0,((z-ctr)/scl-bnd[1])/diff(bnd)))
   bz<-Bz(zg); stopifnot(nrow(bz)==length(zg),all(is.finite(bz)))
   # Test the transformation against the actual TP directional contribution
   # from the training layout (including the correct pairwise sign for each basis).
   ii<-dd$idx; ip<-seq_len(min(nrow(ii),100L)); zmat<-matrix(p$E$TP[ii[ip,,drop=FALSE]],nrow=length(ip))
   bmat<-array(Bz(as.vector(zmat)),dim=c(nrow(zmat),ncol(zmat),3L))
   spr<-sapply(1:3,function(j)spread_rows(matrix(bmat[,,j],nrow=nrow(zmat))))
   uu<-uof(rowMeans(zmat))
   for(st in STAGES){
    f<-fit$fits[[st]]; b<-setNames(drop(fit$transform %*% f$coef),colnames(X))
    constant<-isTRUE(f$constant)||is.finite(f$constant_value)
    converged<-isTRUE(f$converged)
    FCOEF[[length(FCOEF)+1L]]<-data.frame(pool=pool,group=tx,order=k,model=kind,stage=st,term=names(b),coef=unname(b),converged=converged,constant=constant,penalty=if(kind=='Interaction')PEN else NA_real_)
    # The same objective and constraints must reproduce the coherent total fit.
    total_err<-NA_real_; total_coef_err<-NA_real_
    if(st=='total'){
     bc<-coef_for(p,kind,'total',X)
     total_err<-max(abs(plogis(drop(X%*%b))-plogis(drop(X%*%bc))))
     total_coef_err<-max(abs(b-bc));stopifnot(converged,!constant,total_err<1e-6,total_coef_err<1e-5)
    }
    if(kind=='Additive'){
     tl<-paste0('TP_f',1:3);tu<-tl;betaL<-betaU<-unname(b[tl])
     from_design<-drop(X[ip,tl,drop=FALSE]%*%betaL);from_basis<-drop(spr%*%betaL)
    }else{
     tl<-paste0('TP_f',1:3,'_L');tu<-paste0('TP_f',1:3,'_U')
     betaL<-unname(b[tl]);betaU<-unname(b[tu])
     from_design<-drop(X[ip,c(tl,tu),drop=FALSE]%*%c(betaL,betaU))
     from_basis<-(1-uu)*drop(spr%*%betaL)+uu*drop(spr%*%betaU)
    }
    stopifnot(all(is.finite(c(betaL,betaU))),max(abs(from_design-from_basis))<1e-7)
    FITCHECK[[length(FITCHECK)+1L]]<-data.frame(pool=pool,group=tx,order=k,model=kind,stage=st,converged=converged,constant=constant,coherent_total_prediction_error=total_err,coherent_total_coef_error=total_coef_err,directional_reconstruction_error=max(abs(from_design-from_basis)))
    if(!converged||constant){warning('No fitted curve exported: ',pool,' ',tx,' ',k,' ',kind,' ',st);next}
    for(pos in POS){
     # Do not present a target TP position outside training set mean range.
     lx<-log10(TARGET[[pos]]);pos_ok<-lx>=min(dd$raw_position[,'TP'])&&lx<=max(dd$raw_position[,'TP'])
     uu0<-uof(lx);bb<-(1-uu0)*betaL+uu0*betaU
     fv<-drop(bz%*%bb);v<-fv-fv[1]
     if(st=='total')stopifnot(all(diff(v)>=-1e-8))
     # One position-independent additive f; retain one reference curve only.
     if(kind=='Additive'&&pos!='Low')next
     if(kind=='Interaction'&&!pos_ok){warning('Target TP position outside training range: ',pool,' ',tx,' ',k,' ',pos);next}
     SROWS[[length(SROWS)+1L]]<-data.frame(pool=pool,group=tx,order=k,stage=st,model=kind,position=if(kind=='Additive')'Additive'else pos,TP_position_ugL=if(kind=='Additive')NA_real_ else TARGET[[pos]],u=if(kind=='Additive')NA_real_ else uu0,TP_ugL=10^zg,TP_log10=zg,f=fv,f_centred=v,reference_TP=10^zg[1],beta1=bb[1],beta2=bb[2],beta3=bb[3])
    }
   }
   if(kind=='Interaction'){
    # Report support of each fixed mean TP in training TP-position bins. The
    # displayed f is a fitted basis transformation, not a synthetic set or a
    # supported full-range raw-spread prediction. Export counts, not a CI.
    ct<-traincell(list(position=dd$raw_position,spread=dd$raw_spread,idx=dd$idx),'TP',BINS)
    for(pos in POS){lx<-log10(TARGET[[pos]]); cc<-cellfn(rep(lx,BINS),seq(ct$spread_range[1],ct$spread_range[2],length.out=BINS),BINS,ct$position_range,ct$spread_range)$cell
     SUPP[[length(SUPP)+1L]]<-data.frame(pool=pool,group=tx,order=k,position=pos,TP_position_ugL=TARGET[[pos]],n_supported_spread_cells=sum(!is.na(cc)&ct$sets[cc]>=MINSET&ct$sites[cc]>=MINSITE,na.rm=TRUE),position_in_range=lx>=min(dd$raw_position[,'TP'])&&lx<=max(dd$raw_position[,'TP']))
    }
   }
  }
 }
}
S<-bind_rows(SROWS);stopifnot(nrow(S)>0,all(is.finite(S$f_centred)))
write.csv(S,file.path(OUT,'Figure6_B_GDM_spline_values.csv'),row.names=FALSE)
write.csv(bind_rows(FCOEF),file.path(OUT,'Figure6_B_independent_coefficients.csv'),row.names=FALSE)
write.csv(bind_rows(FITCHECK),file.path(OUT,'Figure6_B_fit_checks.csv'),row.names=FALSE)
write.csv(bind_rows(SUPP),file.path(OUT,'Figure6_B_position_support.csv'),row.names=FALSE)
# ------------------ Figure assembly ------------------
ink<-'#252B2D';gridcol<-'#E4E8E8';grey<-'#848B8D'
poscol<-c(Low='#79936C',Intermediate='#246878',High='#C45C19',Additive=grey)
compcol<-c(occupancy='#246878',heterogeneity='#C45C19',total='#615A83')
complab<-c(occupancy='Occupancy-associated',heterogeneity='Abundance heterogeneity',total='Total dissimilarity')
poslab<-c(Low='Low (10)',Intermediate='Intermediate (35)',High='High (100)',Additive='Additive')
base<-theme_minimal(base_size=8.5,base_family='sans')+theme(text=element_text(colour=ink),axis.title=element_blank(),axis.text=element_text(size=7,colour=ink),axis.ticks=element_line(colour=ink,linewidth=.25),axis.ticks.length=unit(1,'mm'),axis.line=element_line(colour=ink,linewidth=.30),panel.grid.major=element_line(colour=gridcol,linewidth=.25),panel.grid.minor=element_blank(),legend.position='none',plot.margin=margin(1,1.5,1,1.5))
yspec<-function(x,symmetric=FALSE){x<-x[is.finite(x)];if(!length(x))x<-c(-.1,.1);r<-range(c(0,x));if(diff(r)<1e-12)r<-c(-.01,.01);if(symmetric)r<-c(-1,1)*max(abs(r));lim<-r+c(-1,1)*.07*diff(r);br<-sort(unique(c(0,breaks_extended(4)(lim))));list(lim=lim,br=br[br>=lim[1]&br<=lim[2]])}
bar_specs<-setNames(lapply(TX,function(tx)yspec(B$change_pp[B$group==tx],TRUE)),TX)
# Panel A: exactly three components at each TP position. Use explicit centers
# and rectangle extents; fix previous clipping at the top and bottom.
centres<-c(Low=8,Intermediate=4.75,High=1.5)
yall<-as.vector(outer(centres,-.8*(seq_along(COMP)-2),'+'))
bar_y_lim<-range(yall)+c(-.48,.48)
BARCHECK<-list()
bar_plot<-function(pool,tx,slot){
 v<-B |> filter(.data$pool==.env$pool,group==tx,.data$slot==.env$slot) |> mutate(y=unname(centres[position])-.8*(match(component,COMP)-2),xmin=pmin(0,change_pp),xmax=pmax(0,change_pp),ymin=y-.30,ymax=y+.30)
 sp<-bar_specs[[tx]]
 stopifnot(nrow(v)%in%c(0L,9L),all(v$ymin>=bar_y_lim[1]),all(v$ymax<=bar_y_lim[2]))
 q<-ggplot(v)+geom_vline(xintercept=0,colour=grey,linewidth=.4)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=component))+scale_fill_manual(values=compcol)
 if(!nrow(v))q<-q+annotate('text',x=0,y=4.75,label='Insufficient common\nenvironmental support',size=2.7,colour=grey)
 q<-q+scale_y_continuous(limits=bar_y_lim,breaks=unname(centres),labels=names(centres),expand=expansion(mult=0))+scale_x_continuous(limits=sp$lim,breaks=sp$br,expand=expansion(mult=0))+base+theme(panel.grid.major.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
 gg<-ggplot_build(q);rd<-gg$data[[2]];stopifnot(nrow(rd)==nrow(v),all(is.finite(rd$xmin)),all(is.finite(rd$xmax)),all(is.finite(rd$ymin)),all(is.finite(rd$ymax)))
 BARCHECK[[length(BARCHECK)+1L]]<<-data.frame(pool=pool,group=tx,slot=slot,n_input=nrow(v),n_rendered=nrow(rd),ymin=bar_y_lim[1],ymax=bar_y_lim[2],xmin=sp$lim[1],xmax=sp$lim[2])
 q
}
# Panel B: nested strips for taxon and order. Each taxon has a separate facet
# block so free_y shares limits across its two orders but not across taxa.
# Invisible anchors enforce the same limits across shared and full pools.
line_block<-function(pool,tx,logaxis=FALSE){
 k<-ORD[[tx]]; xlim<-10^DOM[[tx]]
 decorate<-function(x)x |> mutate(component=factor(stage,levels=STAGES,labels=unname(STLAB)),order_label=factor(paste0('k = ',order),levels=paste0('k = ',k)),taxon_label=unname(TXL[[tx]]),position=factor(position,levels=names(poscol)))
 allv<-S |> filter(group==tx) |> decorate(); v<-allv |> filter(.data$pool==.env$pool)
 anchors<-bind_rows(lapply(STAGES,function(st){sp<-yspec(allv$f_centred[allv$stage==st]);expand.grid(stage=st,order=k,TP_ugL=xlim,f_centred=sp$lim,stringsAsFactors=FALSE)})) |> mutate(position='Additive') |> decorate()
 counts<-v |> group_by(stage,order,position) |> summarise(n=n(),.groups='drop');stopifnot(all(counts$n==300L))
 # Explicit annotation if a model or TP position failed or lacks range support.
 missing<-expand.grid(stage=STAGES,order=k,position=names(poscol),stringsAsFactors=FALSE) |> anti_join(counts,by=c('stage','order','position'))
 if(nrow(missing))warning('Some GDM curves unavailable: ',pool,' ',tx)
 missing_note<-missing |> group_by(stage,order) |> summarise(note=paste('Unavailable:',paste(position,collapse=', ')),.groups='drop') |> mutate(TP_ugL=if(logaxis)10^mean(DOM[[tx]])else mean(xlim),f_centred=0,position='Additive') |> decorate()
 q<-ggplot(v,aes(TP_ugL,f_centred))+geom_blank(data=anchors)+geom_hline(yintercept=0,colour=grey,linewidth=.25)+
   geom_line(data=filter(v,position=='Additive'),colour=grey,linetype='42',linewidth=.55)+
   geom_line(data=filter(v,position!='Additive'),aes(colour=position,group=position),linewidth=.68,lineend='round')+
   scale_colour_manual(values=poscol)+
   ggh4x::facet_nested(rows=vars(component),cols=vars(taxon_label,order_label),scales='free_y',drop=FALSE)+
   scale_y_continuous(breaks=breaks_extended(3),expand=expansion(mult=0))+
   base+theme(panel.spacing=unit(2.4,'mm'),strip.background=element_rect(fill='#F0F2F2',colour=NA),strip.text.x=element_text(size=7.5,face='bold',margin=margin(2,1,2,1)),strip.text.y=element_text(size=5.2,angle=-90,face='bold',margin=margin(2,1,2,1)),axis.text=element_text(size=6.6))
 if(nrow(missing_note))q<-q+geom_text(data=missing_note,aes(label=note),size=2,colour=grey,vjust=-.5)
 if(logaxis){br<-c(1,3,10,30,100,300,1000,3000);br<-br[br>=xlim[1]&br<=xlim[2]];q<-q+scale_x_log10(limits=xlim,breaks=br,labels=label_number(),expand=expansion(mult=.025))}
 else q<-q+scale_x_continuous(limits=xlim,breaks=breaks_extended(3),labels=label_number(),expand=expansion(mult=.025))
 z<-ggplot_build(q);stopifnot(length(z$layout$panel_params)==6L)
 la<-z$layout$layout
 audit<-bind_rows(lapply(seq_len(nrow(la)),function(i){pr<-z$layout$panel_params[[as.integer(la$PANEL[i])]];data.frame(pool=pool,group=tx,component=as.character(la$component[i]),order=as.character(la$order_label[i]),logaxis=logaxis,ymin=pr$y.range[1],ymax=pr$y.range[2],xmin=pr$x.range[1],xmax=pr$x.range[2])}))
 for(st in unique(audit$component)){w<-audit[audit$component==st,];stopifnot(nrow(w)==2,max(abs(w$ymin-w$ymin[1]))<1e-8,max(abs(w$ymax-w$ymax[1]))<1e-8)}
 list(plot=q,audit=audit)
}
label<-function(s,size=9,bold=FALSE,rot=0,hjust=.5,x=.5)textGrob(s,x=unit(x,'npc'),hjust=hjust,rot=rot,gp=gpar(col=ink,fontsize=size,fontface=if(bold)'bold'else'plain'))
legend<-function(cols,labs,lines=FALSE){xs<-if(length(cols)==3)c(.02,.35,.69)else c(.02,.27,.57,.82);kids<-list();for(i in seq_along(cols)){if(lines)mark<-segmentsGrob(unit(xs[i],'npc'),unit(.5,'npc'),unit(xs[i]+.043,'npc'),unit(.5,'npc'),gp=gpar(col=cols[i],lwd=2,lty=if(names(cols)[i]=='Additive')2 else 1))else mark<-rectGrob(x=unit(xs[i]+.024,'npc'),width=unit(.046,'npc'),height=unit(2.6,'mm'),gp=gpar(col=NA,fill=cols[i]));kids<-c(kids,list(mark,label(labs[i],7.5,x=xs[i]+.055,hjust=0)))};do.call(grobTree,kids)}
export<-function(g,stem,h){for(ext in c('pdf','svg','png'))ggsave(file.path(OUT,paste0(stem,'.',ext)),g,width=183,height=h,units='mm',dpi=600,bg='white',device=switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png'))}
ALLAUDIT<-list()
for(pool in POOLS){
 # Same bar grob widths in each order row and each taxon, fixed clipping above.
 bars<-list();for(slot in 1:2)for(tx in TX)bars[[paste(slot,tx)]]<-ggplotGrob(bar_plot(pool,tx,slot))
 stopifnot(length(unique(vapply(bars,function(g)length(g$widths),integer(1))))==1L)
 ww<-do.call(unit.pmax,lapply(bars,function(g)g$widths));for(n in names(bars))bars[[n]]$widths<-ww
 for(slot in 1:2){ix<-paste(slot,TX);hh<-unit.pmax(bars[[ix[1]]]$heights,bars[[ix[2]]]$heights);for(n in ix)bars[[n]]$heights<-hh}
 for(logaxis in c(FALSE,TRUE)){
  br<-lapply(TX,function(tx)line_block(pool,tx,logaxis));ALLAUDIT[[length(ALLAUDIT)+1L]]<-bind_rows(lapply(br,`[[`,'audit'))
  gg<-lapply(br,function(x)ggplotGrob(x$plot))
  # Both nested blocks have the same structural widths/heights (two columns,
  # three rows). Synchronize gutters for exact alignment, keeping independent scales.
  if(length(gg[[1]]$widths)==length(gg[[2]]$widths)){w<-unit.pmax(gg[[1]]$widths,gg[[2]]$widths);for(j in 1:2)gg[[j]]$widths<-w}
  if(length(gg[[1]]$heights)==length(gg[[2]]$heights)){h<-unit.pmax(gg[[1]]$heights,gg[[2]]$heights);for(j in 1:2)gg[[j]]$heights<-h}
  # Standalone GDM-like panel, in addition to the complete synthesis figure.
  gb<-gtable(widths=unit.c(unit(8,'mm'),unit(1,'null'),unit(3,'mm'),unit(1,'null')),heights=unit.c(unit(1,'null'),unit(6,'mm'),unit(5,'mm'),unit(6,'mm')))
  gb<-gtable_add_grob(gb,gg,t=1,l=c(2,4),clip='off');gb<-gtable_add_grob(gb,label('Centred fitted f(TP | position)\n(logit scale)',8,rot=90),t=1,l=1,clip='off')
  gb<-gtable_add_grob(gb,label(if(logaxis)expression('Site TP ('*mu*'g/L, log axis)')else expression('Site TP ('*mu*'g/L)'),9),t=2,l=2,r=4,clip='off')
  gb<-gtable_add_grob(gb,label(expression('Geometric mean TP position ('*mu*'g/L)'),8),t=3,l=2,r=4,clip='off')
  gb<-gtable_add_grob(gb,legend(poscol,poslab,TRUE),t=4,l=2,r=4,clip='off')
  suffix<-if(logaxis)'_logaxis'else''
  export(gb,paste0('Figure6_GDM_splines_',pool,suffix),112)
  # A and B combined: A 2x2; B 3x4; small multiples with real data units.
  g<-gtable(widths=unit.c(unit(8,'mm'),unit(1,'null'),unit(3,'mm'),unit(1,'null')),heights=unit.c(unit(7,'mm'),unit(7,'mm'),unit(5,'mm'),unit(1,'null'),unit(5,'mm'),unit(1,'null'),unit(5,'mm'),unit(6,'mm'),unit(6,'mm'),unit(7,'mm'),unit(2.5,'null'),unit(6,'mm'),unit(5,'mm'),unit(6,'mm')))
  add<-function(o,t,l,b=t,r=l)g<<-gtable_add_grob(g,o,t=t,l=l,b=b,r=r,clip='off')
  add(label('A',11,TRUE),1,1);add(label('Component responses to phosphorus heterogeneity',10,TRUE,hjust=0,x=0),1,2,r=4)
  for(j in 1:2){tx<-TX[j];cl<-2*j;add(label(TXL[tx],9.5,TRUE),2,cl)
   for(slot in 1:2){r<-if(slot==1)3 else 5;add(label(paste0('k = ',ORD[[tx]][slot]),8.5,TRUE),r,cl);add(bars[[paste(slot,tx)]],r+1,cl)}
  }
  add(label('TP position',8,rot=90),4,1,b=6)
  add(label('Change in predicted dissimilarity (percentage points)',8.5),7,2,r=4);add(legend(compcol,complab),8,2,r=4)
  # Explicitly report distinct raw-spread contrasts. k=2-only is not an order comparison.
  note<-paste(vapply(TX,function(tx){z<-selection[selection$group==tx,];paste0(if(tx=='BM')'BM'else'ZP',': ',if(is.na(z$spread_from))'unsupported'else sprintf('%.2f-%.2f',z$spread_from,z$spread_to),if(!z$matched)' (k=2 only)'else'')},character(1)),collapse='; ')
  add(label(paste0('TP spread interval (log10 units): ',note),7.5),9,2,r=4)
  add(label('B',11,TRUE),10,1);add(label('Position-dependent fitted phosphorus transformations',10,TRUE,hjust=0,x=0),10,2,r=4)
  add(gg[[1]],11,2);add(gg[[2]],11,4);add(label('Centred fitted f(TP | position)\n(logit scale)',8,rot=90),11,1)
  add(label(if(logaxis)expression('Site TP ('*mu*'g/L, log axis)')else expression('Site TP ('*mu*'g/L)'),8.5),12,2,r=4)
  add(label(expression('Geometric mean TP position ('*mu*'g/L)'),8),13,2,r=4);add(legend(poscol,poslab,TRUE),14,2,r=4)
  export(g,paste0('Figure6_GDM_synthesis_',pool,suffix),227)
 }
}
AUDIT<-bind_rows(ALLAUDIT);write.csv(AUDIT,file.path(OUT,'Figure6_B_axis_checks.csv'),row.names=FALSE)
write.csv(bind_rows(BARCHECK),file.path(OUT,'Figure6_A_bar_checks.csv'),row.names=FALSE)
message('Finished. Panel A uses coherent predictions; B independent marginal transformations. Inspect support, convergence and all rendered panels.')
