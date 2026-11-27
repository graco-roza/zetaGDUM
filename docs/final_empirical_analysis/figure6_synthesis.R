# Figure 6: eutrophication synthesis. Run from docs/final_empirical_analysis.
# ZG_LIB6=<f7e115f library> Rscript figure6_synthesis.R fig6_interaction/tables figures/figure6_synthesis
# No refitting. Uses the existing all-site coherent TP-only and additive coefficients.
# A: changes in predicted total, occupancy-associated, heterogeneity (percentage points).
# B: mean centred total-stage linear predictor from EXACT same matched TP profiles.
# Low k=2; high k=highest evaluated k with P(shared)>=.5 and common TP contrast.
# All choices are made from environment/support and sharing, not predicted effects.
.libPaths(c(Sys.getenv('ZG_LIB6'),.libPaths()))
suppressPackageStartupMessages({library(zetaGDUM);library(dplyr);library(tidyr);library(ggplot2);library(scales);library(grid);library(gtable)})
stopifnot(requireNamespace('RANN',quietly=TRUE),requireNamespace('svglite',quietly=TRUE))
a<-commandArgs(TRUE); IN<-if(length(a)>=1)a[1]else'fig6_interaction/tables'; OUT<-if(length(a)>=2)a[2]else'figures/figure6_synthesis';dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
PEN<-as.numeric(Sys.getenv('FIG6_PENALTY','1e-6'));stopifnot(PEN%in%c(.01,1e-4,1e-6))
TX<-c('BM','ZP');TXL<-c(BM='EPA macroinvertebrates',ZP='EPA zooplankton');POOLS<-c('shared','full')
POS<-c('Low','Intermediate','High');TARGET<-c(Low=10,Intermediate=35,High=100);COMP<-c('occupancy','heterogeneity','total')
MINREF<-10L;NREF<-1000L;MINWIDTH<-.20;GRID<-seq(0,1.50,by=.05); BINS<-12L; MINSET<-8L; MINSITE<-12L
# New display criterion, separate from original >=.10 informative-order rule.
MINSHARE<-.50
ns<-asNamespace('zetaGDUM');intern<-function(n)get(n,ns,inherits=FALSE)
cellfn<-intern('.component_cells'); traincell<-intern('.component_support_cells'); layoutfn<-intern('.interaction_layout')
orders<-read.csv('tables/evaluated_orders.csv',stringsAsFactors=FALSE) |> filter(community%in%unname(TXL),pool%in%POOLS) |> mutate(group=names(TXL)[match(community,TXL)],eligible=informative & sharing_probability>=MINSHARE)
write.csv(orders,file.path(OUT,'Figure6_order_information.csv'),row.names=FALSE)
# Select candidates on shared pool, reproduce exactly in full-pool sensitivity.
candidates<-setNames(lapply(TX,function(tx)sort(unique(orders$order[orders$group==tx & orders$pool=='shared' & orders$eligible]),decreasing=TRUE)),TX)
stopifnot(all(vapply(candidates,function(x)2%in%x,logical(1))))
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
# Main order determined from shared-pool sharing and support only. Full pool
# must use the identical order and contrast or be explicitly marked unsupported.
SEL<-list();screen<-list()
for(tx in TX){low<-getp('shared',tx,2);chosen<-NA_integer_;iv<-NA_character_
 for(k in candidates[[tx]][candidates[[tx]]>2]){
  high<-getp('shared',tx,k);kk<-intersect(names(low$valid),names(high$valid));b<-best_interval(kk)
  screen[[length(screen)+1L]]<-data.frame(group=tx,order=k,n_common_intervals=length(kk),best_interval=b)
  if(is.na(chosen)&&length(kk)){chosen<-k;iv<-b}
 }
 if(is.na(chosen)){warning(tx,': no common supported higher order; no effect will be invented.');iv<-best_interval(names(low$valid))}
 SEL[[tx]]<-list(high=chosen,interval=iv,k50=max(candidates[[tx]]))
}
write.csv(bind_rows(screen),file.path(OUT,'Figure6_candidate_order_support.csv'),row.names=FALSE)
selection<-bind_rows(lapply(TX,function(tx){s<-SEL[[tx]];ij<-if(is.na(s$interval))c(NA_integer_,NA_integer_)else parse_interval(s$interval);data.frame(group=tx,k50=s$k50,k_selected=s$high,spread_from=GRID[ij[1]],spread_to=GRID[ij[2]],sharing_at_selected=orders$sharing_probability[match(paste(tx,'shared',s$high),paste(orders$group,orders$pool,orders$order))],selection_rule='Highest k>2 with P(shared)>=0.50 and same supported interval as k=2')}))
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
  k<-ks[slot];iv<-ss$interval;available<-!is.na(k)&&!is.na(iv)&&iv%in%names(getp(pool,tx,k)$valid)
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
write.csv(C,file.path(OUT,'Figure6_B_curves.csv'),row.names=FALSE);write.csv(B,file.path(OUT,'Figure6_A_bars.csv'),row.names=FALSE)
# -------------------- journal-style figure, explicit gtable alignment --------------------
ink<-'#252B2D';gridcol<-'#E4E8E8';grey<-'#848B8D';poscol<-c(Low='#79936C',Intermediate='#246878',High='#C45C19');compcol<-c(occupancy='#246878',heterogeneity='#C45C19',total='#615A83')
complab<-c(occupancy='Occupancy-associated',heterogeneity='Abundance heterogeneity',total='Total dissimilarity');poslab<-c(Low='Low (10)',Intermediate='Intermediate (35)',High='High (100)')
base<-theme_minimal(base_size=9,base_family='sans')+theme(text=element_text(colour=ink),axis.title=element_blank(),axis.text=element_text(size=7.8,colour=ink),axis.ticks=element_line(colour=ink,linewidth=.25),axis.ticks.length=unit(1,'mm'),axis.line=element_line(colour=ink,linewidth=.3),panel.grid.major=element_line(colour=gridcol,linewidth=.25),panel.grid.minor=element_blank(),legend.position='none',plot.margin=margin(1,1.5,1,1.5))
yspec<-function(x,symmetric=FALSE){x<-x[is.finite(x)];if(!length(x))x<-c(-.1,.1);r<-range(c(0,x));if(diff(r)<1e-12)r<-c(-.01,.01);if(symmetric)r<-c(-1,1)*max(abs(r));lim<-r+c(-1,1)*.06*diff(r);br<-sort(unique(c(0,breaks_extended(4)(lim))));list(lim=lim,br=br[br>=lim[1]&br<=lim[2]])}
# Scales use the actual model units, identical across order rows and pools within
# each taxon. Taxa have different scales. A's x-axis is signed pp; B's y is logit.
bar_specs<-setNames(lapply(TX,function(tx)yspec(B$change_pp[B$group==tx],TRUE)),TX)
line_specs<-setNames(lapply(TX,function(tx)yspec(C$delta_eta[C$group==tx])),TX)
xs<-setNames(lapply(TX,function(tx){iv<-SEL[[tx]]$interval;if(is.na(iv))c(0,1)else GRID[parse_interval(iv)]}),TX)
bar_plot<-function(pool,tx,slot){
 v<-B |> filter(.data$pool==.env$pool,group==tx,.data$slot==.env$slot) |> mutate(y=8-3.25*(match(position,POS)-1)-.8*(match(component,COMP)-2),xmin=pmin(0,change_pp),xmax=pmax(0,change_pp))
 sp<-bar_specs[[tx]];p<-ggplot(v)+geom_vline(xintercept=0,colour=grey,linewidth=.45)+geom_rect(aes(xmin=xmin,xmax=xmax,ymin=y-.30,ymax=y+.30,fill=component))+scale_fill_manual(values=compcol)
 if(!nrow(v))p<-p+annotate('text',x=0,y=4.75,label='Insufficient common\nenvironmental support',size=3,colour=grey)
 p+scale_y_continuous(limits=c(.4,9.1),breaks=c(8,4.75,1.5),labels=c('Low','Intermediate','High'),expand=expansion(mult=0))+scale_x_continuous(limits=sp$lim,breaks=sp$br,expand=expansion(mult=0))+base+theme(panel.grid.major.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
}
line_plot<-function(pool,tx,slot){
 v<-C |> filter(.data$pool==.env$pool,group==tx,.data$slot==.env$slot);sp<-line_specs[[tx]];xl<-xs[[tx]]
 p<-ggplot(v,aes(spread,delta_eta,colour=position,linetype=model,group=interaction(position,model)))+geom_hline(yintercept=0,colour=grey,linewidth=.35)+geom_line(data=filter(v,model=='Additive'),linewidth=.65)+geom_line(data=filter(v,model=='Interaction'),linewidth=.95)+scale_colour_manual(values=poscol)+scale_linetype_manual(values=c(Additive='42',Interaction='solid'))
 if(!nrow(v))p<-p+annotate('text',x=mean(xl),y=mean(sp$lim),label='Insufficient common\nenvironmental support',size=3,colour=grey)
 p+scale_x_continuous(limits=xl,breaks=breaks_extended(4),expand=expansion(mult=.035))+scale_y_continuous(limits=sp$lim,breaks=sp$br,expand=expansion(mult=0))+base
}
label<-function(s,size=9,bold=FALSE,rot=0,hjust=.5,x=.5)textGrob(s,x=unit(x,'npc'),hjust=hjust,rot=rot,gp=gpar(col=ink,fontsize=size,fontface=if(bold)'bold'else'plain'))
legend<-function(cols,labs,lines=FALSE){xs<-c(.02,.35,.69);kids<-list();for(i in seq_along(cols)){if(lines)mark<-segmentsGrob(unit(xs[i],'npc'),unit(.5,'npc'),unit(xs[i]+.052,'npc'),unit(.5,'npc'),gp=gpar(col=cols[i],lwd=2))else mark<-rectGrob(x=unit(xs[i]+.024,'npc'),width=unit(.046,'npc'),height=unit(2.6,'mm'),gp=gpar(col=NA,fill=cols[i]));kids<-c(kids,list(mark,label(labs[i],8,x=xs[i]+.062,hjust=0)))};do.call(grobTree,kids)}
make_figure<-function(pool){
 plots<-list();audit<-list();for(section in c('A','B'))for(slot in 1:2)for(tx in TX){p<-if(section=='A')bar_plot(pool,tx,slot)else line_plot(pool,tx,slot);b<-ggplot_build(p)$layout$panel_params[[1]];audit[[length(audit)+1L]]<-data.frame(pool=pool,section=section,group=tx,slot=slot,xmin=b$x.range[1],xmax=b$x.range[2],ymin=b$y.range[1],ymax=b$y.range[2]);plots[[paste(section,slot,tx)]]<-ggplotGrob(p)}
 # Match all underlying plot widths, regardless of tick-label widths, plus heights
 # within each pair. Exactly equal panel columns with conventional axes.
 stopifnot(length(unique(vapply(plots,function(g)length(g$widths),integer(1))))==1L)
 ww<-do.call(unit.pmax,lapply(plots,function(g)g$widths));for(i in names(plots))plots[[i]]$widths<-ww
 for(section in c('A','B'))for(slot in 1:2){ix<-paste(section,slot,TX);hh<-unit.pmax(plots[[ix[1]]]$heights,plots[[ix[2]]]$heights);for(i in ix)plots[[i]]$heights<-hh}
 g<-gtable(widths=unit.c(unit(8,'mm'),unit(1,'null'),unit(4,'mm'),unit(1,'null')),heights=unit.c(unit(7,'mm'),unit(7,'mm'),unit(5,'mm'),unit(1.25,'null'),unit(5,'mm'),unit(1.25,'null'),unit(6,'mm'),unit(6,'mm'),unit(4,'mm'),unit(7,'mm'),unit(5,'mm'),unit(1,'null'),unit(5,'mm'),unit(1,'null'),unit(6,'mm'),unit(6,'mm'),unit(5,'mm')))
 add<-function(o,t,l,b=t,r=l)g<<-gtable_add_grob(g,o,t=t,l=l,b=b,r=r,clip='off')
 add(label('A',12,TRUE),1,1);add(label('Component responses to phosphorus heterogeneity',11,TRUE,hjust=0,x=0),1,2,r=4)
 add(label('B',12,TRUE),10,1);add(label('Context-dependent sensitivity to phosphorus heterogeneity',11,TRUE,hjust=0,x=0),10,2,r=4)
 for(j in 1:2){tx<-TX[j];cl<-2*j;add(label(TXL[tx],10.5,TRUE),2,cl)
  for(slot in 1:2){k<-if(slot==1)2 else SEL[[tx]]$high;heading<-if(is.na(k))paste0('k50 = ',SEL[[tx]]$k50,' (unsupported)')else paste0('k = ',k);ra<-if(slot==1)3 else 5;rb<-if(slot==1)11 else 13
   add(label(heading,9,TRUE),ra,cl);add(plots[[paste('A',slot,tx)]],ra+1,cl)
   add(label(heading,9,TRUE),rb,cl);add(plots[[paste('B',slot,tx)]],rb+1,cl)
  }
 }
 add(label('TP position',9,rot=90),4,1,b=6)
 add(label('Change in dissimilarity (percentage points)',9),7,2,r=4);add(legend(compcol,complab),8,2,r=4)
 add(label('TP spread contribution\n(logit total dissimilarity)',9,rot=90),12,1,b=14)
 add(label('TP spread (mean pairwise absolute difference in log10 TP)',9),15,2,r=4)
 add(legend(poscol,poslab,TRUE),16,2,r=4);add(label(expression('TP position ('*mu*'g/L); solid: TP interaction; dashed: additive'),8),17,2,r=4)
 audit<-bind_rows(audit);stopifnot(all(is.finite(as.matrix(audit[,c('xmin','xmax','ymin','ymax')]))))
 for(sec in c('A','B'))for(tx in TX){aa<-audit[audit$section==sec & audit$group==tx,c('xmin','xmax','ymin','ymax')];stopifnot(nrow(aa)==2L,max(abs(unlist(aa[1,],use.names=FALSE)-unlist(aa[2,],use.names=FALSE)))<1e-8)}
 write.csv(audit,file.path(OUT,paste0('Figure6_axis_checks_',pool,'.csv')),row.names=FALSE)
 for(ext in c('pdf','svg','png'))ggsave(file.path(OUT,paste0('Figure6_synthesis_',pool,'.',ext)),g,width=183,height=216,units='mm',dpi=600,bg='white',device=switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png'))
}
for(pool in POOLS)make_figure(pool)
message('Done. Inspect selected orders, panel support and all values before interpretation.')
