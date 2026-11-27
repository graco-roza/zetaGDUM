# Conditional Figure 6: exact matched synthetic site-set profiles.
# Rscript figure6_conditional.R figures/figure6_redesign figures/figure6_conditional
# Run from docs/final_empirical_analysis, with ZG_LIB6 pointing at zetaGDUM f7e115f.
# Fits the original coherent additive model (same settings/seeds as 04_fit.R).
# Uses original observed k-lake covariate/coordinate configurations as references.
# Vary ONLY log10 TP, by shifting/rescaling its within-set profile to the desired
# position and spread. TN, Chla, coordinates and TP shape remain fixed exactly.
# Retain a COMMON set of reference profiles supported at every displayed cell.
# These are synthetic conditional predictions, not observed or causal effects.
.libPaths(c(Sys.getenv('ZG_LIB6'),.libPaths()))
suppressPackageStartupMessages({library(zetaGDUM);library(dplyr);library(tidyr)})
a<-commandArgs(TRUE);input<-if(length(a)>=1)a[1]else'figures/figure6_redesign';out<-if(length(a)>=2)a[2]else'figures/figure6_conditional';dir.create(out,recursive=TRUE,showWarnings=FALSE)
NREF<-as.integer(Sys.getenv('FIG6_NREF','120'));MINREF<-as.integer(Sys.getenv('FIG6_MINREF','10'));DO_FOLDS<-Sys.getenv('FIG6_FOLDS','1')=='1';stopifnot(NREF>=MINREF,MINREF>=2)
POS<-c('Oligotrophic','Mesotrophic','Eutrophic');COMP<-c('total','occupancy','heterogeneity');TX<-c('BM','ZP')
# Optional RANN speeds up nearest-neighbour diagnostics; exact base-R fallback.
nearest<-function(train,query,self=FALSE){
 if(requireNamespace('RANN',quietly=TRUE))return(RANN::nn2(train,query,k=if(self)2 else 1)$nn.dists[,if(self)2 else 1])
 ans<-numeric(nrow(query)); qn<-rowSums(query^2);tn<-rowSums(train^2)
 for(lo in seq(1,nrow(query),by=128)){ii<-lo:min(nrow(query),lo+127);d<-outer(qn[ii],tn,'+')-2*tcrossprod(query[ii,,drop=FALSE],train);if(self)d[cbind(seq_along(ii),ii)]<-Inf;ans[ii]<-sqrt(pmax(0,apply(d,1,min)))};ans
}
hav<-function(x1,y1,x2,y2){r<-pi/180;a<-sin((y2-y1)*r/2)^2+cos(y1*r)*cos(y2*r)*sin((x2-x1)*r/2)^2;2*6371*asin(sqrt(pmin(1,pmax(0,a))))}
set_geo<-function(xy,idx){v<-numeric(nrow(idx)); k<-ncol(idx);for(i in 1:(k-1))for(j in (i+1):k)v<-v+hav(xy[idx[,i],1],xy[idx[,i],2],xy[idx[,j],1],xy[idx[,j],2]);v/choose(k,2)}
raw_stats<-function(des,xy){cbind(des$raw_position,des$raw_spread,geographic_distance=set_geo(xy,des$idx)) |> unname()}
# Training-derived joint support in raw position, spread and geographic distance.
# This supplements the package's per-predictor support diagnostic; not a proof
# of full high-dimensional support. Duplicate statistics removed for NN threshold.
make_support<-function(fit,design,xy){r<-unique(raw_stats(design,xy));mu<-colMeans(r);sig<-apply(r,2,sd);sig[!is.finite(sig)|sig==0]<-1;z<-sweep(sweep(r,2,mu,'-'),2,sig,'/');stopifnot(nrow(z)>1);thr<-as.numeric(quantile(nearest(z,z,TRUE),.95));list(mu=mu,sig=sig,z=z,thr=thr)}
check_support<-function(fit,design,xy,sup){r<-raw_stats(design,xy);z<-sweep(sweep(r,2,sup$mu,'-'),2,sup$sig,'/');nd<-nearest(sup$z,z);pk<-ugdm_component_support(fit,design,bins=12,min_sets=8,min_sites=12)$supported;list(ok=pk&is.finite(nd)&nd<=sup$thr+1e-12,nn=nd,package=pk)}
get_grid<-function(pool,k,tx){stem<-if(k==3)if(pool=='shared')'Figure6_main_redesign'else'FigureS6_full_pool_redesign'else paste0('FigureS6_TP_B5_',pool);x<-read.csv(file.path(input,paste0(stem,'_B_supported_predictions.csv'))) |> filter(group==tx,order==k,model=='additive',component=='total');stopifnot(setequal(unique(x$trophic),POS)); pp<-x |> distinct(trophic,TP_position);stopifnot(nrow(pp)==3); sp<-Reduce(intersect,lapply(split(x$spread,x$trophic),function(z)round(z,8)));stopifnot(length(sp)>=3);expand.grid(trophic=POS,spread=sort(sp),stringsAsFactors=FALSE) |> left_join(pp,by='trophic') |> arrange(trophic,spread)}
# Set-level nuisance profiles selected independently of biological responses.
# Draw from observed set combinations, rank by proximity to median nuisance stats
# (TN, Chla, their spreads and mean geographic distance) to target typical profiles.
reference_sets<-function(E,xy,k){idx<-sample_sets(nrow(E),k,8000,seed=9000+k);idx<-matrix(as.integer(idx),ncol=k); idx<-idx[!duplicated(t(apply(idx,1,sort))),,drop=FALSE];m<-function(v)rowMeans(matrix(v[idx],nrow(idx)));spr<-function(v){r<-numeric(nrow(idx));for(i in 1:(k-1))for(j in (i+1):k)r<-r+abs(v[idx[,i]]-v[idx[,j]]);r/choose(k,2)};tp<-spr(E$TP);nuis<-cbind(m(E$TN),m(E$Chla),spr(E$TN),spr(E$Chla),set_geo(xy,idx));med<-apply(nuis,2,median);sd0<-apply(nuis,2,sd);sd0[sd0==0]<-1;rank<-order(rowSums(sweep(sweep(nuis,2,med,'-'),2,sd0,'/')^2));rank<-rank[tp[rank]>1e-10];# Systematically sample across the central half of nuisance profiles, rather than
 # concentrating exclusively at the median. No outcome is used.
 central<-head(rank,max(NREF,floor(length(rank)/2)));rank<-central[unique(round(seq(1,length(central),length.out=min(NREF,length(central)))))];stopifnot(length(rank)>=MINREF);idx[rank,,drop=FALSE]}
# One synthetic dataset per cell; original indices correspond to fixed profiles.
# Rebuilding the design with the fit as reference evaluates the original spline
# bases at the actual constructed TP concentrations. No design-column averaging.
synthetic<-function(E,xy,ref,pos,spread){k<-ncol(ref);n<-nrow(ref);i<-as.vector(t(ref));en<-E[i,,drop=FALSE];coords<-xy[i,,drop=FALSE];m<-matrix(en$TP,nrow=n,byrow=TRUE);v<-numeric(n);for(a in 1:(k-1))for(b in (a+1):k)v<-v+abs(m[,a]-m[,b]);v<-v/choose(k,2);stopifnot(all(v>0));m<-sweep(sweep(m,1,rowMeans(m),'-'),1,spread/v,'*')+pos;en$TP<-as.vector(t(m));list(env=en,xy=coords,idx=matrix(seq_len(n*k),nrow=n,byrow=TRUE),tp=m)}
fit_one<-function(d,E,xy,k,tr,fit_id,ref,grid){
 # Same 6000 combinations and fixed seed as the existing all-site curves.
 # Deletion fits are sensitivity analyses, not held-out skill estimates.
 et<-E[tr,,drop=FALSE];xt<-xy[tr,,drop=FALSE];it<-sample_sets(length(tr),k,6000,seed=7000+k);des<-ugdm_component_design(et,it,'additive',coordinates=xt,lonlat=TRUE);y<-ugdm_component_response(d$x[tr,,drop=FALSE],it);fit<-ugdm_fit_components(des,y,architecture='coherent');if(!isTRUE(fit$converged))return(NULL);sup<-make_support(fit,des,xt);rows<-list()
 for(j in seq_len(nrow(grid))){g<-grid[j,];sy<-synthetic(E,xy,ref,log10(g$TP_position),g$spread);zz<-ugdm_component_design(sy$env,sy$idx,'additive',coordinates=sy$xy,lonlat=TRUE,reference=fit);stopifnot(max(abs(zz$raw_position[,'TP']-log10(g$TP_position)))<1e-7,max(abs(zz$raw_spread[,'TP']-g$spread))<1e-7);ss<-check_support(fit,zz,sy$xy,sup);pred<-ugdm_predict_components(fit,zz);inside<-apply(sy$tp,1,function(v)all(v>=min(et$TP)&v<=max(et$TP)));stopifnot(all(abs(pred$total-pred$occupancy-pred$heterogeneity)<1e-8));rows[[j]]<-cbind(data.frame(fit_id=fit_id,profile=seq_len(nrow(ref)),trophic=g$trophic,TP_position=g$TP_position,spread=g$spread,supported=ss$ok&inside,package_support=ss$package,inside_TP_range=inside,nn_distance=ss$nn,nn_threshold=sup$thr),pred[COMP])};bind_rows(rows)
}
# Find widest contiguous spread interval with >=3 grid cells and >=MINREF
# profiles jointly supported at every position. Support-only rule, never uses
# responses or predicted effect size. Identical profiles for all curves.
select_common<-function(r,grid){sp<-sort(unique(grid$spread));best<-NULL;for(i in seq_along(sp))for(j in seq(i,length(sp))){if(j-i+1<3)next; x<-r |> filter(spread%in%sp[i:j]);ok<-x |> group_by(profile) |> summarise(ok=n()==3*(j-i+1)&&all(supported),.groups='drop');ids<-ok$profile[ok$ok];if(length(ids)<MINREF)next;width<-sp[j]-sp[i];if(is.null(best)||width>best$width+1e-12||(abs(width-best$width)<1e-12&&length(ids)>length(best$ids)))best<-list(width=width,sp=sp[i:j],ids=ids)};best}
run_group<-function(pool,k,tx,d){E<-data.frame(TP=log10(d$env$TP),TN=log10(d$env$TN),Chla=log10(d$env$Chla));xy<-as.matrix(d$env[,c('lon','lat')]);stopifnot(all(is.finite(as.matrix(E))),all(is.finite(xy)));grid<-get_grid(pool,k,tx);ref<-reference_sets(E,xy,k);all<-fit_one(d,E,xy,k,seq_len(nrow(E)),'all',ref,grid);if(is.null(all))stop('Primary fit failed: ',pool,' ',tx,' k=',k)
 common<-select_common(all,grid);filebase<-file.path(out,paste0(pool,'_',tx,'_k',k));write.csv(all,paste0(filebase,'_all_profile_diagnostics.csv'),row.names=FALSE);write.csv(data.frame(profile=seq_len(nrow(ref)),ref),paste0(filebase,'_reference_sets.csv'),row.names=FALSE)
 if(is.null(common)){warning('No common support at >=3 spread cells and >=',MINREF,' profiles: ',pool,' ',tx,' k=',k);return(list(summary=NULL,diagnostic=data.frame(pool=pool,group=tx,order=k,available=FALSE,n_profiles=0,n_cells=0,n_valid_deletion_fits=0)))}
 summ<-function(z){z |> filter(profile%in%common$ids,spread%in%common$sp) |> group_by(fit_id,trophic,TP_position,spread) |> summarise(n_profiles=n_distinct(profile),supported=all(supported),across(all_of(COMP),mean),.groups='drop') |> mutate(pool=pool,group=tx,order=k,model='additive')}
 res<-list(summ(all));valid<-0L
 if(DO_FOLDS){# Same complete-linkage spatial blocks as 04_fit.R.
   n<-nrow(xy);dm<-matrix(0,n,n);for(i in seq_len(n))dm[i,]<-hav(xy[i,1],xy[i,2],xy[,1],xy[,2]);fold<-cutree(hclust(as.dist(dm),'complete'),5)
   for(f in sort(unique(fold))){z<-tryCatch(fit_one(d,E,xy,k,which(fold!=f),paste0('delete_',f),ref[common$ids,,drop=FALSE],grid |> filter(spread%in%common$sp)),error=function(e){warning(conditionMessage(e));NULL});if(is.null(z))next;z$profile<-common$ids[z$profile];write.csv(z,paste0(filebase,'_delete',f,'_diagnostics.csv'),row.names=FALSE);if(all(z$supported)){res[[length(res)+1]]<-summ(z);valid<-valid+1L}}
 }
 list(summary=bind_rows(res),diagnostic=data.frame(pool=pool,group=tx,order=k,available=TRUE,n_profiles=length(common$ids),n_cells=length(common$sp),spread_from=min(common$sp),spread_to=max(common$sp),n_valid_deletion_fits=valid))
}
DIAG<-list();for(pool in c('shared','full')){p<-file.path(paste0('fig6_',pool),'out','prepared.rds');if(!file.exists(p))stop('Missing ',p,'; run original data setup.');D<-readRDS(p);for(k in c(3,5)){a<-lapply(TX,function(tx)run_group(pool,k,tx,D[[tx]]));s<-bind_rows(lapply(a,'[[','summary'));di<-bind_rows(lapply(a,'[[','diagnostic'));DIAG[[length(DIAG)+1]]<-di;if(nrow(s)){stopifnot(all(abs(s$total-s$occupancy-s$heterogeneity)<1e-8));write.csv(s,file.path(out,sprintf('conditional_%s_k%d.csv',pool,k)),row.names=FALSE);# export deletion-fit ranges, not confidence intervals
   del<-s |> filter(fit_id!='all') |> pivot_longer(all_of(COMP),names_to='component',values_to='pred') |> group_by(group,fit_id,trophic,component) |> arrange(spread,.by_group=TRUE) |> mutate(change_pp=100*(pred-first(pred))) |> ungroup() |> group_by(group,trophic,component,spread) |> summarise(n_fits=n_distinct(fit_id),lo=if(n_fits>=3)min(change_pp)else NA_real_,hi=if(n_fits>=3)max(change_pp)else NA_real_,.groups='drop');write.csv(del,file.path(out,sprintf('conditional_%s_k%d_deletion_ranges.csv',pool,k)),row.names=FALSE)}}};write.csv(bind_rows(DIAG),file.path(out,'conditional_support_summary.csv'),row.names=FALSE)
# Do not interpret failure to establish overlap as zero response. If no common
# support, retain the explicitly labelled binned descriptive figure instead.
