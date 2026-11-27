# Figure 6: response shape and compensation between coherent components.
# Rscript figure6_closing.R figures/figure6_redesign figures/figure6_closing [conditional_dir]
# Run in docs/final_empirical_analysis. Without conditional_dir, plots the
# existing supported cell means (NOT isolated TP effects). With conditional_dir,
# reads the matched-profile predictions from figure6_conditional.R instead.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(grid); library(gtable); library(scales)})
stopifnot(requireNamespace('svglite',quietly=TRUE))
a <- commandArgs(TRUE); input <- if(length(a)>=1) a[1] else 'figures/figure6_redesign'; out <- if(length(a)>=2) a[2] else 'figures/figure6_closing'; cond <- if(length(a)>=3) a[3] else NULL
dir.create(out,recursive=TRUE,showWarnings=FALSE)
TX <- c('BM','ZP'); TXL <- c(BM='EPA macroinvertebrates',ZP='EPA zooplankton')
COMP <- c('total','occupancy','heterogeneity'); CL <- c(total='Total dissimilarity',occupancy='Occupancy-associated',heterogeneity='Abundance heterogeneity')
POS <- c('Oligotrophic','Mesotrophic','Eutrophic'); COL <- setNames(c('#79936C','#246878','#C45C19'),POS); SH <- setNames(c(16,17,15),POS)
ink <- '#252B2D'; gridcol <- '#E5E9E9'; CC <- c(occupancy='#246878',heterogeneity='#C45C19')
base <- theme_minimal(base_size=9,base_family='sans')+theme(text=element_text(colour=ink),axis.title=element_blank(),axis.text=element_text(size=8,colour=ink),axis.ticks=element_line(linewidth=.25,colour=ink),axis.ticks.length=unit(1,'mm'),axis.line=element_line(linewidth=.3,colour=ink),panel.grid.major=element_line(colour=gridcol,linewidth=.25),panel.grid.minor=element_blank(),plot.margin=margin(1.2,1.5,1.2,1.5),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),legend.key.width=unit(7,'mm'),legend.margin=margin(0,0,0,0))
zero <- geom_hline(yintercept=0,colour='#9DA5A7',linewidth=.4)
read_input <- function(pool,k) {
 if(!is.null(cond)) {
   p <- file.path(cond,sprintf('conditional_%s_k%d.csv',pool,k)); if(!file.exists(p)) stop('Missing ',p)
   x <- read.csv(p); x <- x |> filter(model=='additive',fit_id=='all',supported)
   if(!nrow(x)) stop('No supported conditional predictions: ',pool,' k=',k,'; inspect diagnostics, do not silently substitute cell means.')
   x <- x |> pivot_longer(all_of(COMP),names_to='component',values_to='pred') |> mutate(segment=1L)
 } else {
   stem <- if(k==3) if(pool=='shared') 'Figure6_main_redesign' else 'FigureS6_full_pool_redesign' else paste0('FigureS6_TP_B5_',pool)
   p <- file.path(input,paste0(stem,'_B_supported_predictions.csv')); if(!file.exists(p)) stop('Missing ',p)
   x <- read.csv(p) |> filter(model=='additive',order==k)
 }
 stopifnot(all(c('group','order','component','trophic','TP_position','spread','pred','segment') %in% names(x)))
 stopifnot(all(is.finite(x$pred)),all(x$pred>=-1e-8 & x$pred<=1+1e-8),all(x$group %in% TX),all(x$component %in% COMP),all(x$trophic %in% POS))
 if(anyDuplicated(x[c('group','component','trophic','spread')])) stop('Duplicate prediction cells')
 # Each taxon uses exactly the same spread grid across all positions/components.
 x <- bind_rows(lapply(TX,function(tx){ z <- x |> filter(group==tx); if(!nrow(z)) stop('Missing taxon ',tx)
   if(!setequal(unique(z$trophic),POS)||!setequal(unique(z$component),COMP)) stop('Incomplete positions/components for ',pool,' ',tx,' k=',k)
   grids <- split(round(z$spread,8),interaction(z$trophic,z$component,drop=TRUE)); common <- Reduce(intersect,lapply(grids,unique))
   if(length(common)<3) stop('Fewer than three common supported spread cells: ',pool,' ',tx,' k=',k)
   z <- z |> filter(round(spread,8)%in%common)
   # Preserve gaps rather than draw a smooth or interpolate across missing cells.
   z <- z |> group_by(component,trophic) |> arrange(spread,.by_group=TRUE) |> mutate(segment=as.integer(segment)) |> ungroup(); z
 }))
 w <- x |> select(group,trophic,spread,component,pred) |> pivot_wider(names_from=component,values_from=pred)
 stopifnot(all(abs(w$total-w$occupancy-w$heterogeneity)<1e-7))
 x <- x |> group_by(group,component,trophic) |> arrange(spread,.by_group=TRUE) |> mutate(delta_pp=100*(pred-first(pred))) |> ungroup()
 end <- x |> group_by(group,component,trophic) |> arrange(spread,.by_group=TRUE) |> summarise(TP_position=first(TP_position),spread_from=first(spread),spread_to=last(spread),change_pp=last(delta_pp),n_bins=n(),.groups='drop')
 ww <- end |> select(group,trophic,component,change_pp) |> pivot_wider(names_from=component,values_from=change_pp)
 stopifnot(all(abs(ww$total-ww$occupancy-ww$heterogeneity)<1e-5))
 x$trophic <- factor(x$trophic,POS); end$trophic <- factor(end$trophic,POS)
 list(curve=x,end=end)
}
# Pre-specified order 3 primary, 5 sensitivity. Same scale for each taxon
# across components and pools, but independent across taxa, as in Figure 5.
D <- setNames(lapply(c(3,5),function(k) setNames(lapply(c('shared','full'),function(pool){
 if(is.null(cond)) return(read_input(pool,k))
 tryCatch(read_input(pool,k),error=function(e){warning(conditionMessage(e));NULL})
}),c('shared','full'))),c('3','5'))
if(!any(vapply(unlist(D,recursive=FALSE),function(x)!is.null(x),logical(1)))) stop('No complete supported figure. Inspect conditional_support_summary.csv; use explicitly labelled binned figures instead.')
scale_for <- function(k,tx){ v<-unlist(lapply(D[[as.character(k)]],function(d) if(is.null(d)) numeric() else d$curve$delta_pp[d$curve$group==tx])); m<-max(abs(v)); if(m==0) m<-1; r<-c(-1,1)*m*1.10; list(lim=r,br=pretty(r,n=4)[pretty(r,n=4)>=r[1]&pretty(r,n=4)<=r[2]]) }
xs_for <- function(k,tx){ v<-unlist(lapply(D[[as.character(k)]],function(d)if(is.null(d))numeric()else d$curve$spread[d$curve$group==tx])); range(v) }
ys <- function(s) scale_y_continuous(limits=s$lim,breaks=s$br,expand=expansion(mult=0),labels=label_number())
make_curve <- function(d,k,tx,cp,showx){ z<-d$curve |> filter(group==tx,component==cp); s<-scale_for(k,tx)
 p<-ggplot(z,aes(spread,delta_pp,colour=trophic,shape=trophic,group=interaction(trophic,segment)))+zero+geom_line(linewidth=.8)+geom_point(size=1.65)+scale_colour_manual(values=COL,limits=POS,labels=c('Low (~10)','Intermediate (~35)','High (~100)'),drop=FALSE)+scale_shape_manual(values=SH,limits=POS,labels=c('Low (~10)','Intermediate (~35)','High (~100)'),drop=FALSE)+scale_x_continuous(limits=xs_for(k,tx),breaks=breaks_extended(4),expand=expansion(mult=.045))+ys(s)+base+theme(legend.position='none')
 if(!showx) p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank());p }
make_bar <- function(d,k,tx){ z<-d$end |> filter(group==tx); s<-scale_for(k,tx); b<-z |> filter(component!='total'); t<-z |> filter(component=='total');
 ggplot()+zero+geom_col(data=b,aes(trophic,change_pp,fill=component),position='stack',width=.53)+geom_point(data=t,aes(trophic,change_pp),size=2.5,shape=21,fill=ink,colour='white',stroke=.4)+scale_fill_manual(values=CC,limits=names(CC),labels=c('Occupancy-associated','Abundance heterogeneity'))+scale_x_discrete(limits=POS,labels=c('Low','Intermediate','High'))+ys(s)+base+theme(legend.position='none',panel.grid.major.x=element_blank()) }
# Fixed gtable: every row shares the same two column widths and panel edges.
# Equalize ALL gtable width tracks across plots; single-panel plots only.
label <- function(txt,size=8.5,bold=FALSE,rot=0) textGrob(txt,rot=rot,gp=gpar(fontsize=size,fontface=if(bold)'bold'else'plain',col=ink))
legend_grob <- function(kind){ if(kind=='pos') {labs<-c('Low (~10)','Intermediate (~35)','High (~100)'); cols<-unname(COL); sh<-unname(SH); x<-c(.20,.50,.80); children<-list(); for(i in seq_along(x)){children<-c(children,list(segmentsGrob(unit(x[i]-.095,'npc'),unit(.5,'npc'),unit(x[i]-.045,'npc'),unit(.5,'npc'),gp=gpar(col=cols[i],lwd=1.5)),pointsGrob(unit(x[i]-.070,'npc'),unit(.5,'npc'),pch=sh[i],size=unit(1.6,'mm'),gp=gpar(col=cols[i])),textGrob(labs[i],x=unit(x[i]-.035,'npc'),hjust=0,gp=gpar(fontsize=8,col=ink))))}; do.call(grobTree,children) }
 else {grobTree(rectGrob(x=.10,y=.5,width=.04,height=.38,gp=gpar(fill=CC[1],col=NA)),textGrob('Occupancy-associated',x=.13,hjust=0,gp=gpar(fontsize=8,col=ink)),rectGrob(x=.45,y=.5,width=.04,height=.38,gp=gpar(fill=CC[2],col=NA)),textGrob('Abundance heterogeneity',x=.48,hjust=0,gp=gpar(fontsize=8,col=ink)),pointsGrob(.84,.5,pch=16,size=unit(1.8,'mm'),gp=gpar(col=ink)),textGrob('Total',x=.865,hjust=0,gp=gpar(fontsize=8,col=ink)))}}
build <- function(d,k){ ps<-lapply(COMP,function(cp)lapply(TX,function(tx)make_curve(d,k,tx,cp,cp=='heterogeneity'))); ps[[4]]<-lapply(TX,function(tx)make_bar(d,k,tx)); gs<-lapply(ps,function(row)lapply(row,ggplotGrob)); flat<-unlist(gs,recursive=FALSE); # flatten one level only
 nwidth<-vapply(flat,function(g)length(g$widths),integer(1));stopifnot(length(unique(nwidth))==1L)
 ww<-do.call(unit.pmax,lapply(flat,function(g)g$widths)); for(i in seq_along(gs))for(j in 1:2)gs[[i]][[j]]$widths<-ww
 # Matching heights within each paired row; external text and legends have fixed rows.
 for(i in seq_along(gs)){hh<-unit.pmax(gs[[i]][[1]]$heights,gs[[i]][[2]]$heights); for(j in 1:2)gs[[i]][[j]]$heights<-hh}
 h<-unit.c(unit(7,'mm'),unit(1,'null'),unit(1,'null'),unit(1,'null'),unit(6,'mm'),unit(6,'mm'),unit(1.10,'null'),unit(6,'mm'),unit(7,'mm'))
 g<-gtable(widths=unit.c(unit(5,'mm'),unit(9,'mm'),unit(1,'null'),unit(4,'mm'),unit(1,'null')),heights=h)
 add<-function(g,o,t,l,b=t,r=l)gtable_add_grob(g,o,t=t,l=l,b=b,r=r,clip='off')
 for(j in 1:2){c0<-2*j+1;g<-add(g,label(TXL[TX[j]],10,TRUE),1,c0);for(i in 1:3)g<-add(g,gs[[i]][[j]],i+1,c0);g<-add(g,gs[[4]][[j]],7,c0)}
 g<-add(g,label('A',12,TRUE),2,1);g<-add(g,label('B',12,TRUE),7,1)
 for(i in 1:3)g<-add(g,label(paste0(CL[COMP[i]],'\nChange (pp)'),8,TRUE,90),i+1,2)
 g<-add(g,label('Endpoint change (pp)',8,TRUE,90),7,2)
 g<-add(g,label(paste0('TP spread: mean pairwise |difference in log10 TP|  (k = ',k,')'),9),5,3,r=5)
 g<-add(g,legend_grob('pos'),6,3,r=5);g<-add(g,label('TP position',9),8,3,r=5);g<-add(g,legend_grob('comp'),9,3,r=5);g }
export <- function(g,stem){ for(ext in c('pdf','svg','png'))ggsave(file.path(out,paste0(stem,'.',ext)),g,width=183,height=195,units='mm',dpi=600,bg='white',device=switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png')) }
for(k in c(3,5))for(pool in c('shared','full')){ d<-D[[as.character(k)]][[pool]]; if(is.null(d))next; stem<-paste0(if(is.null(cond))'Figure6_closing_binned'else'Figure6_closing_conditional','_',pool,'_k',k);export(build(d,k),stem); write.csv(d$curve,file.path(out,paste0(stem,'_A_curves.csv')),row.names=FALSE);write.csv(d$end,file.path(out,paste0(stem,'_B_endpoints.csv')),row.names=FALSE) }
# Caption: all curves are changes in predicted dissimilarity (percentage points)
# from the lowest common supported spread, separately for each TP position.
# Exact TP positions and endpoints are exported. B bars are the signed endpoint
# changes; the point is their sum. These are not independent contributions to D2.
# Binned mode: other predictors follow their distribution among observed cells.
# Conditional mode: other covariates and within-set TP shape held fixed per profile;
# common support screened before averaging. No causal inference.
