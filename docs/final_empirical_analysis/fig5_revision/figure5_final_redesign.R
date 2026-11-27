# Figure 5 | Added predictive value and order-dependent environmental sensitivity
# From docs/final_empirical_analysis:
# Rscript fig5_revision/figure5_final_redesign.R fig5_revision/tables fig5_revision/figures
# Plotting only: reads committed 50-52 outputs; no refits or new predictions.
# Keeps candidate and existing supplementary figures intact.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(scales)
})
stopifnot(requireNamespace('svglite', quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else 'fig5_revision/tables'
out_dir <- if (length(args) >= 2) args[2] else 'fig5_revision/figures'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink <- '#252B2D'; petrol <- '#246878'; orange <- '#C45C19'
purple <- '#615A83'; grey <- '#848B8D'; grid_col <- '#E5E9E9'
W <- 183; H <- 212
panels <- c('BM_inc','BM_rel','ZP_inc','ZP_rel')
labels_rep <- c(BM_inc='Incidence', BM_rel='Relative abundance',
                ZP_inc='Incidence', ZP_rel='Relative abundance')
contrast_lab <- c(delta_linear_vs_spread='Linear position',
                  delta_nonlinear_vs_spread='Nonlinear position')
contrast_col <- c('Linear position'=petrol, 'Nonlinear position'=orange)
pred_lab <- c(TP='TP', TN='TN', Chla='Chl-a', distance='Geographic distance')
pred_col <- c(TP=petrol, TN=orange, 'Chl-a'=purple, 'Geographic distance'=grey)
tp_col <- c('TP spread'=petrol, 'TP position'=orange)
tp_type <- c('TP spread'='solid', 'TP position'='42')

base <- theme_minimal(base_size=9, base_family='sans') + theme(
  text=element_text(colour=ink),
  axis.text=element_text(size=7.5, colour=ink),
  axis.title=element_text(size=8.2, colour=ink),
  axis.title.y=element_text(margin=margin(r=4)),
  axis.title.x=element_text(margin=margin(t=3)),
  axis.line=element_line(colour=ink, linewidth=.3),
  axis.ticks=element_line(colour=ink, linewidth=.25),
  axis.ticks.length=grid::unit(1,'mm'),
  panel.grid.major=element_line(colour=grid_col, linewidth=.25),
  panel.grid.minor=element_blank(),
  panel.spacing.x=grid::unit(3.5,'mm'),
  strip.background=element_blank(),
  strip.text=element_text(size=8, face='bold', colour=ink, margin=margin(b=4)),
  legend.position='bottom', legend.title=element_blank(),
  legend.text=element_text(size=7.2, colour=ink),
  legend.key.width=grid::unit(8,'mm'), legend.key.height=grid::unit(2.8,'mm'),
  legend.spacing.x=grid::unit(1.5,'mm'), legend.margin=margin(0,0,0,0),
  legend.box.spacing=grid::unit(1,'mm'),
  plot.tag=element_text(size=11, face='bold', colour=ink),
  plot.margin=margin(2,3,1,2))
no_strip <- theme(strip.text=element_blank())
no_x <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
zero <- geom_hline(yintercept=0, colour='#A2AAAC', linewidth=.35)

read_tab <- function(stem,pool) {
  p <- file.path(in_dir,paste0(stem,'_',pool,'.csv'))
  if (!file.exists(p)) stop('Missing committed table: ',p)
  d <- read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(all(d$pool==pool)); d
}
add_panel <- function(d) {
  d <- d |> mutate(panel=factor(case_when(
    community=='EPA macroinvertebrates' & response=='incidence' ~ 'BM_inc',
    community=='EPA macroinvertebrates' & response=='relative' ~ 'BM_rel',
    community=='EPA zooplankton' & response=='incidence' ~ 'ZP_inc',
    community=='EPA zooplankton' & response=='relative' ~ 'ZP_rel',
    TRUE ~ NA_character_),levels=panels))
  stopifnot(!anyNA(d$panel)); d
}
facets <- function() facet_grid(cols=vars(panel),scales='free_x',space='fixed',
  drop=FALSE,labeller=labeller(panel=as_labeller(labels_rep)))
order_breaks <- function(lim) if(max(lim)<=6.5) 2:5 else c(2,4,6,8,10,12)
order_scale <- function() scale_x_continuous(breaks=order_breaks,
  expand=expansion(mult=c(.05,.05)))
# Header is separate from the representation strips; one label per taxon group.
header <- ggplot(data.frame(x=c(.25,.75),
  label=c('EPA macroinvertebrates','EPA zooplankton')),aes(x,.5,label=label))+
  geom_text(size=3.35,fontface='bold',colour=ink)+
  scale_x_continuous(limits=c(0,1),expand=c(0,0))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))+
  theme_void()+theme(plot.margin=margin(0,3,0,2))

prepare <- function(pool) {
  F <- read_tab('sens_paired_folds',pool) |> filter(transform=='log10',contrast %in% names(contrast_lab)) |> add_panel()
  S <- read_tab('sens_paired_summary',pool) |> filter(transform=='log10',contrast %in% names(contrast_lab)) |> add_panel()
  I <- read_tab('sens_importance',pool) |> filter(transform=='log10',model=='nonlinear_const') |> add_panel()
  C <- read_tab('sens_curves_summary',pool) |> filter(transform=='log10',position_model=='nonlinear',predictor=='TP') |> add_panel()
  # Check plotted values and coverage. Never silently fill missing folds/orders.
  stopifnot(nrow(F)>0,nrow(S)>0,nrow(I)>0,nrow(C)>0,
    !anyDuplicated(F[c('panel','order','fold','contrast')]),
    !anyDuplicated(S[c('panel','order','contrast')]),
    !anyDuplicated(I[c('panel','order','predictor','component')]),
    !anyDuplicated(C[c('panel','order','axis')]),
    all(S$n_folds==5), all(I$n_folds==5),
    all(is.finite(F$delta_D2)),all(is.finite(S$median)),all(is.finite(I$abs_mean)),
    all(I$predictor %in% names(pred_lab)),all(I$component %in% c('spread','position','distance')),
    all(C$axis %in% c('spread','position')),all(C$n_feasible==25),
    all(is.finite(C$link_change)),all(is.finite(C$response_change)),
    all(C$x_to>C$x_from),all(C$mean_dissimilarity>=0 & C$mean_dissimilarity<=1),
    all(abs(C$link_change-(C$link_to-C$link_from))<1e-10),
    all(abs(C$response_change-(C$response_to-C$response_from))<1e-10),
    all(abs(C$response_from-plogis(C$link_from))<1e-10),
    all(abs(C$response_to-plogis(C$link_to))<1e-10))
  # The endpoints must be identical across orders AND representations within
  # each community and axis. They differ between taxa and slightly between pools.
  ranges <- C |> group_by(community,axis) |> summarise(
    x_from=min(x_from),x_to=max(x_to),
    range_lo_diff=max(x_from)-min(x_from),range_hi_diff=max(x_to)-min(x_to),
    .groups='drop')
  stopifnot(all(ranges$range_lo_diff<1e-10),all(ranges$range_hi_diff<1e-10))
  ck <- S |> select(panel,order) |> distinct()
  stopifnot(nrow(anti_join(ck,distinct(C,panel,order),by=c('panel','order')))==0,
            nrow(anti_join(distinct(C,panel,order),ck,by=c('panel','order')))==0)
  chk <- F |> group_by(panel,order,contrast) |> summarise(m=median(delta_D2),n=n(),.groups='drop') |>
    inner_join(S,by=c('panel','order','contrast'))
  stopifnot(nrow(chk)==nrow(S),all(chk$n==5),all(abs(chk$m-chk$median)<1e-10))
  F <- F |> mutate(position=factor(unname(contrast_lab[contrast]),levels=unname(contrast_lab)))
  S <- S |> mutate(position=factor(unname(contrast_lab[contrast]),levels=unname(contrast_lab)))
  I <- I |> mutate(pred=factor(unname(pred_lab[predictor]),levels=unname(pred_lab)),
    term=factor(ifelse(component=='position','Position','Spread or distance'),levels=c('Spread or distance','Position')))
  C <- C |> mutate(tp_axis=factor(ifelse(axis=='spread','TP spread','TP position'),levels=names(tp_col)))
  list(F=F,S=S,I=I,C=C,ranges=mutate(ranges,pool=pool))
}
D <- list(shared=prepare('shared'),full=prepare('full'))
# One scale per row, identical across both pools. All extrema are retained.
limits <- function(x) {
  x <- x[is.finite(x)]; stopifnot(length(x)>0)
  r <- range(c(0,x)); if(diff(r)==0) r <- r+c(-.01,.01); r
}
med_rng <- range(unlist(lapply(D,function(d)d$S$median)))
pad <- .15*diff(med_rng)
L <- list(A=limits(unlist(lapply(D,function(d)d$F$delta_D2))),
          Am=c(floor((min(med_rng,0)-pad)/.05)*.05,ceiling((max(med_rng,0)+pad)/.05)*.05),
          B=limits(unlist(lapply(D,function(d)d$I$abs_mean))),
          C=limits(unlist(lapply(D,function(d)d$C$link_change))),
          D=limits(unlist(lapply(D,function(d)d$C$response_change))))
y_scale <- function(lim,labels=waiver()) scale_y_continuous(limits=lim,
  breaks=breaks_extended(4),labels=labels,expand=expansion(mult=c(.045,.065)))

# A: paired held-out improvement. Transparent individual fold points (not CIs),
# slight horizontal dodge purely for legibility. Median lines at true orders.
make_A <- function(d, folds=FALSE, show_x=FALSE) {
  # Main figure: medians only. Supplement (folds=TRUE): all fold values, nothing clipped.
  # Each of the four columns (taxon x representation) has its OWN y range, from its own
  # values (medians; plus fold values in the supplement), zero always included. To keep the
  # panel geometry identical to B-D, values are mapped linearly to [0,1] within each column
  # (layout only; no value is altered) and gridlines/tick labels are drawn at the true values
  # inside the panels.
  rng <- lapply(panels,function(pn){
    v <- d$S$median[d$S$panel==pn]; if(folds) v <- c(v,d$F$delta_D2[d$F$panel==pn])
    r <- range(c(0,v)); pd <- (if(folds) .06 else .12)*diff(r); r+c(-pd,pd) })
  names(rng) <- panels
  tf <- function(y,pn) { lo <- vapply(as.character(pn),function(q)rng[[q]][1],0); hi <- vapply(as.character(pn),function(q)rng[[q]][2],0); (y-lo)/(hi-lo) }
  grd <- do.call(rbind,lapply(panels,function(pn){ r <- rng[[pn]]; b <- breaks_extended(5)(r); sp <- diff(r); b <- b[b>=r[1]+.03*sp & b<=r[2]-.10*sp]
    data.frame(panel=factor(pn,levels=panels),v=b,y=tf(b,pn),label=label_number(accuracy=.01,trim=TRUE)(b)) }))
  Fd <- d$F; Fd$y <- tf(Fd$delta_D2,Fd$panel); Sd <- d$S; Sd$y <- tf(Sd$median,Sd$panel)
  zr <- data.frame(panel=factor(panels,levels=panels),y=tf(rep(0,4),panels))
  p <- ggplot()+geom_hline(data=grd,aes(yintercept=y),colour=grid_col,linewidth=.25)+
    geom_hline(data=zr,aes(yintercept=y),colour='#A2AAAC',linewidth=.35)
  if(folds) p <- p+geom_point(data=Fd,aes(order,y,colour=position,shape=position),
      position=position_dodge(width=.23),size=1.05,alpha=.45,show.legend=FALSE)
  p <- p+
    geom_line(data=Sd,aes(order,y,colour=position,linetype=position),linewidth=.85)+
    geom_point(data=Sd,aes(order,y,colour=position,shape=position),size=1.65)+
    geom_label(data=grd,aes(x=-Inf,y=y,label=label),hjust=-.1,vjust=-.2,size=2.5,colour=ink,fill=NA,linewidth=0,label.padding=grid::unit(.5,'mm'))+
    facets()+order_scale()+scale_y_continuous(limits=c(0,1),breaks=NULL,expand=expansion(mult=0))+
    scale_colour_manual(values=contrast_col)+
    scale_shape_manual(values=c('Linear position'=16,'Nonlinear position'=17))+
    scale_linetype_manual(values=c('Linear position'='solid','Nonlinear position'='42'))+
    guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+
    labs(tag='A',y=expression(atop('Paired held-out '*Delta*D^2,'vs spread only')),x=if(show_x) 'Zeta order, k' else NULL,
      colour=NULL,shape=NULL,linetype=NULL)+base+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.ticks.length.y=grid::unit(0,'mm'),panel.grid.major.y=element_blank())
  if(show_x) p else p+no_x
}
# B: absolute, signed Shapley/Owen allocations from the nonlinear position fit.
make_B <- function(d) {
  ggplot(d$I,aes(order,abs_mean,colour=pred,linetype=term,group=interaction(pred,term)))+zero+
    geom_line(linewidth=.65)+geom_point(size=1.05,show.legend=FALSE)+
    facets()+order_scale()+y_scale(L$B)+scale_colour_manual(values=pred_col)+
    scale_linetype_manual(values=c('Spread or distance'='solid',Position='42'))+
    guides(colour=guide_legend(order=1,nrow=1,override.aes=list(linetype='solid')),
      linetype=guide_legend(order=2,nrow=1,override.aes=list(colour=ink,linetype=c('solid','42'))))+
    labs(tag='B',y=expression('Absolute contribution ('*Delta*D^2*')'),x=NULL,colour=NULL,linetype=NULL)+
    base+no_strip+no_x+theme(legend.box='vertical')
}
# C/D: same fitted conditional contrast, two scales. Unlike the candidate
# curves, the order-specific intercept does not obscure the contrast magnitude.
# D expressed in percentage points, not relative percent change.
make_effect <- function(d,link=TRUE) {
  v <- if(link) 'link_change' else 'response_change'
  yy <- if(link) 'Conditional change in\nlogit dissimilarity' else 'Conditional change in\ndissimilarity (pp)'
  fmt <- if(link) label_number(accuracy=.01,trim=TRUE) else function(x) label_number(accuracy=.1,trim=TRUE)(100*x)
  p <- ggplot(d$C,aes(order,.data[[v]],colour=tp_axis,linetype=tp_axis,shape=tp_axis,group=tp_axis))+zero+
    geom_line(linewidth=.75)+geom_point(size=1.65)+facets()+order_scale()+
    y_scale(if(link)L$C else L$D,labels=fmt)+scale_colour_manual(values=tp_col)+
    scale_linetype_manual(values=tp_type)+scale_shape_manual(values=c('TP spread'=16,'TP position'=17))+
    guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1))+
    labs(tag=if(link)'C' else 'D',y=yy,x='Zeta order, k',colour=NULL,linetype=NULL,shape=NULL)+base+no_strip
  if(link) p+no_x+theme(legend.position='none') else p
}

# Supplemental ceiling diagnostic (not a new analysis). Both TP axes have
# separately computed mean fitted dissimilarity over their common grids.
make_ceiling <- function(d) {
  ggplot(d$C,aes(order,mean_dissimilarity,colour=tp_axis,linetype=tp_axis,shape=tp_axis,group=tp_axis))+
    geom_hline(yintercept=1,colour=grey,linetype='dotted',linewidth=.5)+
    geom_line(linewidth=.7)+geom_point(size=1.5)+facets()+order_scale()+
    scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=c(.01,.04)))+
    scale_colour_manual(values=tp_col)+scale_linetype_manual(values=tp_type)+
    scale_shape_manual(values=c('TP spread'=16,'TP position'=17))+
    labs(y='Mean fitted\ndissimilarity',x='Zeta order, k',colour=NULL,linetype=NULL,shape=NULL)+base
}
stack <- function(plots,hs) wrap_plots(c(list(header),plots),ncol=1,
  heights=grid::unit(c(6,hs),c('mm',rep('null',length(hs)))))
export <- function(p,name,h) {
  for(ext in c('pdf','svg','png')) {
    device <- switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png')
    ggsave(file.path(out_dir,paste0(name,'.',ext)),p,width=W,height=h,units='mm',dpi=600,device=device,bg='white')
  }
}
dump <- function(d,name,suffix) write.csv(d,file.path(out_dir,paste0(name,'_',suffix,'.csv')),row.names=FALSE)
for(pool in names(D)) {
  d <- D[[pool]]
  nm <- if(pool=='shared') 'Figure5_final' else 'FigureS5_final_full_pool'
  fig <- stack(list(make_A(d),make_B(d),make_effect(d,TRUE),make_effect(d,FALSE)),c(1.15,1.20,1,1.08))
  export(fig,nm,H)
  dump(d$F,nm,'A_folds'); dump(d$S,nm,'A_medians'); dump(d$I,nm,'B_importance')
  dump(mutate(d$C,response_change_pp=100*response_change),nm,'CD_TP_contrasts')
  dump(d$ranges,nm,'TP_common_ranges')
  sn <- paste0('FigureS5_final_panelA_folds_',pool)
  export(stack(list(make_A(d,folds=TRUE,show_x=TRUE)),1),sn,80)
  dump(d$F,sn,'fold_values'); dump(d$S,sn,'medians')
  cn <- paste0('FigureS5_final_ceiling_',pool)
  export(stack(list(make_ceiling(d)),1),cn,76)
  dump(d$C,cn,'values')
  message('Saved ',nm,' and ',cn)
}
# Caption notes (for manuscript, not printed on the figure):
# A: small points are individual held-out spatial-fold contrasts, after seed
# averaging, slightly dodged; larger points/lines are medians across five folds.
# Both models are compared to the existing spread-only baseline (different
# fitting procedure); this is a predictive comparison, not a significance test.
# B: mean signed held-out Shapley/Owen absolute allocations for the log10
# nonlinear position model, averaged across seeds within fold and then folds.
# C,D: difference between fitted endpoint predictions, using the same 25-point
# feasible common range within each taxon, representation and TP axis.
# The ranges differ among taxa and slightly between pools; see exported ranges.
# Other design terms held at their medians per order; focal design terms are
# locally averaged as in 51_conditional_curves.R. These are descriptive
# conditional contrasts, not causal effects; fold-level curve uncertainty is
# unavailable. Do not treat missing error bars as evidence of significance.
# D uses percentage points (100 * difference in predicted dissimilarity).
