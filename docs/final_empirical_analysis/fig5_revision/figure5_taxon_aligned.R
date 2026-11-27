# Figure 5 | Added predictive value and order-dependent environmental sensitivity
# From docs/final_empirical_analysis:
# Rscript fig5_revision/figure5_taxon_aligned.R fig5_revision/tables fig5_revision/figures [both|pool]
# Plotting only: reads committed 50-52 outputs; no refits or new predictions.
# Keeps candidate and existing supplementary figures intact.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(grid); library(gtable)
})
stopifnot(requireNamespace('svglite', quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else 'fig5_revision/tables'
out_dir <- if (length(args) >= 2) args[2] else 'fig5_revision/figures'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ink <- '#252B2D'; petrol <- '#246878'; orange <- '#C45C19'
purple <- '#615A83'; grey <- '#848B8D'; grid_col <- '#E5E9E9'
W <- 183; H <- 212
axis_scope <- if (length(args) >= 3) args[3] else "both"
stopifnot(axis_scope %in% c("both", "pool"))
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
no_strip <- theme(strip.text=element_blank(), strip.text.x=element_blank())
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
# ---- Explicit scales and independent taxon plots ----------------------------
# Each taxon receives a two-column facet_grid with FIXED x and y scales.
# No normalization, synthetic axis text, free scales, or nested patchwork.
# Default: same taxon-specific limits in shared and full-pool figures.
# Optional third CLI arg 'pool': each pool has its own limits as well.
taxa <- c('BM','ZP')
taxon_name <- c(BM='EPA macroinvertebrates',ZP='EPA zooplankton')
subset_tx <- function(x, tx) {
  x <- x[substr(as.character(x$panel),1,2)==tx,,drop=FALSE]
  x$panel <- factor(as.character(x$panel), levels=panels[substr(panels,1,2)==tx])
  stopifnot(nrow(x)>0, nlevels(x$panel)==2, all(table(x$panel)>0))
  x
}
metric_values <- function(d, tx, metric) {
  switch(metric,
    A = subset_tx(d$S,tx)$median,
    Af = c(subset_tx(d$S,tx)$median,subset_tx(d$F,tx)$delta_D2),
    B = subset_tx(d$I,tx)$abs_mean,
    C = subset_tx(d$C,tx)$link_change,
    D = 100*subset_tx(d$C,tx)$response_change,
    K = subset_tx(d$C,tx)$mean_dissimilarity,
    stop('Unknown metric: ',metric))
}
scale_spec <- function(pool, tx, metric) {
  if (metric=='K') return(list(lim=c(0,1), breaks=c(0,.25,.5,.75,1)))
  ps <- if(axis_scope=='both') names(D) else pool
  v <- unlist(lapply(ps,function(p) metric_values(D[[p]],tx,metric)),use.names=FALSE)
  stopifnot(length(v)>0,all(is.finite(v)))
  r <- range(c(0,v)); if(diff(r)==0) r <- r+c(-.01,.01)
  pad <- .055*diff(r)
  lim <- r+c(-pad,pad)
  # Keep true values, and always include zero in the visible range and breaks.
  b <- sort(unique(c(0,scales::breaks_extended(n=4)(lim))))
  b <- b[b>=lim[1]-1e-12 & b<=lim[2]+1e-12]
  list(lim=lim,breaks=b)
}
scale_y_tx <- function(pool,tx,metric) {
  sp <- scale_spec(pool,tx,metric)
  scale_y_continuous(limits=sp$lim,breaks=sp$breaks,
    labels=scales::label_number(accuracy=NULL,trim=TRUE),
    expand=if(metric=='K') expansion(mult=c(.015,.045)) else expansion(mult=0))
}
scale_x_tx <- function(tx) {
  r <- if(tx=='BM') c(2,12) else c(2,5)
  b <- if(tx=='BM') c(2,4,6,8,10,12) else 2:5
  scale_x_continuous(limits=r,breaks=b,expand=expansion(mult=c(.06,.06)))
}
facet_tx <- function() facet_grid(cols=vars(panel),scales='fixed',space='fixed',
  drop=TRUE,labeller=labeller(panel=as_labeller(labels_rep)))
# Both taxon plots use identical gtable structure. Axis titles and tags are
# separate grobs in a common left gutter, not hidden titles of right-hand plots.
style_tx <- base+theme(
  axis.title=element_blank(),plot.margin=margin(1,0,1,0),
  panel.spacing.x=unit(3,'mm'),
  legend.position='bottom',legend.box.spacing=unit(0,'mm'),
  strip.text.x=element_text(size=8,face='bold',margin=margin(b=3)))
model_shape <- c('Linear position'=16,'Nonlinear position'=17)
model_type <- c('Linear position'='solid','Nonlinear position'='42')
tp_shape <- c('TP spread'=16,'TP position'=17)

make_plot <- function(pool,tx,metric,show_strip=FALSE,show_x=FALSE) {
  d <- D[[pool]]
  if (metric %in% c('A','Af')) {
    f <- subset_tx(d$F,tx); s <- subset_tx(d$S,tx)
    p <- ggplot()+zero
    if(metric=='Af') p <- p+geom_point(data=f,
      aes(order,delta_D2,colour=position,shape=position),
      size=1,alpha=.42,show.legend=FALSE)
    p <- p+geom_line(data=s,aes(order,median,colour=position,linetype=position,
      group=position),linewidth=.78)+
      geom_point(data=s,aes(order,median,colour=position,shape=position),size=1.5)+
      scale_colour_manual(values=contrast_col,limits=names(contrast_col),drop=FALSE)+
      scale_linetype_manual(values=model_type,limits=names(contrast_col),drop=FALSE)+
      scale_shape_manual(values=model_shape,limits=names(contrast_col),drop=FALSE)+
      guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1))
  } else if(metric=='B') {
    i <- subset_tx(d$I,tx)
    p <- ggplot(i,aes(order,abs_mean,colour=pred,linetype=term,
      group=interaction(pred,term)))+zero+
      geom_line(linewidth=.65)+geom_point(size=1.05,show.legend=FALSE)+
      scale_colour_manual(values=pred_col,limits=names(pred_col),drop=FALSE)+
      scale_linetype_manual(values=c('Spread or distance'='solid',Position='42'))+
      guides(colour=guide_legend(order=1,nrow=1,override.aes=list(linetype='solid')),
        linetype=guide_legend(order=2,nrow=1,override.aes=list(colour=ink,linetype=c('solid','42'))))+
      theme(legend.box='vertical')
  } else {
    cc <- subset_tx(d$C,tx)
    # Panel D uses actual percentage-point units; all source values preserved.
    cc$plot_y <- switch(metric,C=cc$link_change,D=100*cc$response_change,K=cc$mean_dissimilarity)
    p <- ggplot(cc,aes(order,plot_y,colour=tp_axis,linetype=tp_axis,shape=tp_axis,group=tp_axis))+
      (if(metric=='K') geom_hline(yintercept=1,colour=grey,linetype='dotted',linewidth=.5) else zero)+
      geom_line(linewidth=.75)+geom_point(size=1.5)+
      scale_colour_manual(values=tp_col,limits=names(tp_col),drop=FALSE)+
      scale_linetype_manual(values=tp_type,limits=names(tp_col),drop=FALSE)+
      scale_shape_manual(values=tp_shape,limits=names(tp_col),drop=FALSE)+
      guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1))
  }
  p <- p+facet_tx()+scale_x_tx(tx)+scale_y_tx(pool,tx,metric)+
    labs(x=NULL,y=NULL,colour=NULL,shape=NULL,linetype=NULL)+style_tx
  if(!show_strip) p <- p+no_strip
  if(!show_x) p <- p+no_x
  p
}

# ---- Exact grid layout ------------------------------------------------------
# All component gtables are given the same width vector, including the left
# axis gutter. Both taxon plot regions are equally wide. Panel columns have
# exactly equal null weights, so all four plotting rectangles align in every
# row and the taxon headings sit over the actual panel centres.
# One separate legend per metric row, shared across both taxa.
get_legend <- function(p) {
  g <- ggplotGrob(p+theme(legend.position='bottom'))
  j <- which(g$layout$name=='guide-box-bottom')
  if(!length(j)) j <- which(g$layout$name=='guide-box') # older ggplot2
  stopifnot(length(j)==1,!inherits(g$grobs[[j]],'zeroGrob'))
  g$grobs[[j]]
}
panel_cols <- function(g) {
  z <- g$layout[grepl('^panel(-[0-9]+-[0-9]+)?$',g$layout$name),,drop=FALSE]
  stopifnot(nrow(z)==2,all(z$l==z$r))
  sort(z$l)
}
panel_rows <- function(g) {
  z <- g$layout[grepl('^panel(-[0-9]+-[0-9]+)?$',g$layout$name),,drop=FALSE]
  stopifnot(nrow(z)==2,all(z$t==z$b),length(unique(z$t))==1)
  unique(z$t)
}
combine_max <- function(u) do.call(grid::unit.pmax,u)
text_center <- function(label,size=8.2,bold=FALSE,rot=0) textGrob(label,
  rot=rot,gp=gpar(fontsize=size,fontfamily='sans',fontface=if(bold)'bold' else 'plain',col=ink))
row_tag <- function(label) textGrob(label,x=unit(0,'npc'),y=unit(1,'npc')-unit(1,'mm'),
  just=c('left','top'),gp=gpar(fontsize=11,fontface='bold',fontfamily='sans',col=ink))
row_labels <- list(A=expression(atop('Paired held-out '*Delta*D^2,'vs spread only')),
  Af=expression(atop('Paired held-out '*Delta*D^2,'vs spread only')),
  B=expression('Absolute contribution ('*Delta*D^2*')'),
  C='Conditional change in\nlogit dissimilarity',
  D='Conditional change in\ndissimilarity (pp)',K='Mean fitted\ndissimilarity')

build_figure <- function(pool,metrics,main=FALSE) {
  nr <- length(metrics)
  P <- lapply(seq_len(nr),function(i) setNames(lapply(taxa,function(tx)
    make_plot(pool,tx,metrics[i],show_strip=i==1,show_x=i==nr)),taxa))
  G <- lapply(P,function(pair) lapply(pair,function(p) ggplotGrob(p+theme(legend.position='none'))))
  flat <- unlist(G,recursive=FALSE)
  # Check layout signatures before aligning. This intentionally fails rather
  # than quietly accepting different panel positions across ggplot versions.
  pc <- panel_cols(flat[[1]])
  stopifnot(all(vapply(flat,function(g) length(g$widths)==length(flat[[1]]$widths)&&
      identical(panel_cols(g),pc),logical(1))))
  widths <- combine_max(lapply(flat,function(g) g$widths))
  widths[pc] <- unit(rep(1,length(pc)),'null')
  # Harmonize heights of the two grobs in each row; never match heights across
  # rows, since only the first has facet headers and only the last has x ticks.
  for(i in seq_len(nr)) {
    stopifnot(length(G[[i]][[1]]$heights)==length(G[[i]][[2]]$heights),
      identical(panel_rows(G[[i]][[1]]),panel_rows(G[[i]][[2]])))
    hh <- combine_max(lapply(G[[i]],function(g) g$heights))
    hh[panel_rows(G[[i]][[1]])] <- unit(1,'null')
    for(j in seq_along(taxa)) {
      G[[i]][[j]]$widths <- widths
      G[[i]][[j]]$heights <- hh
    }
  }
  # Outer grid: tag | metric title | macroinvertebrates | gap | zooplankton.
  # Use equal physical widths for the two taxon blocks. Shared gtables ensure
  # equal y-axis gutters, facet gaps and actual panel widths, not just equal
  # outer plot widths. No nesting of patchwork layouts is needed.
  ow <- unit.c(unit(5,'mm'),unit(10,'mm'),unit(1,'null'),unit(2.5,'mm'),unit(1,'null'),unit(2,'mm'))
  oh <- unit(6,'mm'); plot_row <- integer(nr); legend_row <- rep(NA_integer_,nr)
  weights <- if(main) c(1.10,1.13,1,1.07) else rep(1,nr)
  for(i in seq_len(nr)) {
    plot_row[i] <- length(oh)+1L
    oh <- unit.c(oh,unit(weights[i],'null'))
    # Main C and D share one legend below D. B needs two compact guide rows.
    if(!main || metrics[i]!='C') {
      legend_row[i] <- length(oh)+1L
      oh <- unit.c(oh,unit(if(metrics[i]=='B') 11 else 6,'mm'))
    }
    if(i<nr) oh <- unit.c(oh,unit(2,'mm'))
  }
  # Shared x label placed immediately after the last plot and before its legend.
  # Add the row after the final plot, shifting the final legend down one row.
  pos <- plot_row[nr]
  oh <- unit.c(oh[seq_len(pos)],unit(5,'mm'),oh[seq.int(pos+1,length(oh))])
  legend_row[!is.na(legend_row)&legend_row>pos] <- legend_row[!is.na(legend_row)&legend_row>pos]+1L
  g <- gtable(widths=ow,heights=oh,name='figure5_taxon_grid')
  add <- function(g,gro,t,l,b=t,r=l,name='element') gtable_add_grob(g,gro,t=t,l=l,b=b,r=r,clip='off',name=name)
  # Header uses the very same width vector as the taxon plot below it.
  for(j in seq_along(taxa)) {
    h <- gtable(widths=widths,heights=unit(1,'null'))
    h <- gtable_add_grob(h,text_center(taxon_name[[taxa[j]]],size=9.5,bold=TRUE),
      t=1,l=min(pc),r=max(pc),clip='off')
    g <- add(g,h,1,if(j==1) 3 else 5,name=paste0('header_',taxa[j]))
  }
  for(i in seq_len(nr)) {
    rr <- plot_row[i]; m <- metrics[i]
    if(main) g <- add(g,row_tag(m),rr,1,name=paste0('tag_',m))
    g <- add(g,text_center(row_labels[[m]],rot=90),rr,2,name=paste0('metric_',m))
    for(j in seq_along(taxa)) g <- add(g,G[[i]][[j]],rr,if(j==1) 3 else 5,
      name=paste0('plot_',m,'_',taxa[j]))
    if(!is.na(legend_row[i])) g <- add(g,get_legend(P[[i]][[1]]),legend_row[i],3,r=5,name=paste0('legend_',m))
  }
  g <- add(g,text_center('Zeta order, k'),pos+1,3,r=5,name='x_title')
  list(grob=g,plots=P,common_widths=widths,panel_columns=pc)
}

# ---- Validate actual axis scales and export --------------------------------
check_scales <- function(pool,tx,metric,p) {
  b <- ggplot_build(p); pp <- b$layout$panel_params
  stopifnot(length(pp)==2,
    isTRUE(all.equal(pp[[1]]$y.range,pp[[2]]$y.range)),
    isTRUE(all.equal(pp[[1]]$y$get_breaks(),pp[[2]]$y$get_breaks())),
    isTRUE(all.equal(pp[[1]]$x.range,pp[[2]]$x.range)))
  v <- metric_values(D[[pool]],tx,metric)
  sp <- scale_spec(pool,tx,metric)
  stopifnot(min(v)>=sp$lim[1]-1e-10,max(v)<=sp$lim[2]+1e-10)
  data.frame(pool=pool,taxon=tx,metric=metric,y_min=sp$lim[1],y_max=sp$lim[2],
    built_y_min=pp[[1]]$y.range[1],built_y_max=pp[[1]]$y.range[2],
    y_breaks=paste(sp$breaks,collapse=';'),n_facets=length(pp),
    all_values_in_range=TRUE,axis_scope=axis_scope)
}
checks <- bind_rows(lapply(names(D),function(pool) bind_rows(lapply(taxa,function(tx)
  bind_rows(lapply(c('A','Af','B','C','D','K'),function(m)
    check_scales(pool,tx,m,make_plot(pool,tx,m))))))))
if(axis_scope=='both') {
  q <- checks |> group_by(taxon,metric) |> summarise(n_range=n_distinct(y_min,y_max,y_breaks),.groups='drop')
  stopifnot(all(q$n_range==1))
}
write.csv(checks,file.path(out_dir,'Figure5_taxon_axis_checks.csv'),row.names=FALSE)
export <- function(g,name,h) {
  for(ext in c('pdf','svg','png')) {
    dev <- switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png')
    ggsave(file.path(out_dir,paste0(name,'.',ext)),plot=g,width=W,height=h,
      units='mm',dpi=600,device=dev,bg='white')
  }
}
dump <- function(d,name,suffix) write.csv(d,file.path(out_dir,paste0(name,'_',suffix,'.csv')),row.names=FALSE)
for(pool in names(D)) {
  d <- D[[pool]]
  nm <- if(pool=='shared') 'Figure5_final' else 'FigureS5_final_full_pool'
  fig <- build_figure(pool,c('A','B','C','D'),main=TRUE)
  export(fig$grob,nm,H)
  dump(d$F,nm,'A_folds'); dump(d$S,nm,'A_medians'); dump(d$I,nm,'B_importance')
  dump(mutate(d$C,response_change_pp=100*response_change),nm,'CD_TP_contrasts')
  dump(d$ranges,nm,'TP_common_ranges')
  sn <- paste0('FigureS5_final_panelA_folds_',pool)
  export(build_figure(pool,'Af')$grob,sn,86)
  dump(d$F,sn,'fold_values'); dump(d$S,sn,'medians')
  cn <- paste0('FigureS5_final_ceiling_',pool)
  export(build_figure(pool,'K')$grob,cn,80)
  dump(d$C,cn,'values')
  message('Saved ',nm,' and supplements; axis scope = ',axis_scope)
}
# Caption: rows have independent y scales for each taxon, but both
# representations of the same taxon share the same scale. With default
# axis_scope='both', shared and full-pool plots use identical scales as well.
# A: lines/points = medians of five paired spatial-fold contrasts, seed means
# within fold. Individual fold contrasts in supplement. B: mean signed
# Shapley/Owen allocations, log10 nonlinear-position fit. C/D: fitted endpoint
# contrasts over the common feasible TP range within each taxon; ranges differ
# across taxa and pools. D in percentage points. Conditional contrasts have no
# fold-level uncertainty. No refitting or new predictions performed.
