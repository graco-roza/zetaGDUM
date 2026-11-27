# Figure 6: environmental context and coherent multisite dissimilarity components
# Run from docs/final_empirical_analysis/:
# Rscript figure6_taxon_aligned.R fig6_shared/tables fig6_full/tables figures/figure6_redesign
# Plotting only. Reads outputs of 05_tables.R and 06_surface_tables.R.
# No raw data, zetaGDUM library, RDS files, new predictions or model refits required.
# Main: default interaction penalty (0.01), order 3 surfaces (pre-specified).
# Supplement: order 5, weaker penalty, validation fold ranges, penalty sensitivity.
# Dependencies: ggplot2, dplyr, tidyr, scales, grid, gtable, svglite.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(grid); library(gtable)
})
stopifnot(requireNamespace('svglite', quietly=TRUE))
args <- commandArgs(trailingOnly=TRUE)
paths <- c(shared=if(length(args)>=1) args[1] else 'fig6_shared/tables',
           full=if(length(args)>=2) args[2] else 'fig6_full/tables')
out_dir <- if(length(args)>=3) args[3] else 'figures/figure6_redesign'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

ink <- '#252B2D'; petrol <- '#246878'; orange <- '#C45C19'
purple <- '#615A83'; grey <- '#848B8D'; sage <- '#79936C'; grid_col <- '#E5E9E9'
W <- 183; H <- 198
TX <- c('BM','ZP'); TX_NAME <- c(BM='EPA macroinvertebrates',ZP='EPA zooplankton')
ORD <- list(BM=c(2,3,4,5,6,8,10,12),ZP=2:5)
COMP <- c('total','occupancy','heterogeneity')
COMP_LAB <- c(total='Total',occupancy='Occupancy-\nassociated',heterogeneity='Abundance\nheterogeneity')
PRED <- c('TP','TN','Chla','geographic_distance')
PRED_LAB <- c(TP='TP',TN='TN',Chla='Chl-a',geographic_distance='Geographic distance')
PRED_COL <- setNames(c(petrol,orange,purple,grey),PRED)
POS <- c('Oligotrophic','Mesotrophic','Eutrophic')
POS_LAB <- c(Oligotrophic='~10',Mesotrophic='~35',Eutrophic='~100')
POS_COL <- setNames(c(sage,petrol,orange),POS)
MOD <- c('additive','interaction','interaction_pen1e-4')
MOD_LAB <- c(additive='Additive',interaction='Explicit interaction',
             'interaction_pen1e-4'='Explicit interaction (0.0001)')
MOD_TYPE <- c(additive='solid',interaction='42','interaction_pen1e-4'='13')
GAIN <- c('position','interaction')
GAIN_LAB <- c(position='Add position',interaction='Add explicit interaction')
GAIN_COL <- c(position=petrol,interaction=orange)
GAIN_TYPE <- c(position='solid',interaction='42'); GAIN_SHAPE <- c(position=16,interaction=17)
PEN <- c('0.01','0.0001','0.000001')
PEN_COL <- setNames(c(petrol,orange,purple),PEN)
PEN_TYPE <- setNames(c('solid','42','13'),PEN)

base <- theme_minimal(base_size=8.5,base_family='sans')+theme(
  text=element_text(colour=ink),axis.title=element_blank(),
  axis.text=element_text(size=7,colour=ink),axis.line=element_line(colour=ink,linewidth=.3),
  axis.ticks=element_line(colour=ink,linewidth=.25),axis.ticks.length=unit(1,'mm'),
  panel.grid.major=element_line(colour=grid_col,linewidth=.25),panel.grid.minor=element_blank(),
  panel.spacing.x=unit(2,'mm'),strip.background=element_blank(),
  strip.text.x=element_text(size=7.3,face='bold',margin=margin(b=3)),
  plot.margin=margin(1,0,1,0),legend.position='bottom',legend.title=element_blank(),
  legend.text=element_text(size=7,colour=ink),legend.key.width=unit(6,'mm'),
  legend.key.height=unit(2.8,'mm'),legend.spacing.x=unit(1.2,'mm'),
  legend.margin=margin(0,0,0,0),legend.box.spacing=unit(0,'mm'))
no_strip <- theme(strip.text=element_blank(),strip.text.x=element_blank())
zero <- geom_hline(yintercept=0,colour='#A2AAAC',linewidth=.35)

read_tab <- function(pool,file,required) {
  p <- file.path(paths[[pool]],file)
  if(!file.exists(p)) stop('Missing ',p,'; run existing 05_tables.R / 06_surface_tables.R in the relevant fig6 directory (using existing RDS outputs). Do not refit just to plot.')
  message('Reading ',p)
  z <- read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(required %in% names(z))) stop(p,' lacks columns: ',paste(setdiff(required,names(z)),collapse=', '))
  z
}
# Preserve the actual order grid. Missing results stay NA rather than becoming zero.
check_key <- function(x,key,what) {
  if(anyDuplicated(x[key])) stop('Duplicate keys in ',what)
  invisible(x)
}
check_groups <- function(d,comp='component') {
  stopifnot(all(d$group %in% TX),all(d[[comp]] %in% COMP))
  for(tx in TX) stopifnot(all(d$order[d$group==tx] %in% ORD[[tx]]))
}
prepare <- function(pool) {
  I <- read_tab(pool,'importance_predictor.csv',c('group','chla','penalty','order','response','predictor','absolute_mean','abs_lo','abs_hi','n_valid_folds')) |>
    filter(chla,abs(as.numeric(penalty)-.01)<1e-12) |> rename(component=response)
  check_groups(I); stopifnot(all(I$predictor %in% PRED))
  check_key(I,c('group','order','component','predictor'),'importance')
  V <- read_tab(pool,'validation_gains.csv',c('group','chla','order','component','position_med','position_lo','position_hi','interaction_med','interaction_lo','interaction_hi','n_folds')) |> filter(chla)
  check_groups(V); check_key(V,c('group','order','component'),'validation gains')
  G <- bind_rows(
    transmute(V,group,order,component,effect='position',med=position_med,lo=position_lo,hi=position_hi,n_folds),
    transmute(V,group,order,component,effect='interaction',med=interaction_med,lo=interaction_lo,hi=interaction_hi,n_folds))
  P <- read_tab(pool,'interaction_penalty_sensitivity.csv',c('group','order','component','pen','gain_median','gain_min','gain_max','folds_positive')) |>
    mutate(penalty=case_when(abs(as.numeric(pen)-.01)<1e-12~'0.01',abs(as.numeric(pen)-.0001)<1e-12~'0.0001',abs(as.numeric(pen)-.000001)<1e-12~'0.000001',TRUE~NA_character_))
  check_groups(P); stopifnot(!anyNA(P$penalty)); check_key(P,c('group','order','component','penalty'),'penalty sensitivity')
  S <- read_tab(pool,'panelC_predictions.csv',c('group','order','model','trophic','TP_position','spread','total','occupancy','heterogeneity','training_sets','training_sites','evaluated_sets')) |>
    filter(model %in% MOD)
  stopifnot(all(S$group %in% TX),all(S$trophic %in% POS),all(is.finite(S$spread)),all(is.finite(S$TP_position)))
  check_key(S,c('group','order','model','trophic','spread'),'TP prediction cells')
  # Coherent decomposition: enforce closure and nonnegative component predictions.
  stopifnot(all(is.finite(as.matrix(S[COMP]))),all(as.matrix(S[COMP])>=-1e-8),all(as.matrix(S[COMP])<=1+1e-8),
    all(abs(S$total-S$occupancy-S$heterogeneity)<1e-7))
  if(any(I$n_valid_folds<5,na.rm=TRUE)||any(G$n_folds<5,na.rm=TRUE)) warning(pool,': fewer than five valid folds for some estimates; see plotted-value CSVs')
  list(I=I,G=G,P=P,S=S)
}
D <- lapply(names(paths),prepare); names(D) <- names(paths)

# Match the supported spread cells across EVERY displayed model and EVERY
# additive-selected trophic position. Never extrapolate or silently switch order.
# The upstream script already chose supported TP-position bins near 10,35,100.
# These are target labels, not exact concentrations or strict trophic classes.
select_cells <- function(pool,tx,k,models) {
  s <- D[[pool]]$S |> filter(group==tx,order==k,model %in% models)
  selected <- POS[POS %in% s$trophic[s$model=='additive']]
  available <- expand.grid(model=models,trophic=selected,stringsAsFactors=FALSE)
  if(nrow(available)) {
    cells <- lapply(seq_len(nrow(available)),function(i) {
      sort(unique(round(s$spread[s$model==available$model[i]&s$trophic==available$trophic[i]],6)))
    })
    common <- Reduce(intersect,cells)
  } else common <- numeric()
  good <- length(selected)>=2 && length(common)>=2
  if(good) {
    keep <- s |> filter(trophic %in% selected,round(spread,6) %in% common)
    # Each matched model and position must refer to the same actual TP position.
    pp <- keep |> group_by(trophic) |> summarise(r=diff(range(TP_position)),.groups='drop')
    stopifnot(all(pp$r<1e-7),nrow(keep)==length(common)*length(selected)*length(models))
    # Use original coordinates in figures/CSV. Round only the matching key.
    long <- keep |> pivot_longer(all_of(COMP),names_to='component',values_to='pred') |>
      arrange(component,model,trophic,spread)
    # Split segments when an intervening exported additive spread cell is
    # absent from the common grid. Lines connect supported summaries only;
    # the original full bin index is unavailable in this CSV.
    all_x <- sort(unique(round(D[[pool]]$S$spread[D[[pool]]$S$group==tx&D[[pool]]$S$order==k&D[[pool]]$S$model=='additive'],6)))
    long <- long |> group_by(component,model,trophic) |> arrange(spread,.by_group=TRUE) |>
      mutate(segment=cumsum(c(TRUE,diff(match(round(spread,6),all_x))>1L))) |> ungroup()
  } else {
    warning(pool,' ',tx,' k=',k,': no common supported range with >=2 positions and >=2 spread cells. Panel will say so.')
    long <- s[0,] |> pivot_longer(all_of(COMP),names_to='component',values_to='pred') |> mutate(segment=integer())
  }
  pos <- if(good) long |> distinct(trophic,TP_position) else data.frame(trophic=selected,TP_position=rep(NA_real_,length(selected)))
  info <- data.frame(pool=pool,group=tx,order=k,models=paste(models,collapse=';'),
    n_positions=length(selected),n_common_spread_bins=length(common),
    common_spread_min=if(length(common)) min(common) else NA_real_,
    common_spread_max=if(length(common)) max(common) else NA_real_,
    plotted=good,positions=paste(paste0(pos$trophic,'=',format(round(pos$TP_position,2),trim=TRUE)),collapse=';'))
  list(data=long,info=info,ok=good)
}
# Predeclared main k=3. k=5 and weaker penalty are sensitivities, not selected
# because their curves look more favourable. The default comparison is primary.
B_specs <- list(B3=list(k=3,models=c('additive','interaction')),
                B5=list(k=5,models=c('additive','interaction')),
                BL3=list(k=3,models=c('additive','interaction_pen1e-4')),
                BL5=list(k=5,models=c('additive','interaction_pen1e-4')))
for(pool in names(D)) {
  D[[pool]]$B <- lapply(B_specs,function(sp) setNames(lapply(TX,function(tx) select_cells(pool,tx,sp$k,sp$models)),TX))
}

# A/C/penalty: independent taxon-specific scales, fixed across the three
# components AND across the shared and full pools. B: natural dissimilarity
# scale (0-1); same for all components, to preserve the coherent decomposition.
# No within-panel normalization. Explicit scale limits/breaks checked below.
metric_values <- function(pool,tx,metric) {
  d <- D[[pool]]
  if(metric=='A') return(d$I$absolute_mean[d$I$group==tx])
  if(metric=='C') return(d$G$med[d$G$group==tx])
  if(metric=='Cf') return(unlist(d$G[d$G$group==tx,c('med','lo','hi')],use.names=FALSE))
  if(metric=='P') return(d$P$gain_median[d$P$group==tx])
  if(metric %in% names(B_specs)) return(d$B[[metric]][[tx]]$data$pred)
  stop('Unknown metric: ',metric)
}
yspec <- function(tx,metric) {
  if(metric %in% names(B_specs)) return(list(lim=c(0,1),breaks=c(0,.25,.5,.75,1)))
  v <- unlist(lapply(names(D),function(p) metric_values(p,tx,metric)),use.names=FALSE)
  v <- v[is.finite(v)]; stopifnot(length(v)>0)
  r <- range(c(0,v)); if(diff(r)==0) r <- r+c(-.01,.01)
  lim <- r+c(-1,1)*.06*diff(r)
  br <- sort(unique(c(0,breaks_extended(n=4)(lim)))); br <- br[br>=lim[1]&br<=lim[2]]
  list(lim=lim,breaks=br)
}
xspec <- function(tx,metric) {
  if(!metric %in% names(B_specs)) return(list(lim=range(ORD[[tx]]),breaks=if(tx=='BM')c(2,6,12)else 2:5))
  v <- unlist(lapply(names(D),function(p) D[[p]]$B[[metric]][[tx]]$data$spread))
  if(!length(v)) return(list(lim=c(0,1),breaks=c(0,.5,1)))
  lim <- range(v); if(diff(lim)==0) lim <- lim+c(-.01,.01)
  br <- breaks_extended(n=3)(lim); br <- br[br>=lim[1]&br<=lim[2]]
  if(length(br)<2) br <- lim
  list(lim=lim,breaks=br)
}
facet_tx <- function() facet_grid(cols=vars(component),scales='fixed',space='fixed',drop=FALSE,
  labeller=labeller(component=as_labeller(COMP_LAB)))
component_factor <- function(x) mutate(x,component=factor(component,levels=COMP))
# Fill missing orders with NA so no line bridges an absent evaluated result.
complete_orders <- function(x,tx,series) {
  lev <- switch(series,predictor=PRED,effect=GAIN,penalty=PEN,stop("Unknown series: ",series))
  grid <- expand.grid(order=ORD[[tx]],component=COMP,series_value=lev,stringsAsFactors=FALSE)
  names(grid)[3] <- series
  left_join(grid,x,by=c('order','component',series)) |> component_factor()
}
make_plot <- function(pool,tx,metric,show_strip=FALSE) {
  d <- D[[pool]]
  if(metric=='A') {
    x <- d$I |> filter(group==tx) |> complete_orders(tx,'predictor')
    p <- ggplot(x,aes(order,absolute_mean,colour=predictor,group=predictor))+zero+
      geom_line(linewidth=.65,na.rm=TRUE)+geom_point(size=1.1,na.rm=TRUE)+
      scale_colour_manual(values=PRED_COL,limits=PRED,labels=PRED_LAB,drop=FALSE)+
      guides(colour=guide_legend(nrow=1))
  } else if(metric %in% names(B_specs)) {
    b <- d$B[[metric]][[tx]]; x <- component_factor(b$data); mod <- B_specs[[metric]]$models
    p <- ggplot(x,aes(spread,pred,colour=trophic,linetype=model,group=interaction(trophic,model,segment)))+
      geom_line(linewidth=.65,na.rm=TRUE)+geom_point(size=.85,na.rm=TRUE)+
      scale_colour_manual(values=POS_COL,limits=POS,labels=POS_LAB,drop=FALSE,name=expression('Target TP ('*mu*g~L^-1*')'))+
      scale_linetype_manual(values=MOD_TYPE,limits=mod,labels=MOD_LAB[mod],drop=FALSE)+
      guides(colour=guide_legend(order=1,nrow=1,override.aes=list(linetype='solid')),
        linetype=guide_legend(order=2,nrow=1,override.aes=list(colour=ink)))+
      theme(legend.box='vertical',legend.title=element_text(size=7))
    if(!b$ok) {
      # Explicitly annotate every facet; never fabricate a line to fill a gap.
      stub <- data.frame(component=factor(COMP,levels=COMP),spread=.5,pred=.5,label='Insufficient\ncommon\nsupport')
      # Use the same x range as the other pool if it has supported observations.
      stub$spread <- mean(xspec(tx,metric)$lim)
      p <- p+geom_text(data=stub,aes(spread,pred,label=label),inherit.aes=FALSE,size=2.3,colour=grey)
    }
  } else if(metric %in% c('C','Cf')) {
    x <- d$G |> filter(group==tx) |> complete_orders(tx,'effect')
    p <- ggplot(x,aes(order,med,colour=effect,linetype=effect,shape=effect,group=effect))+zero
    if(metric=='Cf') p <- p+geom_linerange(aes(ymin=lo,ymax=hi),alpha=.42,linewidth=.35,
      show.legend=FALSE,na.rm=TRUE)
    p <- p+geom_line(linewidth=.72,na.rm=TRUE)+geom_point(size=1.35,na.rm=TRUE)+
      scale_colour_manual(values=GAIN_COL,limits=GAIN,labels=GAIN_LAB,drop=FALSE)+
      scale_linetype_manual(values=GAIN_TYPE,limits=GAIN,labels=GAIN_LAB,drop=FALSE)+
      scale_shape_manual(values=GAIN_SHAPE,limits=GAIN,labels=GAIN_LAB,drop=FALSE)+
      guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1),shape=guide_legend(nrow=1))
  } else if(metric=='P') {
    x <- d$P |> filter(group==tx) |> complete_orders(tx,'penalty')
    p <- ggplot(x,aes(order,gain_median,colour=penalty,linetype=penalty,group=penalty))+zero+
      geom_line(linewidth=.65,na.rm=TRUE)+geom_point(size=1.1,na.rm=TRUE)+
      scale_colour_manual(values=PEN_COL,limits=PEN,drop=FALSE,name='Interaction penalty')+
      scale_linetype_manual(values=PEN_TYPE,limits=PEN,drop=FALSE,name='Interaction penalty')+
      theme(legend.title=element_text(size=7))
  } else stop('Unknown metric: ',metric)
  ys <- yspec(tx,metric); xs <- xspec(tx,metric)
  p <- p+facet_tx()+scale_y_continuous(limits=ys$lim,breaks=ys$breaks,
    labels=label_number(accuracy=NULL,trim=TRUE),expand=expansion(mult=if(metric %in% names(B_specs))c(.025,.035)else 0))+
    scale_x_continuous(limits=xs$lim,breaks=xs$breaks,expand=expansion(mult=.06))+
    labs(colour=NULL,linetype=NULL,shape=NULL)+base
  # Explicit titles for B and penalty legends (base otherwise hides titles).
  if(metric %in% names(B_specs)) p <- p+labs(colour=expression('Target TP ('*mu*g~L^-1*')'))+theme(legend.title=element_text(size=7))
  if(metric=='P') p <- p+labs(colour='Interaction penalty',linetype='Interaction penalty')+theme(legend.title=element_text(size=7))
  if(!show_strip) p <- p+no_strip
  p
}

# ---- Fixed gtable layout, adapted from the successful Figure 5 renderer ----
panel_cols <- function(g) {
  z <- g$layout[grepl('^panel(-[0-9]+-[0-9]+)?$',g$layout$name),,drop=FALSE]
  stopifnot(nrow(z)==3,all(z$l==z$r)); sort(z$l)
}
panel_row <- function(g) {
  z <- g$layout[grepl('^panel(-[0-9]+-[0-9]+)?$',g$layout$name),,drop=FALSE]
  stopifnot(nrow(z)==3,all(z$t==z$b),length(unique(z$t))==1); unique(z$t)
}
max_units <- function(x) do.call(grid::unit.pmax,x)
text_center <- function(label,size=8,bold=FALSE,rot=0) textGrob(label,rot=rot,
  gp=gpar(fontsize=size,fontfamily='sans',fontface=if(bold)'bold' else 'plain',col=ink))
row_tag <- function(label) textGrob(label,x=unit(0,'npc'),y=unit(1,'npc')-unit(1,'mm'),
  just=c('left','top'),gp=gpar(fontsize=11,fontface='bold',fontfamily='sans',col=ink))
get_legend <- function(p) {
  g <- ggplotGrob(p+theme(legend.position='bottom'))
  j <- which(g$layout$name=='guide-box-bottom'); if(!length(j)) j <- which(g$layout$name=='guide-box')
  if(length(j)!=1||inherits(g$grobs[[j]],'zeroGrob')) return(NULL)
  g$grobs[[j]]
}
row_lab <- function(m) {
  if(m=='A') return(expression('Held-out contribution ('*Delta*D^2*')'))
  if(m %in% names(B_specs)) return('Predicted dissimilarity')
  if(m %in% c('C','Cf')) return(expression('Paired held-out '*Delta*D^2))
  expression('Interaction gain ('*Delta*D^2*')')
}
x_lab <- function(m) {
  if(m %in% names(B_specs)) return(bquote('TP spread: mean pairwise '*abs(Delta*log[10](TP))*'  (k = '*.(B_specs[[m]]$k)*')'))
  'Zeta order, k'
}
# B uses independent x axes, but fixed y range 0--1. A/C use separate taxon
# y scales. All rows and both blocks are aligned at the exact panel edges.
build_figure <- function(pool,metrics,main=FALSE) {
  nr <- length(metrics)
  P <- lapply(seq_len(nr),function(i) setNames(lapply(TX,function(tx)
    make_plot(pool,tx,metrics[i],show_strip=i==1)),TX))
  G <- lapply(P,function(pair) lapply(pair,function(p) ggplotGrob(p+theme(legend.position='none'))))
  flat <- unlist(G,recursive=FALSE); pc <- panel_cols(flat[[1]])
  stopifnot(all(vapply(flat,function(g) length(g$widths)==length(flat[[1]]$widths)&&identical(panel_cols(g),pc),logical(1))))
  ww <- max_units(lapply(flat,function(g) g$widths)); ww[pc] <- unit(rep(1,3),'null')
  for(i in seq_len(nr)) {
    stopifnot(length(G[[i]][[1]]$heights)==length(G[[i]][[2]]$heights),identical(panel_row(G[[i]][[1]]),panel_row(G[[i]][[2]])))
    hh <- max_units(lapply(G[[i]],function(g)g$heights)); hh[panel_row(G[[i]][[1]])] <- unit(1,'null')
    for(j in seq_along(TX)) {G[[i]][[j]]$widths <- ww; G[[i]][[j]]$heights <- hh}
  }
  # Tag and y title gutters are identical across rows. Both taxon blocks are
  # equally wide, and each contains three equal-width plotting rectangles.
  ow <- unit.c(unit(5,'mm'),unit(10,'mm'),unit(1,'null'),unit(2.5,'mm'),unit(1,'null'),unit(2,'mm'))
  oh <- unit(6,'mm'); pr <- xr <- lr <- integer(nr)
  for(i in seq_len(nr)) {
    m <- metrics[i]
    pr[i] <- length(oh)+1L; oh <- unit.c(oh,unit(if(m %in% names(B_specs))1.12 else 1,'null'))
    xr[i] <- length(oh)+1L; oh <- unit.c(oh,unit(5,'mm'))
    lr[i] <- length(oh)+1L; oh <- unit.c(oh,unit(if(m %in% names(B_specs))12 else 6,'mm'))
    if(i<nr) oh <- unit.c(oh,unit(2,'mm'))
  }
  g <- gtable(widths=ow,heights=oh,name='figure6_taxon_grid')
  add <- function(g,gro,t,l,b=t,r=l,name='element') gtable_add_grob(g,gro,t=t,l=l,b=b,r=r,clip='off',name=name)
  for(j in seq_along(TX)) {
    h <- gtable(widths=ww,heights=unit(1,'null'))
    h <- gtable_add_grob(h,text_center(TX_NAME[[TX[j]]],size=9.4,bold=TRUE),t=1,l=min(pc),r=max(pc),clip='off')
    g <- add(g,h,1,if(j==1)3 else 5,name=paste0('header_',TX[j]))
  }
  for(i in seq_len(nr)) {
    m <- metrics[i]
    if(main) g <- add(g,row_tag(LETTERS[i]),pr[i],1,name=paste0('tag_',m))
    g <- add(g,text_center(row_lab(m),rot=90),pr[i],2,name=paste0('metric_',m))
    for(j in seq_along(TX)) g <- add(g,G[[i]][[j]],pr[i],if(j==1)3 else 5,name=paste0('plot_',m,'_',TX[j]))
    g <- add(g,text_center(x_lab(m)),xr[i],3,r=5,name=paste0('x_',m))
    key <- NULL
    # Prefer a taxon with supported data for the curve legend.
    candidates <- if(m %in% names(B_specs)) which(vapply(D[[pool]]$B[[m]],function(z)z$ok,logical(1))) else seq_along(TX)
    for(j in candidates) {key <- get_legend(P[[i]][[j]]); if(!is.null(key)) break}
    if(is.null(key)) key <- text_center('No supported comparison available',size=7)
    g <- add(g,key,lr[i],3,r=5,name=paste0('legend_',m))
  }
  g
}

# ---- Programmatic validation ------------------------------------------------
check_scales <- function(pool,tx,m) {
  p <- make_plot(pool,tx,m); pp <- ggplot_build(p)$layout$panel_params
  stopifnot(length(pp)==3)
  for(j in 2:3) stopifnot(isTRUE(all.equal(pp[[1]]$y.range,pp[[j]]$y.range)),
    isTRUE(all.equal(pp[[1]]$y$get_breaks(),pp[[j]]$y$get_breaks())),
    isTRUE(all.equal(pp[[1]]$x.range,pp[[j]]$x.range)),
    isTRUE(all.equal(pp[[1]]$x$get_breaks(),pp[[j]]$x$get_breaks())))
  ys <- yspec(tx,m); xs <- xspec(tx,m); v <- metric_values(pool,tx,m); v <- v[is.finite(v)]
  inrange <- all(v>=ys$lim[1]-1e-8&v<=ys$lim[2]+1e-8); stopifnot(inrange)
  data.frame(pool=pool,taxon=tx,metric=m,n_facets=length(pp),
    x_min=xs$lim[1],x_max=xs$lim[2],x_breaks=paste(xs$breaks,collapse=';'),
    y_min=ys$lim[1],y_max=ys$lim[2],y_breaks=paste(ys$breaks,collapse=';'),
    built_y_min=pp[[1]]$y.range[1],built_y_max=pp[[1]]$y.range[2],
    n_finite_values=length(v),all_values_in_range=inrange)
}
checks <- list()
for(pool in names(D)) for(tx in TX) for(m in c('A','C','Cf','P',names(B_specs))) {
  checks[[length(checks)+1L]] <- check_scales(pool,tx,m)
}
checks <- bind_rows(checks)
check_pools <- checks |> group_by(taxon,metric) |> summarise(n_scales=n_distinct(x_min,x_max,x_breaks,y_min,y_max,y_breaks),.groups='drop')
stopifnot(all(check_pools$n_scales==1))
write.csv(checks,file.path(out_dir,'Figure6_axis_checks.csv'),row.names=FALSE)

# ---- Export -----------------------------------------------------------------
export <- function(g,name,height) {
  for(ext in c('pdf','svg','png')) {
    dev <- switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png')
    ggsave(file.path(out_dir,paste0(name,'.',ext)),plot=g,width=W,height=height,
      units='mm',dpi=600,device=dev,bg='white')
  }
}
dump <- function(x,name,suffix) write.csv(x,file.path(out_dir,paste0(name,'_',suffix,'.csv')),row.names=FALSE)
for(pool in names(D)) {
  d <- D[[pool]]; nm <- if(pool=='shared')'Figure6_main_redesign' else 'FigureS6_full_pool_redesign'
  export(build_figure(pool,c('A','B3','C'),main=TRUE),nm,H)
  dump(d$I,nm,'A_importance'); dump(d$G,nm,'C_validation_gains')
  # Main B and supplementary surfaces retain exact actual TP positions, model
  # labels, support counts and predicted components for all displayed cells.
  for(m in names(B_specs)) {
    s <- d$B[[m]]; v <- bind_rows(lapply(s,`[[`,'data')); info <- bind_rows(lapply(s,`[[`,'info'))
    stem <- if(m=='B3') nm else paste0('FigureS6_TP_',m,'_',pool)
    if(m!='B3') export(build_figure(pool,m),stem,85)
    dump(v,stem,'B_supported_predictions'); dump(info,stem,'B_support_diagnostics')
    # Numerical contrast over the identical supported endpoints within each
    # taxon, order, model and position; no new model predictions.
    slopes <- v |> group_by(group,order,component,model,trophic) |> arrange(spread,.by_group=TRUE) |>
      summarise(TP_position=first(TP_position),n_bins=n(),spread_from=first(spread),spread_to=last(spread),
        response_from=first(pred),response_to=last(pred),delta=last(pred)-first(pred),
        slope=ifelse(last(spread)>first(spread),(last(pred)-first(pred))/(last(spread)-first(spread)),NA_real_),.groups='drop')
    dump(slopes,stem,'B_endpoint_contrasts')
  }
  sn <- paste0('FigureS6_validation_fold_ranges_',pool)
  export(build_figure(pool,'Cf'),sn,83); dump(d$G,sn,'values')
  pn <- paste0('FigureS6_interaction_penalty_',pool)
  export(build_figure(pool,'P'),pn,83); dump(d$P,pn,'values')
  message('Saved ',nm,' and supplements')
}
# Caption notes (not printed on figure):
# A: mean signed held-out predictor allocations from the coherent interaction
# model, default penalty 0.01; aggregated over seeds and spatial folds, using
# log10 TP, TN and Chl-a plus geographic distance. Allocations are not tests.
# B: supported fitted surface cells, order 3, near target geometric-mean TP
# positions of 10/35/100 µg/L (actual values in B_support_diagnostics and CSV).
# Curves compare additive and explicit interaction (0.01) on common spread
# bins across all displayed positions and both models. Surface curves are
# model-based summaries of supported cells, not exact fixed-covariate partial
# effects or causal effects. The additive model can have context-dependent
# response-scale slopes through the link and coherent component construction.
# C: median paired spatial-fold D2 differences after averaging combination
# repeats within fold. Position = additive minus spread-only; interaction =
# explicit interaction minus additive. Fold ranges are supplementary, not CIs.
# A/C: y scales shared by the three components within a taxon and across pools;
# independent between taxa. B: common natural 0--1 dissimilarity scale.
# No model refit or new prediction performed; all original analysis files kept.
