# Figure 6: synthetic eutrophication figure -- TP-specific fitted spline shapes.
# Run from docs/final_empirical_analysis:
# Rscript figure6_eutrophication_splines.R fig6_interaction/tables figures/figure6_splines
# Optional exact spread-scale profiles (no refit): FIG6_RAW_SPREAD=1 ZG_LIB6=... Rscript ...
# Inputs: outputs of 61_spline_shapes.R at commit 2726d73, package f7e115f.
# Main candidate: heterogeneity-share stage; companion: total-dissimilarity stage.
# These are fitted transformations, not causal effects or confidence intervals.
# Model-scale share is logit(H/T) in the coherent architecture, not logit(H).
# Pure spline f(z; p) has SITE TP (z) on x, not raw set spread. The optional
# spread panels evaluate mean pairwise signed spline differences for a fixed
# equally spaced log-TP configuration and screen TP spread-position support.
# Additive f is position independent; additive raw-spread effect need not be.
# Penalty fixed a priori for display; all three previously tested penalties saved.
# No fitting, selection by shape, p-values, or changes to the package.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(grid); library(gtable)
})
stopifnot(requireNamespace('svglite',quietly=TRUE))
args <- commandArgs(TRUE)
in_dir <- if(length(args)>=1) args[1] else 'fig6_interaction/tables'
out_dir <- if(length(args)>=2) args[2] else 'figures/figure6_splines'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
MAIN_PEN <- as.numeric(Sys.getenv('FIG6_PENALTY','1e-6'))
PENS <- c(1e-2,1e-4,1e-6)
stopifnot(length(MAIN_PEN)==1,MAIN_PEN %in% PENS)
TX <- c('BM','ZP'); POOLS <- c('shared','full'); STAGES <- c('share','total')
ORD <- list(BM=c(2L,10L),ZP=c(2L,5L)) # predeclared low/high contrast
TXL <- c(BM='EPA macroinvertebrates',ZP='EPA zooplankton')
POS <- c('Low','Intermediate','High'); TARGET <- c(Low=10,Intermediate=35,High=100)
ink <- '#252B2D'; grey <- '#848B8D'; gridcol <- '#E5E9E9'
COL <- c(Low='#79936C',Intermediate='#246878',High='#C45C19',Additive=grey)
LAB <- c(Low='Low (10)',Intermediate='Intermediate (35)',High='High (100)',Additive='Additive')
PANEL <- data.frame(group=rep(TX,each=2),order=unlist(ORD,use.names=FALSE))
pen_tag <- function(x) format(x,scientific=TRUE,trim=TRUE)
base <- theme_minimal(base_size=10,base_family='sans')+theme(
  text=element_text(colour=ink),axis.title=element_blank(),
  axis.text=element_text(size=8,colour=ink),axis.line=element_line(linewidth=.35,colour=ink),
  axis.ticks=element_line(linewidth=.3,colour=ink),axis.ticks.length=unit(1,'mm'),
  panel.grid.major=element_line(colour=gridcol,linewidth=.27),panel.grid.minor=element_blank(),
  legend.position='none',plot.margin=margin(1.5,2,1.5,1.5))
read_tab <- function(pool,stem) {
  p <- file.path(in_dir,paste0(stem,'_',pool,'.csv'))
  if(!file.exists(p)) stop('Missing ',p,'; use committed 61_spline_shapes.R outputs.')
  z <- read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(nrow(z)>0,all(z$pool==pool)); z
}
S <- bind_rows(lapply(POOLS,read_tab,stem='interaction_spline_values')) |>
  semi_join(PANEL,by=c('group','order')) |> filter(position %in% POS,
    model=='additive' | variant=='interaction_TP')
CF <- bind_rows(lapply(POOLS,read_tab,stem='interaction_coefficients_allsites')) |>
  semi_join(PANEL,by=c('group','order')) |> filter(model=='additive' | variant=='interaction_TP')
stopifnot(all(c('TP_log10','TP_ugL','f','f_centred','beta1','beta2','beta3','stage','position','u') %in% names(S)),
 !anyDuplicated(S[c('pool','group','order','model','stage','position','TP_log10')]),
 all(is.finite(S$f)),all(is.finite(S$TP_log10)),all(S$TP_ugL>0),
 max(abs(log10(S$TP_ugL)-S$TP_log10))<1e-8,
 all(abs(S$TP_position_ugL-unname(TARGET[S$position]))<1e-8))
# No silent use of a nonconverged/constant fit. Inspect these diagnostics if a stop occurs.
stopifnot(all(CF$converged %in% TRUE),!any(CF$stage_constant %in% TRUE))
check <- S |> group_by(pool,group,order,model,stage,position) |>
  arrange(TP_log10,.by_group=TRUE) |> summarise(n=n(),
  err=max(abs(f_centred-(f-first(f)))),.groups='drop')
stopifnot(all(check$n==60),all(check$err<1e-8))
# Check additive equivalence for all three positions before reducing it to ONE line.
ad <- S |> filter(model=='additive') |> group_by(pool,group,order,stage,TP_log10) |>
  summarise(n=n(),err=diff(range(f)),.groups='drop')
stopifnot(all(ad$n==3),all(ad$err<1e-8))
# Common concentration interval across the two orders, models and both pools,
# separately for each taxon. Keep the evaluated knots/coordinates (no smoothing).
X <- S |> group_by(group,pool,order,model,stage,position) |>
  summarise(lo=min(TP_log10),hi=max(TP_log10),.groups='drop') |>
  group_by(group) |> summarise(lo=max(lo),hi=min(hi),.groups='drop')
stopifnot(all(X$hi>X$lo))
clip_curve <- function(z){
  z <- arrange(z,TP_log10); lo<-z$lo[1]; hi<-z$hi[1]
  xx<-sort(unique(c(lo,z$TP_log10[z$TP_log10>lo & z$TP_log10<hi],hi)))
  ff<-approx(z$TP_log10,z$f,xout=xx,ties='ordered')$y
  data.frame(x=10^xx,y=ff-ff[1],TP_log10=xx,
    f_uncentred=ff,reference_TP=10^lo,TP_position_ugL=z$TP_position_ugL[1],
    u=z$u[1],beta1=z$beta1[1],beta2=z$beta2[1],beta3=z$beta3[1])
}
CURVES <- S |> left_join(X,by='group') |>
  group_by(pool,group,order,model,variant,penalty,stage,position) |>
  group_modify(~clip_curve(.x)) |> ungroup() |>
  mutate(series=ifelse(model=='additive','Additive',position)) |>
  filter(model!='additive' | position=='Low')
# Main spline selection fixed across taxa/orders; no cherry-picking of panels.
choose <- function(z,pen,stage){ z |> filter(.data$stage==.env$stage,
  model=='additive' | (variant=='interaction_TP' & abs(penalty-pen)<1e-12)) }
check_panel <- function(z,raw=FALSE){
  if(!nrow(z)) return(invisible(NULL))
  if(!raw) stopifnot(setequal(z$series,names(COL)))
  else stopifnot(n_distinct(z$position)==3,n_distinct(z$kind)==2)
  invisible(NULL)
}
axis_spec <- function(v){
  v<-v[is.finite(v)]; if(!length(v))return(list(lim=c(-.1,.1),breaks=c(-.1,0,.1)))
  r<-range(c(0,v)); if(diff(r)<1e-10)r<-c(-.01,.01)
  lim<-r+c(-1,1)*.07*diff(r); b<-sort(unique(c(0,breaks_extended(4)(lim))))
  list(lim=lim,breaks=b[b>=lim[1]&b<=lim[2]])
}
numfmt <- function(x){# sufficient decimals even for small fitted effects
  label_number(accuracy=if(max(abs(x),na.rm=TRUE)<.1).001 else if(max(abs(x),na.rm=TRUE)<1).01 else .1,trim=TRUE)(x)
}
# Plot widths and row heights explicitly equalised; avoid patchwork drift.
label <- function(s,size=10,bold=FALSE,rot=0) textGrob(s,rot=rot,
  gp=gpar(col=ink,fontsize=size,fontface=if(bold)'bold'else'plain'))
legend_grob <- function(raw=FALSE){
  if(!raw){
    xs<-c(.02,.275,.56,.81); labs<-unname(LAB); cols<-unname(COL)
    kids<-list(); for(i in seq_along(xs))kids<-c(kids,list(
      segmentsGrob(unit(xs[i],'npc'),unit(.5,'npc'),unit(xs[i]+.045,'npc'),unit(.5,'npc'),gp=gpar(col=cols[i],lwd=2,lty=if(i==4)2 else 1)),
      textGrob(labs[i],x=unit(xs[i]+.055,'npc'),hjust=0,gp=gpar(col=ink,fontsize=8))))
  }else{
    xs<-c(.02,.285,.59); kids<-list();for(i in 1:3)kids<-c(kids,list(
      segmentsGrob(unit(xs[i],'npc'),unit(.5,'npc'),unit(xs[i]+.045,'npc'),unit(.5,'npc'),gp=gpar(col=COL[i],lwd=2)),
      textGrob(LAB[i],x=unit(xs[i]+.055,'npc'),hjust=0,gp=gpar(col=ink,fontsize=8))))
  }; do.call(grobTree,kids)
}
make_panel <- function(z,ys,xlim,raw=FALSE){
  check_panel(z,raw)
  p<-ggplot(z,aes(x,y))+geom_hline(yintercept=0,colour='#A2AAAC',linewidth=.35)
  if(nrow(z)){
    if(!raw){p<-p+geom_line(data=filter(z,series!='Additive'),aes(colour=series,group=series),linewidth=1.05,lineend='round')+
      geom_line(data=filter(z,series=='Additive'),colour=grey,linetype='42',linewidth=.75)+
      scale_colour_manual(values=COL,limits=names(COL))}
    else {p<-p+geom_line(aes(colour=position,linetype=kind,group=interaction(position,kind)),linewidth=.85)+
      scale_colour_manual(values=COL[POS],limits=POS)+scale_linetype_manual(values=c(Additive='42',Interaction='solid'))}
  }else p<-p+annotate('text',x=mean(xlim),y=mean(ys$lim),label='No common supported\nspread interval',size=3,colour=grey)
  if(!raw){
    br<-c(1,3,10,30,100,300,1000,3000,10000); br<-br[br>=xlim[1]&br<=xlim[2]]
    # Limit tick density without removing data.
    if(length(br)>5)br<-br[seq(1,length(br),by=2)]
    p<-p+scale_x_log10(limits=xlim,breaks=br,labels=label_number(big.mark=','),expand=expansion(mult=.03))
  }else p<-p+scale_x_continuous(limits=xlim,breaks=breaks_extended(4),expand=expansion(mult=.035))
  p+scale_y_continuous(limits=ys$lim,breaks=ys$breaks,labels=numfmt,expand=expansion(mult=0))+base
}
assemble <- function(z,pool,stage,raw=FALSE){
  # same y scale for both orders AND pools within each taxon; independent taxa.
  # Both figures use actual model units, not normalised values.
  gs<-vector('list',4); checks<-list()
  for(j in 1:2)for(i in 1:2){tx<-TX[j];k<-ORD[[tx]][i];ii<-(i-1)*2+j
    zd<-z |> filter(group==tx); zp<-zd |> filter(.data$pool==.env$pool,order==k)
    ys<-axis_spec(zd$y)
    xl<-if(raw){if(nrow(zd))range(zd$x)else c(0,1)}else 10^unlist(X[X$group==tx,c('lo','hi')],use.names=FALSE)
    p<-make_panel(zp,ys,xl,raw); b<-ggplot_build(p)$layout$panel_params[[1]]
    stopifnot(all(abs(b$y.range-ys$lim)<1e-8),all(zp$y>=ys$lim[1]&zp$y<=ys$lim[2]))
    gs[[ii]]<-ggplotGrob(p)
    checks[[ii]]<-data.frame(pool=pool,group=tx,order=k,stage=stage,raw_spread=raw,
      y_min=ys$lim[1],y_max=ys$lim[2],x_min=xl[1],x_max=xl[2],n_values=nrow(zp))
  }
  stopifnot(length(unique(vapply(gs,function(x)length(x$widths),integer(1))))==1)
  ww<-do.call(unit.pmax,lapply(gs,function(g)g$widths));for(i in 1:4)gs[[i]]$widths<-ww
  for(ix in list(1:2,3:4)){hh<-unit.pmax(gs[[ix[1]]]$heights,gs[[ix[2]]]$heights);for(i in ix)gs[[i]]$heights<-hh}
  g<-gtable(widths=unit.c(unit(11,'mm'),unit(1,'null'),unit(5,'mm'),unit(1,'null')),
    heights=unit.c(unit(7,'mm'),unit(5,'mm'),unit(1,'null'),unit(5,'mm'),unit(5,'mm'),unit(1,'null'),unit(7,'mm'),unit(5,'mm'),unit(7,'mm'),unit(if(raw)5 else 0,'mm')))
  add<-function(o,t,l,b=t,r=l)g<<-gtable_add_grob(g,o,t=t,l=l,b=b,r=r,clip='off')
  for(j in 1:2){col<-2*j;add(label(TXL[j],11,TRUE),1,col)
    for(i in 1:2){ii<-(i-1)*2+j; rr<-if(i==1)2 else 5
      add(grobTree(label(paste0('k = ',ORD[[TX[j]]][i]),10,TRUE),
        textGrob(LETTERS[ii],x=unit(.015,'npc'),hjust=0,gp=gpar(col=ink,fontsize=11,fontface='bold'))),rr,col)
      add(gs[[ii]],rr+1,col)}}
  yl<-paste0(if(raw)'Centred TP spread contribution'else 'Centred TP I-spline contribution',
    if(stage=='share')'\n(logit heterogeneity share)'else '\n(logit total dissimilarity)')
  add(label(yl,9,FALSE,90),3,1,b=6)
  add(label(if(raw)'TP spread (mean pairwise absolute difference in log10 TP)'else expression('Site TP ('*mu*'g/L, logarithmic scale)'),10),7,2,r=4)
  add(label(expression('Geometric mean TP position ('*mu*'g/L)'),9),8,2,r=4)
  add(legend_grob(raw),9,2,r=4)
  if(raw)add(label('Solid: interaction     Dashed: additive',8),10,2,r=4)
  list(g=g,checks=bind_rows(checks))
}
export <- function(z,stage,pen,raw=FALSE){
  for(pool in POOLS){
    res<-assemble(z,pool,stage,raw)
    stem<-paste0('Figure6_TP_',if(raw)'spread_'else'spline_',stage,'_',pool,'_pen',pen_tag(pen))
    for(ext in c('pdf','svg','png')) ggsave(file.path(out_dir,paste0(stem,'.',ext)),res$g,
      width=183,height=139,units='mm',dpi=600,bg='white',device=switch(ext,pdf=grDevices::cairo_pdf,svg=svglite::svglite,png='png'))
    write.csv(filter(z,.data$pool==.env$pool),file.path(out_dir,paste0(stem,'_values.csv')),row.names=FALSE)
    write.csv(res$checks,file.path(out_dir,paste0(stem,'_axis_checks.csv')),row.names=FALSE)
  }
}
for(st in STAGES)for(pen in PENS)export(choose(CURVES,pen,st),st,pen)
write.csv(X,file.path(out_dir,'Figure6_common_site_TP_domain.csv'),row.names=FALSE)
# Export primary names without altering or replacing the prior Figure 6 results.
for(st in STAGES)for(pool in POOLS)for(ext in c('pdf','svg','png')){
  stem<-paste0('Figure6_TP_spline_',st,'_',pool,'_pen',pen_tag(MAIN_PEN))
  file.copy(file.path(out_dir,paste0(stem,'.',ext)),file.path(out_dir,paste0('Figure6_candidate_',st,'_',pool,'.',ext)),overwrite=TRUE)
}
# OPTIONAL exact raw-spread version: reconstruct only the original training DESIGN,
# NOT the fit. Read the already fitted spline coefficients. Same basis/sets as 61.
# The reference profile is equally spaced in log10 TP at each order (a hypothetical
# standardised profile), not a mean prediction over observed site combinations.
# TP-only support screen counts distinct observed sets and sites in 12x12 bins.
# This is bivariate TP support, not proof of full multivariate overlap.
if(Sys.getenv('FIG6_RAW_SPREAD','0')=='1'){
  .libPaths(c(Sys.getenv('ZG_LIB6'),.libPaths()));stopifnot(requireNamespace('zetaGDUM',quietly=TRUE))
  ns<-asNamespace('zetaGDUM'); getint<-function(n)get(n,envir=ns,inherits=FALSE)
  raw_out<-list(); support_out<-list()
  for(pool in POOLS){
    d<-readRDS(file.path(paste0('fig6_',pool),'out','prepared.rds'))
    for(tx in TX)for(k in ORD[[tx]]){
      e<-data.frame(TP=log10(d[[tx]]$env$TP),TN=log10(d[[tx]]$env$TN),Chla=log10(d[[tx]]$env$Chla))
      xy<-as.matrix(d[[tx]]$env[,c('lon','lat')]);idx<-zetaGDUM::sample_sets(nrow(e),k,6000,seed=7000+k)
      dd<-zetaGDUM::ugdm_component_design(e,idx,'additive',coordinates=xy,lonlat=TRUE)
      tb<-dd$basis;spl<-tb$splines[[match('TP',dd$vars)]];ctr<-tb$center[['TP']];scl<-tb$scale[['TP']];bnd<-attr(spl,'Boundary.knots')
      Bz<-function(z)unclass(predict(spl,pmin(pmax((z-ctr)/scl,bnd[1]),bnd[2])))[,1:3,drop=FALSE]
      # Verify against the committed fitted function values before constructing profiles.
      z0<-S |> filter(.data$pool==.env$pool,group==tx,order==k,model=='additive',stage=='total',position=='Low')
      beta0<-as.numeric(unlist(z0[1,c('beta1','beta2','beta3')],use.names=FALSE))
      stopifnot(max(abs(drop(Bz(z0$TP_log10)%*%beta0)-z0$f))<1e-7)
      tr<-getint('.component_support_cells')(list(position=dd$raw_position,spread=dd$raw_spread,idx=dd$idx),'TP',12L)
      v<-seq(-1,1,length.out=k);v<-v/mean(abs(outer(v,v,'-'))[upper.tri(matrix(0,k,k))])
      s<-seq(max(0,tr$spread_range[1]),tr$spread_range[2],length.out=301)
      ok<-sapply(TARGET,function(p){
        lp<-log10(p); cell<-getint('.component_cells')(rep(lp,length(s)),s,12L,tr$position_range,tr$spread_range)
        inside<-sapply(s,function(t)all(lp+t*v>=min(e$TP)&lp+t*v<=max(e$TP)))
        !is.na(cell$cell)&tr$sets[cell$cell]>=8&tr$sites[cell$cell]>=12&inside
      });good<-apply(ok,1,function(x)all(x %in% TRUE))
      rr<-rle(good);en<-cumsum(rr$lengths);be<-en-rr$lengths+1;run<-which(rr$values&rr$lengths>=3)
      if(length(run))run<-run[which.max(rr$lengths[run])]
      ss<-if(length(run))s[be[run]:en[run]]else numeric()
      # Require at least two distinct TP-spread cells, not just a narrow part of one bin.
      if(length(ss)){cc<-getint('.component_cells')(rep(log10(35),length(ss)),ss,12L,tr$position_range,tr$spread_range);if(n_distinct(cc$spread_bin)<2)ss<-numeric()}
      support_out[[length(support_out)+1]]<-data.frame(pool=pool,group=tx,order=k,n_grid=length(ss),
        spread_from=if(length(ss))min(ss)else NA_real_,spread_to=if(length(ss))max(ss)else NA_real_)
      if(!length(ss))next
      pp<-which(upper.tri(matrix(0,k,k)),arr.ind=TRUE)
      for(st in STAGES)for(pen in PENS)for(pos in POS)for(kind in c('Additive','Interaction')){
        q<-choose(S,pen,st) |> filter(.data$pool==.env$pool,group==tx,order==k,position==pos,
          if(kind=='Additive')model=='additive'else variant=='interaction_TP')
        stopifnot(nrow(q)==60);beta<-as.numeric(unlist(q[1,c('beta1','beta2','beta3')],use.names=FALSE))
        f<-function(s0){z<-log10(TARGET[pos])+s0*v;b<-Bz(z);dif<-abs(b[pp[,2],,drop=FALSE]-b[pp[,1],,drop=FALSE]);sum(colMeans(dif)*beta)}
        yy<-vapply(ss,f,numeric(1));raw_out[[length(raw_out)+1]]<-data.frame(pool=pool,group=tx,order=k,stage=st,
          display_pen=pen,position=pos,kind=kind,series=pos,x=ss,y=yy-yy[1],uncentred=yy,TP_position_ugL=TARGET[pos])
      }
    }
  }
  write.csv(bind_rows(support_out),file.path(out_dir,'Figure6_raw_spread_support.csv'),row.names=FALSE)
  if(length(raw_out)){R<-bind_rows(raw_out);for(st in STAGES)for(pen in PENS)export(filter(R,stage==st,display_pen==pen),st,pen,TRUE)}
  else warning('No supported raw-spread panels; spline figures are still available.')
}
message('Figure 6 candidates saved. Inspect both stages and all penalty sensitivities before manuscript selection.')
