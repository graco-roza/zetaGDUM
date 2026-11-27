# Figure 4 | BCI: local distinctiveness and multisite retention
# Run from docs/final_empirical_analysis:
# Rscript figure4_local_redesign.R tables/fig4 figures
# Uses only committed CSVs; no new fits, sampling or package installation.
# White halos behind the map labels keep Q1/Q2 legible on dark tiles.
# Graphical edits only (no data, calculation, selected-quadrat, statistic or axis-scale changes):
#  - maps enlarged to the width of the C-D columns (left/right margins and row heights below were calibrated by measuring the rendered panel edges), 2:1 aspect preserved;
#  - 200 m scale bar on the first map of each map row; 'Occurrence' relabelled 'Incidence'; symbols added to the Q1/Q2 legend of C; short leader lines for the Q1/Q2 labels of D; "Spearman's rho" with the rho symbol.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)})
args <- commandArgs(TRUE)
in_dir <- if(length(args)) args[1] else 'tables/fig4'
out_dir <- if(length(args)>1) args[2] else 'figures'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
# Layout constant (mm): panel height of the two map rows (= half the width of a C-D panel, 76.2 mm). Maps use coord_equal(), so the 2:1 aspect and the square 20 m quadrats are exact by construction.
FIG_W <- 183; FIG_H <- 201; hA <- 38.1; hB <- 38.1
d <- read.csv(file.path(in_dir,'site_metrics.csv'), check.names=FALSE) |> filter(set=='1985_tree')
stopifnot(nrow(d)==1250L, !anyDuplicated(d$quad), sum(d$habitat=='swamp')==30L)
qids <- c('x005_y024','x014_y012')
stopifnot(all(qids %in% d$quad))
d <- d |> mutate(LCBDa=LCBD_hel*1000, LCBDo=LCBD_hel_occ*1000,
                  group=ifelse(habitat=='swamp','Swamp','Other quadrats'),
                  pair=ifelse(quad==qids[1],'Q1',ifelse(quad==qids[2],'Q2',NA_character_)),
                  x=x_m+10,y=y_m+10)
rank_pct <- function(x) 100*(rank(x,ties.method='average')-1)/(length(x)-1)
ink <- '#252B2D'; petrol <- '#246878'; orange <- '#C45C19'; purple <- '#615A83'; grey <- '#AAB1B3'
qcol <- c(Q1=petrol,Q2=purple); qshape <- c(Q1=21,Q2=24)
reprs <- c('Relative abundance','Incidence')
base <- theme_minimal(base_size=9,base_family='sans')+theme(text=element_text(colour=ink), axis.text=element_text(size=8,colour=ink),axis.title=element_text(size=8.5),panel.grid.major=element_line(colour='#E5E9E9',linewidth=.25),panel.grid.minor=element_blank(),axis.line=element_line(colour=ink,linewidth=.3),axis.ticks=element_line(colour=ink,linewidth=.25),axis.ticks.length=grid::unit(1,'mm'),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=7.5),legend.key.width=grid::unit(8,'mm'),legend.box.spacing=grid::unit(1,'mm'),legend.margin=margin(0,0,0,0),plot.margin=margin(3,4,1,3),plot.tag=element_text(size=11,face='bold'))
# Within-representation percentile ranks remove the need for four colourbars.
# Higher ranks ALWAYS indicate higher values of the metric; raw metrics appear in D.
mapdat <- bind_rows(transmute(d,quad,habitat,pair,x,y,representation=reprs[1],metric='LCBD',rank=rank_pct(LCBDa)),
                    transmute(d,quad,habitat,pair,x,y,representation=reprs[2],metric='LCBD',rank=rank_pct(LCBDo)),
                    transmute(d,quad,habitat,pair,x,y,representation=reprs[1],metric='Retention',rank=rank_pct(T_rel)),
                    transmute(d,quad,habitat,pair,x,y,representation=reprs[2],metric='Retention',rank=rank_pct(T_occ)))
mapdat$representation <- factor(mapdat$representation,reprs)
# 200 m scale bar (first map of each row), on a faint white backing so it stays legible on the tiles.
sb <- data.frame(representation=factor(reprs[1],reprs),x0=45,x1=245,y0=38)
sb_bg <- data.frame(representation=factor(reprs[1],reprs),xmin=22,xmax=268,ymin=14,ymax=82)
map_plot <- function(metric_name,show_titles=FALSE){
  z <- filter(mapdat,.data$metric==.env$metric_name)
  p <- ggplot(z,aes(x,y))+geom_tile(aes(fill=rank),width=20,height=20)+
    geom_tile(data=filter(z,habitat=='swamp'),width=20,height=20,fill=NA,colour=orange,linewidth=.30)+
    geom_rect(data=sb_bg,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill='white',alpha=.7)+
    geom_segment(data=sb,aes(x=x0,xend=x1,y=y0,yend=y0),inherit.aes=FALSE,colour=ink,linewidth=.5,lineend='butt')+
    geom_text(data=sb,aes(x=(x0+x1)/2,y=y0+20,label='200 m'),inherit.aes=FALSE,colour=ink,size=2.4,vjust=0)+
    geom_point(data=filter(z,!is.na(pair)),aes(shape=pair),size=2.5,stroke=.65,fill='white',colour=ink,show.legend=FALSE)+
    geom_label(data=filter(z,!is.na(pair)),aes(x=x+ifelse(pair=='Q1',62,0),y=y+ifelse(pair=='Q1',0,-66),label=pair),colour=ink,fill=scales::alpha('white',.85),linewidth=0,label.padding=unit(1,'pt'),size=2.6,fontface='bold')+
    facet_wrap(vars(representation),nrow=1)+coord_equal(expand=FALSE,clip='off')+
    scale_shape_manual(values=qshape)+scale_fill_gradientn(colours=c('#F0F4EF','#BCD2C6','#6A9F9E',petrol,'#183D48'),limits=c(0,100),breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'),name='Within-representation percentile',guide=guide_colourbar(title.position='left',barwidth=grid::unit(45,'mm'),barheight=grid::unit(2.6,'mm')))+
    labs(x=NULL,y=NULL)+theme_void(base_size=9)+theme(strip.text=element_text(face='bold',size=9,margin=margin(b=5)),panel.spacing.x=grid::unit(5,'mm'),legend.position='bottom',legend.title=element_text(size=7.5),legend.text=element_text(size=7.5),plot.tag=element_text(face='bold',size=11),legend.box.spacing=grid::unit(1.2,'mm'),legend.margin=margin(0,0,0,0),plot.margin=margin(3,4,1,3))
  if(!show_titles) p<-p+theme(strip.text=element_blank())
  p
}
A <- map_plot('LCBD',TRUE)+labs(tag='A')+theme(legend.position='none')
B <- map_plot('Retention')+labs(tag='B')
# Recover exact focal profiles at the already evaluated orders from exported expected dissimilarity.
# E[mean set richness | i] - ED_inc = expected shared richness; divide by focal richness.
orders <- c(1,2,3,4,5,8,12,20,40,80,160,320,640)
stopifnot(all(paste0('ED_rel_',orders[-1]) %in% names(d)), all(paste0('ED_inc_',orders[-1]) %in% names(d)))
n <- nrow(d); total_r <- sum(d$richness)
profiles <- bind_rows(lapply(orders,function(k){
  zrel <- if(k==1) rep(1,n) else 1-d[[paste0('ED_rel_',k)]]
  zocc <- if(k==1) rep(1,n) else ((d$richness+(k-1)*(total_r-d$richness)/(n-1))/k-d[[paste0('ED_inc_',k)]])/d$richness
  bind_rows(data.frame(quad=d$quad,pair=d$pair,k=k,representation=reprs[1],z=zrel),data.frame(quad=d$quad,pair=d$pair,k=k,representation=reprs[2],z=zocc))
}))
stopifnot(all(is.finite(profiles$z)),all(profiles$z>=-1e-8 & profiles$z<=1+1e-8))
profiles <- profiles |> mutate(z=pmax(0,pmin(1,z)),representation=factor(representation,reprs))
env <- profiles |> group_by(representation,k) |> summarise(lo=quantile(z,.05),med=median(z),hi=quantile(z,.95),.groups='drop')
pairprof <- filter(profiles,!is.na(pair))
qlab <- c(Q1='Q1: slope, 73 species',Q2='Q2: mixed, 28 species')
C <- ggplot()+geom_ribbon(data=env,aes(k,ymin=lo,ymax=hi),fill=grey,alpha=.25)+geom_line(data=env,aes(k,med),colour='#82898B',linetype='dashed',linewidth=.45)+
 geom_line(data=pairprof,aes(k,z,colour=pair),linewidth=.8)+geom_point(data=pairprof,aes(k,z,colour=pair,shape=pair),size=1.4)+
 facet_wrap(vars(representation),nrow=1)+scale_x_log10(breaks=c(1,2,5,20,80,640))+scale_y_continuous(labels=label_percent(accuracy=1),breaks=seq(0,1,.25),limits=c(0,1.02))+
 scale_colour_manual(values=qcol,labels=qlab)+scale_shape_manual(values=c(Q1=16,Q2=17),labels=qlab)+
 guides(colour=guide_legend(override.aes=list(linetype='solid',size=1.8)))+
 labs(x='Zeta order, k',y='Focal composition\nretained',tag='C')+base+theme(strip.text=element_blank(),panel.spacing.x=grid::unit(5,'mm'),legend.key.width=grid::unit(9,'mm'))
# Retention and LCBD values on their original scales; only maps are percentile ranked.
scatter <- bind_rows(transmute(d,quad,group,pair,representation=reprs[1],LCBD=LCBDa,T=T_rel),transmute(d,quad,group,pair,representation=reprs[2],LCBD=LCBDo,T=T_occ)) |> mutate(representation=factor(representation,reprs))
rho <- scatter |> group_by(representation) |> summarise(r=cor(LCBD,T,method='spearman'),.groups='drop') |> mutate(label=sprintf('paste("Spearman\'s ", rho, " = ", %s)',formatC(round(r,2),format='f',digits=2)))
# Q1/Q2 labels: short leader lines from the true point coordinates to a nearby label position (offsets in axis units, kept inside the data range so the axis scales are unchanged).
off <- data.frame(representation=factor(c(reprs[1],reprs[1],reprs[2],reprs[2]),reprs),pair=c('Q1','Q2','Q1','Q2'),dx=c(-.30,.32,.10,.10),dy=c(-.60,.20,-.50,-.55))
lab <- filter(scatter,!is.na(pair)) |> left_join(off,by=c('representation','pair')) |> mutate(lx=LCBD+dx,ly=T+dy,sx=LCBD+.78*dx,sy=T+.78*dy)
D <- ggplot(scatter,aes(LCBD,T))+geom_point(data=filter(scatter,group=='Other quadrats'),aes(colour=group),size=.6,alpha=.6)+geom_point(data=filter(scatter,group=='Swamp'),aes(colour=group),size=1.25,alpha=.9)+
 geom_segment(data=lab,aes(x=LCBD,y=T,xend=sx,yend=sy),inherit.aes=FALSE,colour=ink,linewidth=.3)+
 geom_point(data=filter(scatter,!is.na(pair)),aes(shape=pair,fill=pair),colour=ink,size=2.7,stroke=.5,show.legend=FALSE)+
 geom_text(data=lab,aes(x=lx,y=ly,label=pair,colour=pair),size=2.6,fontface='bold',show.legend=FALSE,check_overlap=FALSE)+
 geom_text(data=rho,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.04,vjust=1.2,size=2.6,colour=ink,parse=TRUE)+facet_wrap(vars(representation),nrow=1,scales='free')+
 scale_shape_manual(values=qshape)+scale_fill_manual(values=qcol)+scale_colour_manual(values=c(qcol,`Other quadrats`="#A5ABAD",Swamp=orange),breaks=c("Other quadrats","Swamp"))+labs(x=expression('LCBD ('*10^-3*')'),y='Characteristic\nretention order, T',tag='D')+base+theme(strip.text=element_blank(),panel.spacing.x=grid::unit(5,'mm'))
# Use subtitles within maps rather than an extra panel or repeated headers.
A <- A+labs(subtitle='Local contribution to beta diversity (LCBD)')+theme(plot.subtitle=element_text(size=8.5,colour=ink,margin=margin(b=2),hjust=0))
B <- B+labs(subtitle='Characteristic retention order (T)')+theme(plot.subtitle=element_text(size=8.5,colour=ink,margin=margin(b=2),hjust=0))
fig <- A/B/C/D+plot_layout(heights=unit(c(hA,hB,1.25,1.28),c('mm','mm','null','null')))
name <- 'Figure4_main_redesign'
ggsave(file.path(out_dir,paste0(name,'.pdf')),fig,width=FIG_W,height=FIG_H,units='mm',device=grDevices::cairo_pdf,bg='white')
ggsave(file.path(out_dir,paste0(name,'.png')),fig,width=FIG_W,height=FIG_H,units='mm',dpi=600,bg='white')
ggsave(file.path(out_dir,paste0(name,'.svg')),fig,width=FIG_W,height=FIG_H,units='mm',device=svglite::svglite,bg='white')
write.csv(mapdat,file.path(out_dir,paste0(name,'_maps.csv')),row.names=FALSE)
write.csv(env,file.path(out_dir,paste0(name,'_profile_envelope.csv')),row.names=FALSE)
write.csv(pairprof,file.path(out_dir,paste0(name,'_profiles.csv')),row.names=FALSE)
write.csv(scatter,file.path(out_dir,paste0(name,'_scatter.csv')),row.names=FALSE)
