source("docs/experimental/spatiotemporal/phase4b_annual_metacommunities/scripts/common.R")
library(data.table)
library(ggplot2)
library(patchwork)
theme_set(theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold",size=15,margin=margin(b=5)),
  plot.subtitle=element_text(colour="#475569",size=10,margin=margin(b=10)),
  plot.caption=element_text(colour="#475569",hjust=0), legend.position="bottom",
  panel.spacing=grid::unit(0.85,"lines"), plot.margin=margin(12,17,12,12)))
colours <- c(expansion="#007A78",contraction="#B14724",redistribution="#6654A4",numerical="#536780")
source_note <- "BioTIME 604 · seven traps, 2004–2013 · below the empirical validation minimum"
k_order <- function(k) factor(k,levels=1:6,labels=paste("Spatial order",1:6))
save_plot <- function(plot,name,w=9,h=5.5) {
  ggsave(file.path(figures,paste0(name,".pdf")),plot,width=w,height=h,device=cairo_pdf,bg="white")
  ggsave(file.path(figures,paste0(name,".png")),plot,width=w,height=h,device=ragg::agg_png,dpi=300,bg="white")
}
diverge <- scale_fill_gradient2(low="#117878",mid="#f7f8fa",high="#b9482d",midpoint=0)
colour_process <- scale_colour_manual(values=colours,name=NULL)
fill_process <- scale_fill_manual(values=colours,name=NULL)

coverage <- read_result("sampling_effort")
membership <- read_result("balanced_panel_membership")
coverage <- merge(coverage,membership[,.(site,retained)],by="site")
coverage[, label := paste0(sprintf("%04d", as.integer(site)),ifelse(retained," *",""))]
coverage[, symbol := fcase(retained & year%between%c(2004,2013) & status=="observed community","Retained",
                           status=="observed community","Archived outside panel",
                           default="No BioTIME observation")]
p <- ggplot(coverage,aes(year,label,fill=symbol))+geom_tile(color="white",linewidth=.17)+
  scale_fill_manual(values=c(Retained="#007A78","Archived outside panel"="#9DB5BF","No BioTIME observation"="#ECF0F3"),name=NULL)+
  scale_x_continuous(breaks=seq(1993,2016,3))+
  labs(title="Annual moth monitoring coverage",subtitle="Panel membership is determined by coverage alone",x=NULL,y="Trap code",
       caption="* Retained in the archive illustration. Trap nights and seasonal dates are unavailable.")
save_plot(p,"01-monitoring")

spatial <- read_result("spatial")[k>=2 & response %in% c("zeta_presence","incidence","turnover","nestedness")]
spatial[, metric:=factor(response,levels=c("zeta_presence","incidence","turnover","nestedness"),
  labels=c("Shared richness","Incidence dissimilarity","Full-intersection turnover","Full-intersection nestedness"))]
p <- ggplot(spatial,aes(year,value,colour=factor(k)))+geom_line(linewidth=.75)+geom_point(size=1.2)+
  facet_wrap(~metric,scales="free_y",ncol=2)+scale_colour_viridis_d(option="D",begin=.14,end=.86)+
  labs(title="Spatial diversity through time",subtitle="Every spatial combination in the fixed seven-trap panel",
       x=NULL,y=NULL,colour="Spatial order",caption=source_note)
save_plot(p,"02-spatial-trajectories",9,6.3)

annual <- read_result("annual")[k>=2 & response=="incidence"]
annual[, order:=k_order(k)]
net <- annual[, .(net=sum(net)),by=.(final,order)]
p <- ggplot(annual,aes(final,net,fill=process))+geom_hline(yintercept=0,colour="#8792A0")+
  geom_col(width=.72)+geom_point(data=net,aes(final,net),inherit.aes=FALSE,size=1.45)+
  facet_wrap(~order,ncol=3)+fill_process+
  labs(title="Annual occupancy-process contributions",subtitle="Points show total change. Negative contributions homogenise communities.",
       x="Final year of transition",y="Incidence dissimilarity change",caption=source_note)
q <- ggplot(annual[k==2],aes(final,net,colour=process))+
  geom_hline(yintercept=0,colour="#8792A0")+geom_linerange(aes(ymin=net-spatial_block_se,
    ymax=net+spatial_block_se),position=position_dodge(width=.55))+
  geom_point(position=position_dodge(width=.55),size=1.65)+colour_process+guides(colour="none")+
  labs(x="Final year of transition",y="Pairwise contribution",subtitle="Bars show ±1 spatial-block deletion sensitivity")
save_plot((p/q)+plot_layout(heights=c(2.05,1),guides="collect") & theme(legend.position="bottom"),
          "03-annual-contributions",10,8.4)

cum <- read_result("cumulative")[response=="incidence"&k>=2]
cum[,order:=k_order(k)]
p <- ggplot(cum,aes(year,cumulative_net,colour=process))+geom_hline(yintercept=0,colour="#8792A0")+
  geom_line(linewidth=.85)+geom_point(size=1.3)+facet_wrap(~order,ncol=3)+colour_process+
  labs(title="Cumulative annual process contributions",subtitle="Each curve sums consecutive-year signed changes on fixed support",
       x=NULL,y="Cumulative incidence dissimilarity change",caption=source_note)
save_plot(p,"04-cumulative-processes",10,6)

gross <- melt(cum[,.(year,k,order,process,cumulative_homogenising,cumulative_differentiating)],
              id.vars=c("year","k","order","process"),variable.name="direction",value.name="magnitude")
gross[, direction := factor(direction,levels=c("cumulative_homogenising","cumulative_differentiating"),
                            labels=c("Homogenising","Differentiating"))]
p <- ggplot(gross,aes(year,magnitude,colour=process,linetype=direction))+
  geom_line(linewidth=.85)+facet_wrap(~order,ncol=3)+colour_process+
  labs(title="Gross changes accumulate in both directions",subtitle="Positive and negative effects are summed separately before cancellation",
       x=NULL,y="Cumulative magnitude",linetype=NULL,caption=source_note)
save_plot(p,"05-cumulative-gross",10,6)

cmp <- read_result("comparison")[response=="incidence"&k>=2]
cmp[, order:=k_order(k)]
long <- melt(cmp,id.vars=c("order","process"),measure.vars=c("direct","accumulated"),
             variable.name="estimand",value.name="value")
long[, estimand := factor(estimand,labels=c("Direct endpoint","Accumulated annual"))]
p <- ggplot(long,aes(value,process,colour=process,shape=estimand))+geom_vline(xintercept=0,colour="#8792A0")+
  geom_segment(data=cmp,aes(x=direct,xend=accumulated,y=process,yend=process),inherit.aes=FALSE,colour="#8792A0")+
  geom_point(size=2.55)+facet_wrap(~order,ncol=3)+colour_process+guides(colour="none")+
  labs(title="Endpoint attribution misses the annual path",subtitle="Both paths sum to the same observed endpoint diversity difference",
       x="Incidence dissimilarity change",y=NULL,shape=NULL,caption=source_note)
save_plot(p,"06-endpoint-versus-annual",10,6)

local_reorg <- read_result("annual")[k==1 & response=="site_temporal_presence" & process=="expansion",.(final,local=endpoint_after)]
spatial_change <- annual[k==2,.(spatial=sum(net)),by=final]
reorg <- merge(local_reorg,spatial_change,by="final")
p <- ggplot(reorg,aes(spatial,local,colour=final))+
  geom_vline(xintercept=0,colour="#8792A0")+geom_point(size=3)+
  geom_text(aes(label=final),vjust=-.8,size=3.1,show.legend=FALSE)+
  scale_colour_viridis_c(begin=.18,end=.78,guide="none")+
  scale_y_continuous(limits=c(0,.302))+
  labs(title="Annual reorganisation despite small net spatial changes",
       subtitle="Temporal distance and spatial change use distinct incidence estimands",
       x="Annual change in pairwise spatial Sørensen",y="Mean local temporal Sørensen",
       caption=source_note)
save_plot(p,"07-temporal-reorganisation",9,5.5)

g <- read_result("window")[response=="incidence" & k>=2]
g[, text := ifelse(supported==0,"NA",sprintf("%+.4f",net))]
g[,process:=factor(process,levels=processes)]
p <- ggplot(g,aes(m,k,fill=net))+geom_tile(colour="white",linewidth=.6)+
  geom_text(aes(label=text),size=3.0)+facet_wrap(~process,nrow=1)+diverge+
  scale_x_continuous(breaks=1:5)+scale_y_continuous(breaks=2:6)+
  labs(title="Occupancy contributions across spatial and temporal orders",
       subtitle="Earliest aligned windows start in 2004 and 2009. Each survey pair is five years apart.",
       x="Temporal order",y="Spatial order",fill="Contribution",caption=source_note)
save_plot(p,"08-orders",10,4.6)

ab <- read_result("count_sensitivity")[k>=2 & response %in% c("incidence","abundance") & process=="expansion"]
ab[,metric:=factor(response,levels=c("incidence","abundance"),labels=c("Incidence","Relative trap catch"))]
p <- ggplot(ab,aes(final,after,colour=metric))+geom_line(linewidth=.85)+geom_point(size=1.3)+
  facet_wrap(~k_order(k),ncol=3)+scale_colour_manual(values=c(Incidence="#007A78","Relative trap catch"="#B14724"))+
  labs(title="Incidence and relative-catch dissimilarity",subtitle="Trap counts are unstandardised for trap nights and season length",
       x=NULL,y="Spatial dissimilarity",colour=NULL,caption=source_note)
save_plot(p,"09-count-sensitivity",10,6)

sim <- read_result("simulation_ecological")[n==60]
sim[, scenario := factor(scenario,levels=c("null","expansion","contraction","opposing","reversal","shocks",
                                          "stable_redistribution","detection"),
  labels=c("Null","Expansion","Contraction","Opposing","Reversal","AR(1) shocks","Stable spatial","Missed detections"))]
long <- melt(sim,id.vars=c("scenario","replicate"), measure.vars=c("mean_annual_local","mean_annual_latent_local",
  "endpoint_local","exact_pair_change"),variable.name="measure")
long[,measure:=factor(measure,labels=c("Annual observed local","Annual latent local",
  "Direct endpoint local","Exact pairwise spatial change"))]
s <- long[,.(average=mean(value),lo=quantile(value,.1),hi=quantile(value,.9)),by=.(scenario,measure)]
p <- ggplot(s,aes(scenario,average))+geom_hline(yintercept=0,colour="#9BABB7")+
  geom_linerange(aes(ymin=lo,ymax=hi),colour="#007A78")+geom_point(colour="#007A78",size=1.7)+
  facet_wrap(~measure,ncol=2,scales="free_y")+coord_flip()+
  labs(title="Continuous trajectories reveal distinct validation signals",
       subtitle="Points are replicate means. Bars show the 10th–90th percentile across panels.",
       x=NULL,y="Dissimilarity or change",caption="60 sites, 200 taxa, 30 years per primary simulation panel.")
save_plot(p,"10-continuous-validation",10,7)
