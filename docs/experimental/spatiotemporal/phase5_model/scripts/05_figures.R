source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
library(ggplot2)
pal <- c(expansion="#11766C",contraction="#BD4F28",redistribution="#7059A3",numerical="#A37222")
theme_set(theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
    legend.position="bottom",plot.title.position="plot",plot.caption.position="plot",
    strip.text=element_text(face="bold"),axis.title=element_text(face="bold")))
save <- function(p,name,width=8.5,height=5) {
  ggsave(file.path(phase,"figures",paste0(name,".pdf")),p,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(phase,"figures",paste0(name,".png")),p,width=width,height=height,dpi=260,bg="white")
}
C <- read_result("empirical_coverage")
C <- subset(C,dataset!="moths_balanced")
C$task <- factor(C$task,levels=c("spatial","local","local_observed","rectangle","change","change_joint","shared"),
                 labels=c("Spatial","Local, 1 year","Local, observed","Rectangle","Pair change","Window change","Shared turnover"))
save(ggplot(C,aes(task,dataset,fill=eligible_observation_fraction))+
  geom_tile(colour="white",linewidth=.8)+geom_text(aes(label=sprintf("%.1f%%",100*eligible_observation_fraction)),size=3.2)+
  scale_fill_gradient(low="#e5f1f2",high="#4ab0b3",limits=c(0,1),name="Observed records used")+
  labs(x=NULL,y=NULL,title="Observation coverage by response",caption="Eligibility is computed before combination sampling and response denominators.")+
  theme(axis.text.x=element_text(angle=20,hjust=1)),"01_coverage",9,4)

T <- read_result("empirical_trajectory")
T <- subset(T,metric=="incidence")
save(ggplot(T,aes(time,value,colour=panel,linetype=panel))+
  geom_line(linewidth=.65)+facet_wrap(~dataset,ncol=2)+
  scale_colour_manual(values=c("Phase IV balanced"="#975F37","response specific"="#16677B"))+
  labs(x="Survey year",y="Mean spatial Sørensen",colour=NULL,linetype=NULL,
       title="Spatial dissimilarity across observed communities",
       caption="Pairwise means use all observed pairs. The moth fixed-panel series is shown in its own panel."),
  "02_spatial_trajectory",9,6)

L <- subset(read_result("empirical_transitions"),task%in%c("local","local_observed")&metric=="site_temporal_presence")
L$contrast <- ifelse(L$task=="local","One-year contrasts","Consecutive observed surveys")
save(ggplot(L,aes(final,endpoint,colour=contrast))+
  geom_line(data=subset(L,task=="local"),linewidth=.65)+geom_point(size=1.3,alpha=.7)+
  facet_wrap(~dataset,ncol=2)+scale_colour_manual(values=c("#C46335","#247888"))+
  labs(x="End year",y="Mean local temporal Sørensen",colour=NULL,title="Turnover between local surveys",
       caption="Consecutive observed surveys may span gaps. Each point uses its available site-times."),
  "03_local_turnover",9,6)

P <- subset(read_result("empirical_predictions"),task=="change")
P <- aggregate(cbind(obs,pred)~dataset+initial+process,P,mean)
PP <- rbind(transform(P,value=obs,series="Observed"),transform(P,value=pred,series="Held-out prediction"))
save(ggplot(PP,aes(initial,value,colour=process,linetype=series))+
  geom_hline(yintercept=0,colour="grey60")+geom_line(linewidth=.6)+
  facet_wrap(~dataset,ncol=1)+scale_colour_manual(values=pal)+
  labs(x="Start year",y="Mean contribution (Sørensen units)",colour="Process",linetype=NULL,
       title="Observed and predicted process contributions",
       caption="Predictions come from temporal blocks. Only common held-out rows enter both series."),
  "04_process_predictions",9,7)

curves <- read_result("empirical_curves")
save(ggplot(curves,aes(change_sd_per_site,contribution,colour=process,group=interaction(process,fold)))+
  geom_hline(yintercept=0,colour="grey70")+geom_line(alpha=.63,linewidth=.45)+
  facet_grid(dataset~predictor)+scale_colour_manual(values=pal)+
  labs(x="Directional change per site (training SD)",y="Predicted change (Sørensen units)",colour="Process",
       title="Partial environmental response curves",
       caption="Each line uses one training fold. Other design columns are held at their training medians."),
  "05_response_curves",10,7)

I <- subset(read_result("empirical_importance"),task=="change" & method=="temporal")
BP <- bind_rows(lapply(split(I,I$dataset),function(i) cbind(dataset=i$dataset[1],st5_block_stability(i))))
write_result(BP,"empirical_block_stability")
IP <- subset(BP,feature_type=="predictor" & response%in%names(pal))
save(ggplot(IP,aes(feature,mean_fold_importance,colour=response))+
  geom_hline(yintercept=0,colour="grey70")+geom_pointrange(aes(ymin=lower,ymax=upper),
    position=position_dodge(width=.65),linewidth=.35)+facet_wrap(~dataset,ncol=1)+scale_colour_manual(values=pal)+
  labs(x=NULL,y="Held-out importance",colour="Process",title="Predictor importance for spatial change",
       caption="Points average foldwise predictive Shapley values. Ranges resample whole held-out time blocks.")+
  theme(axis.text.x=element_text(angle=15,hjust=1)),"06_predictor_importance",9,7)

Ig <- subset(read_result("empirical_importance"),fold=="pooled" & method=="temporal" &
  feature_type=="group" & response=="total" & task%in%c("local","spatial","rectangle","change"))
Ig$feature <- gsub("_"," ",Ig$feature)
save(ggplot(Ig,aes(feature,importance,fill=task))+
  geom_hline(yintercept=0,colour="grey70")+geom_col(position="dodge")+
  facet_wrap(~dataset,ncol=1)+scale_fill_brewer(palette="Dark2")+
  labs(x=NULL,y="Held-out fractional loss improvement",fill="Response",title="Spatial and temporal feature importance",
       caption="Feature-group Shapley values use fixed test observations within each response.")+
  theme(axis.text.x=element_text(angle=22,hjust=1)),"07_feature_importance",9,7)

S <- subset(read_result("simulation_importance"),fold=="pooled" & response%in%c("joint","expansion","contraction","redistribution"))
S <- S[is.finite(S$importance),]
ranks <- bind_rows(lapply(split(S,interaction(S$scenario,S$replicate,S$response,drop=TRUE)),function(x)
  transform(x,rank=rank(-importance,ties.method="average"),top=rank(-importance,ties.method="min")==1)))
top <- aggregate(top~scenario+feature+response,ranks,mean)
write_result(top,"simulation_rank_stability")
save(ggplot(top,aes(feature,top,fill=scenario))+
  geom_col(position="dodge")+facet_wrap(~response,ncol=2)+
  scale_fill_manual(values=c(multiple="#16677B",null="#A9926B",reference="#79A9AD"))+
  scale_y_continuous(labels=scales::percent,limits=c(0,1))+
  labs(x=NULL,y="Frequency ranked first",fill="Simulation",title="Predictor rank stability across independent panels",
       caption="Reference panels use independent seeds and a larger uniform combination budget."),
  "08_rank_recovery",8.5,6)

M <- subset(read_result("empirical_models"),method=="temporal" &
  task%in%c("spatial","local","change","rectangle","crossed") & mode!="ridge")
save(ggplot(M,aes(model,rmse,colour=model))+
  geom_point(position=position_jitter(width=.12,height=0,seed=3),size=1.9,alpha=.78)+
  facet_grid(dataset~task,scales="free_y")+scale_colour_brewer(palette="Dark2")+
  labs(x=NULL,y="Held-out endpoint RMSE",colour="Model",title="Predictive performance against baseline models",
       caption="Each point is a held-out time block. Vertical scales vary across assemblages.")+
  theme(axis.text.x=element_blank()),"09_model_performance",10,6.5)

miss <- read_result("simulation_missing")
metrics <- rbind(transform(miss,value=relative_eligible_pairs,metric="Supported pair transitions"),
                 transform(miss,value=rmse_all_minus_matched,metric="All–matched discrepancy (RMSE)"),
                 transform(miss,value=endpoint_rmse,metric="Held-out endpoint RMSE"))
save(ggplot(metrics,aes(factor(rate),value,colour=scheme))+
  geom_boxplot(outlier.alpha=.32,position=position_dodge2(width=.8,preserve="single"))+
  facet_wrap(~metric,scales="free_y",ncol=1)+scale_colour_brewer(palette="Dark2")+
  labs(x="Fraction of observations removed",y=NULL,colour="Missingness",title="Sensitivity to incomplete monitoring",
       caption="Simulation panels share the same complete-panel reference within a replicate."),
  "10_incomplete_sampling",8.5,8)
