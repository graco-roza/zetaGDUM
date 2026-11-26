source("docs/phase5/common.R")
root <- "docs/phase5/results"
stack <- function(pattern, name) {
  files <- list.files("docs/phase5/runs", pattern, full.names=TRUE)
  target <- file.path(root,paste0(name,".csv.gz"))
  if (!length(files)) return(read.csv(target,na.strings=c("NA","")))
  a <- do.call(rbind,lapply(files,read.csv,na.strings=c("NA","")))
  con <- gzfile(target,"wt",compression=9)
  on.exit(close(con))
  write.csv(a,con,row.names=FALSE)
  a
}
summarise <- function(a, key, fn) {
  ix <- split(seq_len(nrow(a)),interaction(a[key],drop=TRUE))
  out <- do.call(rbind,lapply(ix,function(i) cbind(a[i[1L],key,drop=FALSE],
    as.data.frame(as.list(fn(a[i, ,drop=FALSE]))))))
  out[do.call(order,unname(as.list(out[key]))), ,drop=FALSE]
}
save <- function(a,name) write_study(a,root,name)
cal <- stack("^calibration_[0-9]+", "calibration")
di <- stack("^diagnostics_[0-9]+", "posterior_diagnostics")
integ <- stack("^calibration_integration", "integration_sensitivity")
invisible(stack("^diagnostics_integration", "integration_diagnostics"))
gen <- stack("^generating_", "generating_effects")
couple <- stack("^coupling_", "coupling")
paths <- stack("^path_cost_", "path_cpu")
time <- stack("^equal_time", "matched_cpu")
boundary <- stack("^boundary", "active_boundary")
invisible(stack("^retention_stress", "retention_stress"))
adaptive <- stack("^adaptive_confirmation_", "adaptive_confirmation")
invisible(stack("^adaptive_pilots_", "adaptive_pilots"))
emp <- stack("^empirical_points_", "empirical_points")
delete <- stack("^empirical_uncertainty_", "empirical_uncertainty")
lc <- stack("^empirical_local_", "empirical_lczd")
invisible(stack("^empirical_config_", "empirical_configurations"))
status <- stack("^order_status_", "order_status")
ecology <- read.csv("docs/phase4/results/empirical_profiles.csv")
support_names <- c("community","order","expected_richness","richness_variance","sharing_lower","sharing_upper",
                   "sharing_estimate","abundance_zeta","effective_incidence_species","effective_abundance_species",
                   "incidence_retention","abundance_retention","retention_break")
save(merge(status,ecology[ecology$component=="total",support_names],
           by=c("community","order"),all.x=TRUE),"order_limits_with_ecological_support")
invisible(stack("^budget_[NS]", "empirical_budgets"))
invisible(stack("^budget_pilot_", "empirical_budget_pilots"))

cal$valid <- is.finite(cal$lo) & is.finite(cal$hi)
cal$covered <- ifelse(cal$valid, cal$truth >= cal$lo & cal$truth <= cal$hi, NA)
cal$decision <- ifelse(cal$method=="jackknife", cal$lo>0, cal$p_positive >= .95)
cal$kind <- ifelse(grepl("_change$",cal$contrast),"change","level")
caltab <- summarise(cal,c("n","signal","method","contrast","level"),function(x) c(
  replicates=nrow(x),available=sum(x$valid),coverage=mean(x$covered,na.rm=TRUE),
  posterior_mean_bias=mean(x$posterior_mean-x$truth,na.rm=TRUE),
  point_bias=mean(x$point-x$truth),prob_positive=mean(x$p_positive,na.rm=TRUE),
  decision=mean(x$decision,na.rm=TRUE),prob_practical=mean(x$p_practical,na.rm=TRUE)))
save(caltab,"calibration_by_cell")
dtab <- summarise(di,c("n","signal","method"),function(x) c(reps=nrow(x),unavailable=sum(!is.na(x$reason)),
  ESS=median(x$ess,na.rm=TRUE),maximum_invalid_proposal_share=max(c(0,x$rejected),na.rm=TRUE),
  seconds=median(x$seconds,na.rm=TRUE)))
save(dtab,"posterior_diagnostics_by_cell")
fixed <- subset(cal,signal %in% c("null","weak","strong"))
pooled <- summarise(fixed,c("method","level"),function(x) c(
  coverage=mean(x$covered,na.rm=TRUE),availability=mean(x$valid)))
save(pooled,"coverage_summary")
num <- function(a) c(available=sum(is.finite(a$estimate)),
  RMSE=sqrt(mean((a$estimate-a$reference)^2,na.rm=TRUE)),
  variance=var(a$estimate,na.rm=TRUE),median_CPU=median(a$cpu))
coup <- summarise(couple,c("case","component","budget","sampler","coupled","comparison","effect"),num)
save(coup,"coupling_summary")
save(summarise(time,c("sampler","coupled","M","effect"),num),"matched_cpu_summary")
save(summarise(gen,c("n","case","model","effect","order"),function(x) c(
  RMSE=sqrt(mean((x$point-x$truth)^2,na.rm=TRUE)),
  bias=mean(x$point-x$truth,na.rm=TRUE),
  scaled_condition=median(x$condition,na.rm=TRUE),converged=mean(x$converged),
  prediction_RMSE=mean(x$pred_RMSE))),"generating_summary")
stopifnot(max(abs(lc$observed-lc$baseline-lc$spread-lc$position-lc$residual),na.rm=TRUE) < 1e-10)

color <- c(cp="#ad5b32",jackknife="#414b56",pseudo="#076f70",jel="#5869a8",etel="#9480ad")
labels <- c("Product-weight bootstrap","Site jackknife","Centred pseudovalue bootstrap","Joint Bayesian JEL","Joint ETEL")
png(file.path(root,"posterior_calibration.png"),1450,780,res=145)
par(mfrow=c(1,2),mar=c(4.2,4.2,2.6,1),mgp=c(2.65,.7,0))
plot(c(.49,.96),c(.45,1.01),type="n",xlab="Nominal interval probability",ylab="Observed coverage",
  main="Fixed generating parameters")
abline(0,1,col="#b0b0b0",lty=2)
for (method in names(color)) {
  a <- pooled[pooled$method==method,]
  lines(a$level,a$coverage,type="b",col=color[method],pch=16,lwd=2)
}
legend("bottomright",labels,col=color,pch=16,bty="n",cex=.8)
a <- unique(subset(cal,signal=="prior" & level==.95,c(method,replicate,contrast,posterior_cdf_truth)))
plot(c(0,1),c(0,1),type="n",xlab="Posterior CDF at the generating effect",
  ylab="Empirical cumulative fraction",main="Prior-predictive check")
abline(0,1,col="#b0b0b0",lty=2)
for (method in intersect(names(color),unique(a$method))) {
  z <- sort(na.omit(a$posterior_cdf_truth[a$method==method]))
  if(length(z)) lines(z,seq_along(z)/length(z),type="s",col=color[method],lwd=2)
}
invisible(dev.off())

a <- subset(coup,budget==2400 & comparison=="change")
key <- c("case","component","sampler","effect")
q <- merge(a[a$coupled,],a[!a$coupled,],by=key,suffixes=c("_paired","_independent"))
q$ratio <- q$RMSE_paired/q$RMSE_independent
order <- c("spread","position","importance","shapley_spread","shapley_position",
           "relative_position","local_spread","local_position","local_predicted","local_residual")
label <- c("Spread","Position","Predictor importance","Spread Shapley","Position Shapley",
           "Position Shapley share","Local spread","Local position","Local prediction","Local residual")
png(file.path(root,"paired_numerical_precision.png"),1170,910,res=145)
par(mar=c(4.7,10.7,3.4,1),mgp=c(2.85,.75,0))
plot(c(.012,9),c(.4,10.7),type="n",log="x",yaxt="n",ylab="",
  xlab="Paired change RMSE: coupled / independent",
  main="Coupling helps some summaries and hurts others")
axis(2,at=seq_along(label),labels=label,las=1)
abline(v=1,col="#626262",lty=2)
cols <- c(site="#076f70",species="#ad5b32",uniform="#5869a8")
for (j in seq_along(order)) {
  a <- q[q$effect==order[j]&is.finite(q$ratio)&q$ratio>0,]
  y <- j + .082*(match(a$sampler,names(cols))-2)
  points(a$ratio,y,pch=16,col=adjustcolor(cols[a$sampler],.65),cex=.83)
}
legend("topright",c("Site-balanced","Species-informed","Uniform"),col=cols,
  pch=16,bty="n",cex=.8)
invisible(dev.off())

a <- subset(delete,component=="relative" & is.na(previous) & effect=="position")
a <- do.call(rbind,lapply(split(a,a$community),function(x)x[x$order==max(x$order),]))
community <- sort(unique(a$community))
png(file.path(root,"empirical_site_uncertainty.png"),1290,780,res=145)
par(mar=c(4.8,10,3.3,1),mgp=c(2.75,.75,0))
plot(range(c(a$jack_lo,a$jack_hi),finite=TRUE),c(.35,length(community)+.6),type="n",yaxt="n",
  xlab="One-site IQR position effect (logit)",ylab="",
  main="Site-deletion intervals at the middle manuscript order")
axis(2,seq_along(community),labels=community,las=1,cex.axis=.85)
abline(v=0,lty=2,col="#a7a7a7")
cc <- c(A="#454f66",B="#087b79",C="#ac592f")
for(i in seq_len(nrow(a))) {
  j <- match(a$community[i],community) + .20*(match(a$model[i],names(cc))-2)
  col <- cc[a$model[i]]
  segments(a$jack_lo[i],j,a$jack_hi[i],j,col=adjustcolor(col,.55),lwd=1)
  segments(a$pseudo_lo[i],j,a$pseudo_hi[i],j,col=col,lwd=2.7)
  points(a$point[i],j,col=col,pch=16)
}
legend("bottomright",c("A","B","C"),col=cc,pch=16,bty="n",horiz=TRUE)
invisible(dev.off())
