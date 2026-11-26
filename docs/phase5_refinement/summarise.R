source("docs/phase5_refinement/common.R")
root <- "docs/phase5_refinement/results"
stack <- function(pattern, name, optional=FALSE) {
  paths <- list.files("docs/phase5_refinement/runs",pattern,full.names=TRUE)
  out <- file.path(root,paste0(name,".csv.gz"))
  if (!length(paths)) {
    if (file.exists(out)) return(read.csv(out,na.strings=c("","NA")))
    if (optional) return(NULL)
    stop("Missing run data for ",name)
  }
  a <- do.call(rbind,lapply(paths,read.csv,na.strings=c("","NA")))
  con <- gzfile(out,"wt",compression=9)
  on.exit(close(con))
  write.csv(a,con,row.names=FALSE)
  a
}
tabulate_by <- function(a, key, fn) {
  g <- split(seq_len(nrow(a)),interaction(a[key],drop=TRUE))
  z <- do.call(rbind,lapply(g,function(i)cbind(a[i[1L],key,drop=FALSE],
    as.data.frame(as.list(fn(a[i,,drop=FALSE]))))))
  z[do.call(order,unname(as.list(z[key]))),,drop=FALSE]
}
save <- function(a,name) write_study(a,root,name)
rmse <- function(x) sqrt(mean(x^2,na.rm=TRUE))
noise <- stack("^noise_", "noise")
cover <- stack("^coverage_", "coverage")
exact <- stack("^exact_[0-9]", "exact_paired")
exact_low <- stack("^exact_lower", "exact_lower")
curve <- stack("^curvature_", "curvature")
samplers <- stack("^equal_cpu_", "equal_cpu")
budget <- stack("^emp_budget_", "empirical_budget",TRUE)
empirical <- stack("^emp_uncert_", "empirical_uncertainty",TRUE)
posterior <- stack("^posterior_(calibration|reference|mixture)_.*_[0-9]+-[0-9]+|^posterior_long_reference_", "posterior",TRUE)
diag_post <- stack("^postdiag_(calibration|reference|mixture)_.*_[0-9]+-[0-9]+|^postdiag_long_reference_", "posterior_diagnostics",TRUE)
invisible(stack("^hull_.*_[0-9]+-[0-9]+", "hull_stress",TRUE))
invisible(stack("^kernel_timing", "kernel_timing",TRUE))

n <- tabulate_by(noise,c("n","signal","effect","method","sampler","M"),function(x)c(
  replicates=nrow(x),median_fraction=median(x$mc_fraction,na.rm=TRUE),
  median_MC=median(x$MC_var,na.rm=TRUE),median_JK=median(x$JK_var,na.rm=TRUE),
  median_pooled_JK=median(x$pooled_var,na.rm=TRUE),cpu=median(x$cpu,na.rm=TRUE)))
save(n,"noise_summary")
cal <- function(a)c(replicates=nrow(a),coverage=mean(a$covered),width=median(a$width),
                    bias=mean(a$estimate-a$truth),RMSE=rmse(a$estimate-a$truth))
save(tabulate_by(cover,c("n","signal","sampler","M","method","effect","level"),cal),
     "coverage_summary")
save(tabulate_by(exact_low,c("n","signal","method","effect","level"),cal),"exact_lower_summary")
save(tabulate_by(exact,c("n","signal","sampler","M","effect"),function(a)c(
  finite_point_RMSE=rmse(a$point_error),median_pseudo_error=median(a$pseudo_RMSE),
  median_pooled_error=median(a$pooled_pseudo_RMSE))),"exact_paired_summary")
save(tabulate_by(curve,c("n","kind","model","order","effect"),function(a)c(
  bias=mean(a$estimate-a$truth,na.rm=TRUE),error_SD=sd(a$estimate-a$truth,na.rm=TRUE),
  estimate_SD=sd(a$estimate,na.rm=TRUE),RMSE=rmse(a$estimate-a$truth),
  curvature=median(a$curvature),convergence=mean(a$converged))),"curvature_summary")
save(tabulate_by(samplers,c("case","sampler","M","effect"),function(a)c(
  spread_RMSE=rmse(a$spread_error),position_RMSE=rmse(a$position_error),
  median_score_error=median(a$score_error_ratio),median_CPU=median(a$CPU),
  target_CPU=median(a$target_CPU),high_sharing_lower=a$high_sharing_lower[1],
  high_sharing_upper=a$high_sharing_upper[1])),"equal_cpu_summary")

if (!is.null(posterior)) {
  p <- posterior[posterior$mode=="calibration",]
  if (nrow(p)) {
  for (q in c(50,80,90,95)) {
    a <- p
    a$covered <- a[[paste0("l",q)]]<=a$truth & a[[paste0("h",q)]]>=a$truth
    a$level <- q/100
    if (q==50) wide <- a else wide <- rbind(wide,a)
  }
  ptab <- tabulate_by(wide,c("signal","M","method","effect","level"),function(a) {
    use <- which(a$stable %in% TRUE)
    total <- with(diag_post, sum(mode=="calibration" & signal==a$signal[1] &
                                  M==a$M[1] & method==a$method[1]))
    c(total=total,returned=nrow(a),stable=length(use),coverage=mean(a$covered[use]),
      bias=mean(a$mean[use]-a$truth[use]),power=mean(a$p[use]>=.95),
      practical_decision=mean(a$practical[use]>=.95))
  })
  save(ptab,"posterior_coverage_stable")
  }
  sample <- subset(posterior,mode %in% c("reference","mixture","long_reference"))
  main <- subset(sample,method=="mcmc" & stable)
  keys <- c("signal","r","M","effect")
  paired <- merge(sample[sample$method!="mcmc",],main[,c(keys,"p","practical","mean","sd","l95","h95")],
    by=keys,suffixes=c("","_mcmc"))
  paired$p_error <- abs(paired$p-paired$p_mcmc)
  paired$practical_error <- abs(paired$practical-paired$practical_mcmc)
  paired$mean_diff_sd <- abs(paired$mean-paired$mean_mcmc)/paired$sd_mcmc
  paired$quantile_diff_sd <- pmax(abs(paired$l95-paired$l95_mcmc),abs(paired$h95-paired$h95_mcmc))/paired$sd_mcmc
  save(paired,"posterior_reference_comparison")
}

cols <- c(uniform="#356182",site="#a64c30")
png(file.path(root,"site_pseudovalue_noise.png"),1430,670,res=145)
par(mfrow=c(1,3),mar=c(4.1,4.45,2.7,.9),mgp=c(2.7,.75,0))
p <- subset(n,n==70 & signal=="weak" & method=="parameter")
for (type in c("lower","higher","paired")) {
  plot(c(900,50000),c(0,1),type="n",log="x",xlab="Combinations per pool",
    ylab="Numerical fraction of jackknife variance",main=paste(tools::toTitleCase(type),"contrasts"))
  for (s in names(cols)) for (e in c("spread","position")) {
    a <- p[p$sampler==s & p$effect==paste(type,e,sep="_"),]
    lines(a$M,a$median_fraction,type="b",lwd=2,pch=if(e=="spread") 16 else 17,
      lty=if(e=="spread") 1 else 2,col=cols[s])
  }
  if (type=="lower") legend("topright",c("Uniform","Site-balanced"),pch=16,col=cols,bty="n",cex=.82)
  if (type=="higher") legend("topright",c("Spread","Position"),lty=1:2,pch=c(16,17),bty="n",cex=.82)
}
invisible(dev.off())

if (!is.null(posterior)) {
  png(file.path(root,"jel_importance_efficiency.png"),1290,720,res=145)
  par(mfrow=c(1,2),mar=c(4.6,4.75,2.8,1),mgp=c(3,.8,0))
  e <- subset(diag_post, mode %in% c("reference","mixture","long_reference") & method!="mcmc" & is.finite(ess))
  names_method <- c(adaptive="Modal",chain_mixture="Chain mixture",legacy="Legacy")
  if (nrow(e)) {
    e$method <- names_method[e$method]
    boxplot(ess~method,data=e,col=c("#d5e2f0","#c4e2dc","#e7d8cd"),log="y",las=1,
      main="Importance ESS",ylab="Effective weighted draws")
  } else {
    plot.new()
    title("Importance comparison unavailable")
  }
  if (nrow(paired)) {
    v <- paired[is.finite(paired$p_error),]
    v$method <- names_method[v$method]
    boxplot(p_error~method,data=v,col=c("#d5e2f0","#c4e2dc","#e7d8cd"),las=1,
      main="Agreement with eligible MCMC runs",ylab="Absolute error in positive probability")
  } else {
    plot.new()
    title("No eligible MCMC comparison")
  }
  invisible(dev.off())
}
