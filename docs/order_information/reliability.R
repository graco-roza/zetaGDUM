# Reliability is derived only from model fits, independent holds, references and
# repeated fit integrations. Screening covariates are deliberately absent here.
source("docs/order_information/design.R")
read_gz <- function(path) read.csv(gzfile(path))
cal <- do.call(rbind,lapply(list.files(output,"_calibration[.]csv[.]gz$",full.names=TRUE),read_gz))
ref <- do.call(rbind,lapply(list.files(output,"^reference_.*[.]csv$",full.names=TRUE),read.csv))
ref <- ref[ref$converged,]
refkeys <- c("name","order","response")
groups <- split(ref,ref[refkeys],drop=TRUE)
targets <- do.call(rbind,lapply(groups,function(d) {
  data.frame(d[1,refkeys], truth_gradient=mean(d$gradient), truth_secondary=mean(d$secondary),
    truth_difference=mean(d$difference), truth_D2=mean(d$D2),
    reference_gap_sd=sd(d$difference), reference_n=nrow(d))
}))
cal <- merge(cal,targets,by=refkeys,sort=FALSE)
key <- c("name","n","order","split","response")
delta <- .02
summarise <- function(d) {
  truth <- c(gradient=d$truth_gradient[1],secondary=d$truth_secondary[1],difference=d$truth_difference[1])
  valid <- is.finite(d$difference_lo95) & d$converged
  metric <- function(v) {
    observed <- d[[v]]
    lo80 <- d[[paste0(v,"_lo80")]]
    hi80 <- d[[paste0(v,"_hi80")]]
    lo95 <- d[[paste0(v,"_lo95")]]
    hi95 <- d[[paste0(v,"_hi95")]]
    a <- truth[[v]]
    # A failed interval contributes no coverage to the unconditional rate.
    c(bias=mean(observed-a,na.rm=TRUE),rmse=sqrt(mean((observed-a)^2,na.rm=TRUE)),
      coverage80=mean(valid & lo80<=a & hi80>=a,na.rm=FALSE),
      coverage95=mean(valid & lo95<=a & hi95>=a,na.rm=FALSE),
      width95=mean((hi95-lo95)[valid],na.rm=TRUE),
      mc_rmse=sqrt(mean(d[[paste0("mc_",v)]]^2,na.rm=TRUE)),
      null_non_negligible=mean(valid & (lo95>delta | hi95< -delta)))
  }
  effects <- lapply(names(truth),metric)
  names(effects) <- names(truth)
  meaningful <- abs(truth["difference"])>delta
  strong_difference <- abs(truth["difference"])>=2*delta
  direction <- sign(truth["difference"])
  positive <- d$p_meaningful >= .95
  negative <- d$p_reverse_meaningful >= .95
  leading <- if (meaningful) mean(valid & sign(d$difference)==direction) else NA_real_
  power <- if (meaningful) mean(valid & if(direction>0) positive else negative) else NA_real_
  false_superiority <- if (!meaningful) mean(valid & (positive | negative)) else
    mean(valid & if(direction>0) negative else positive)
  enough_prediction <- median(d$heldout,na.rm=TRUE)>=.025 &&
    mean(d$converged & is.finite(d$heldout) & d$heldout>0)>=.75 &&
    mean(d$converged)>=.95
  adequate <- function(m) isTRUE(abs(m["bias"])<=.03 && m["rmse"]<=.06 &&
    m["coverage80"]>=.70 && m["coverage95"]>=.85 && m["width95"]<=.20 &&
    m["mc_rmse"]<=.04)
  checks <- vapply(effects,adequate,logical(1))
  near_null <- names(truth)[names(truth)!="difference" & abs(truth)<=delta]
  meaningful_support <- enough_prediction && mean(valid)>=.90 && all(checks[setdiff(names(checks),near_null)]) &&
    isTRUE(false_superiority<=.10) && (!meaningful || isTRUE(leading>=.85)) &&
    (!strong_difference || isTRUE(power>=.80))
  strict <- meaningful_support && all(checks)
  regime <- if (!enough_prediction) "unsupported" else if (strict) "inferentially_supported" else "prediction_only"
  failure <- if (strict) "none" else if (!enough_prediction) "predictive_transfer_or_fit" else
    if (meaningful_support) "near_null_importance" else "meaningful_precision_or_calibration"
  base <- data.frame(d[1,key], replicates=nrow(d), available=mean(valid),
    converged=mean(d$converged), repeat_converged=mean(d$repeat_converged),
    heldout_median=median(d$heldout,na.rm=TRUE), heldout_positive=mean(d$heldout>0,na.rm=TRUE),
    truth_gradient=d$truth_gradient[1],truth_secondary=d$truth_secondary[1],
    truth_difference=d$truth_difference[1], truth_D2=d$truth_D2[1],
    reference_gap_sd=d$reference_gap_sd[1],reference_n=d$reference_n[1],
    meaningful, leading, superiority_power=power, false_superiority,
    seconds=median(d$seconds), regime, meaningful_support, failure)
  for(v in names(effects)) base[paste0(v,"_",names(effects[[v]]))] <- as.list(effects[[v]])
  base
}
out <- do.call(rbind,lapply(split(cal,cal[key],drop=TRUE),summarise))
write_result(out,file.path(output,"reliability.csv"))
write_result(targets,file.path(output,"population_targets.csv"))
