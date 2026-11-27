source("docs/final_gdum/design.R")
read_gz <- function(f) read.csv(gzfile(f))
files <- list.files(outdir,"_(pilot|principal)[.]csv[.]gz$",full.names=TRUE)
sim <- do.call(rbind,lapply(files,function(f) {
  a <- read_gz(f)
  if (!"calibration_slope" %in% names(a)) a$calibration_slope <- NA_real_
  data.frame(phase=if(grepl("_pilot",f)) "pilot" else "principal", a)
}))
references <- do.call(rbind,lapply(file.path(outdir,paste0("reference_",mechanisms$name,".csv")),read.csv))
keys <- c("scenario","order","model","architecture","response","predictor")
valid <- references[is.finite(references$absolute) & references$converged,]
ref <- aggregate(absolute~scenario+order+model+architecture+response+predictor+reference,valid,sum)
truth <- do.call(rbind,lapply(split(ref,ref[keys],drop=TRUE),function(x) {
  data.frame(x[1,keys], target=mean(x$absolute), target_range=diff(range(x$absolute)),
    reference_replicates=nrow(x))
}))
atomic_csv(truth,"reference_importance.csv")
sim$response <- sim$component
a <- merge(sim,truth,by=keys,all.x=TRUE)
key <- c(keys,"phase","n")
summary <- function(x) {
  med <- function(v) if(all(is.na(v))) NA_real_ else median(v,na.rm=TRUE)
  mean_num <- function(v) if(all(is.na(v))) NA_real_ else mean(v,na.rm=TRUE)
  delta <- x$absolute-x$target
  data.frame(x[1,key], replicates=nrow(x),convergence=mean(x$converged),
    valid_importance=mean(is.finite(x$absolute)), positive_D2=mean_num(x$D2_observed>0),
    median_D2=med(x$D2_observed),median_oracle_D2=med(x$D2_oracle),observed_RMSE=med(x$RMSE),
    oracle_RMSE=med(x$oracle_RMSE),mean_calibration_bias=mean_num(x$bias),
    calibration_slope=med(x$calibration_slope), closure_RMSE=med(x$closure_RMSE),
    numerical_prediction_RMSE=med(x$numerical_prediction_RMSE),
    rank_deficient=mean(x$min_rank<x$parameters),condition=med(x$max_condition),
    min_interaction_residual=med(x$min_interaction_residual),
    target=med(x$target),target_reference_range=med(x$target_range),
    importance_bias=mean_num(delta),importance_RMSE=sqrt(mean_num(delta^2)),
    median_importance=med(x$absolute),fit_seconds=med(x$fit_seconds),total_seconds=med(x$total_seconds))
}
res <- do.call(rbind,lapply(split(a,a[key],drop=TRUE),summary))
atomic_csv(res,"model_summary.csv")
keys_pair <- c("phase","scenario","n","order","model","architecture","response","replicate")
q <- a[a$predictor=="phosphorus",c(keys_pair,"absolute","target")]
r <- a[a$predictor=="secondary",c(keys_pair,"absolute","target")]
match <- merge(q,r,by=keys_pair)
key_pair <- setdiff(keys_pair,"replicate")
unit <- lapply(split(match,match[key_pair],drop=TRUE),function(x) {
  h <- x$absolute.x-x$absolute.y
  truth <- x$target.x-x$target.y
  ok <- is.finite(h) & is.finite(truth)
  leading <- ok & abs(truth)>.02
  data.frame(x[1,key_pair],true_difference=mean(truth[ok]),recovered=if(any(leading))
    mean(sign(h[leading])==sign(truth[leading])) else NA_real_,
    difference_bias=mean((h-truth)[ok]),difference_RMSE=sqrt(mean((h-truth)[ok]^2)),
    n_valid=sum(ok),n_meaningful=sum(leading))
})
atomic_csv(do.call(rbind,unit),"contrast_summary.csv")
files <- list.files(outdir,"^interval_.*[.]csv[.]gz$",full.names=TRUE)
if(length(files)) {
 u <- do.call(rbind,lapply(files,function(f) {
  a <- read_gz(f)
  if (!"assessment" %in% names(a)) a$assessment <- "surface"
  a
 }))
 spec <- c("scenario","order","mode","assessment","estimand")
 interval <- function(x) {
  m <- function(v) mean(v,na.rm=TRUE)
  coverage <- function(lo,hi,target) mean(is.finite(lo) & lo<=target & hi>=target, na.rm=TRUE)
  valid <- is.finite(x$lo95)
  null_truth <- abs(x$truth)<=.02
  data.frame(x[1,spec],replicates=nrow(x),interval_available=mean(valid),
   draws=x$draws[1],truth=x$truth[1],mean_active=m(x$active),
   point_bias=m(x$point-x$truth),finite_bias=m(x$finite_target-x$truth),
   integration_bias=m(x$point-x$finite_target),
   integration_RMSE=sqrt(m((x$point-x$repeat_point)^2)/2),
   bootstrap_center_bias=m(x$bootstrap_mean-x$point),
   point_RMSE=sqrt(m((x$point-x$truth)^2)),corrected_bias=m(x$corrected-x$truth),
   corrected_RMSE=sqrt(m((x$corrected-x$truth)^2)),
   percentile80=coverage(x$lo80,x$hi80,x$truth),
   percentile95=coverage(x$lo95,x$hi95,x$truth),
   basic80=coverage(x$basic_lo80,x$basic_hi80,x$truth),
   basic95=coverage(x$basic_lo95,x$basic_hi95,x$truth),
   percentile95_finite=coverage(x$lo95,x$hi95,x$finite_target),
   basic95_finite=coverage(x$basic_lo95,x$basic_hi95,x$finite_target),
   width95=m(x$hi95-x$lo95),false_percentage_declaration=if(all(null_truth))
    mean(valid & (x$lo95>.02 | x$hi95< - .02)) else NA_real_,
   false_basic_declaration=if(all(null_truth))
    mean(valid & (x$basic_lo95>.02 | x$basic_hi95< - .02)) else NA_real_)
 }
 results <- do.call(rbind,lapply(split(u,u[spec],drop=TRUE),interval))
 atomic_csv(results,"interval_summary.csv")
}
emp <- do.call(rbind,lapply(list.files(outdir,"^emp_validation_.*[.]csv$",full.names=TRUE),read.csv))
if(nrow(emp)) {
 ekey <- c("id","community","order","model","architecture","component")
 v <- do.call(rbind,lapply(split(emp,emp[ekey],drop=TRUE),function(x) {
  data.frame(x[1,ekey],folds=nrow(x),converged=sum(x$converged),
    pooled_D2=1-sum(x$deviance)/sum(x$null_deviance), positive_folds=sum(x$D2>0,na.rm=TRUE),
    median_RMSE=median(x$RMSE), mean_bias=mean(x$bias),support=mean(x$supported_fraction),
    closure_RMSE=median(x$closure_RMSE))
 }))
 atomic_csv(v,"empirical_summary.csv")
}
