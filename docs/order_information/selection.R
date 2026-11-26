source("docs/order_information/design.R")
profiles <- do.call(rbind,lapply(list.files(output,"^profile_.*[.]csv[.]gz$",full.names=TRUE),
                               function(f) read.csv(gzfile(f))))
profiles <- profiles[profiles$response %in% c("incidence","relative"),]
profiles$response[profiles$response=="incidence"] <- "total"
reliability <- read.csv(file.path(output,"reliability.csv"))
key <- c("name","n","split","response","order")
data <- merge(profiles,reliability[c(key,"regime")],by=key,sort=FALSE)
data$rcond <- ifelse(is.finite(data$design_condition),1/data$design_condition,0)
data$geometry <- sqrt(data$position_variance*data$spread_variance)
data$break_rule <- with(data,order<=retention_breakpoint)
features <- c("zeta_remaining","retention_rate","response_variance","relative_information","geometry","rcond")
cases <- do.call(rbind,lapply(split(data,data[key],drop=TRUE),function(a) {
  measures <- vapply(a[features],function(x) if(all(is.na(x))) NA_real_ else median(x,na.rm=TRUE), numeric(1))
  data.frame(a[1,c(key,"regime")],as.list(measures),
    breakpoint=if (all(is.na(a$break_rule))) NA else mean(a$break_rule,na.rm=TRUE)>=.5)
}))
rule_types <- c(zeta_decline="zeta_remaining",retention="retention_rate",
                absolute_variance="response_variance",relative_variance="relative_information",
                environmental="geometry")
rates <- function(p,y) {
  c(false_positive=if(any(!y)) mean(p[!y]) else NA_real_,
    excluded=if(any(y)) mean(!p[y]) else NA_real_)
}
grid <- function(x, n=11L) unique(c(-Inf,quantile(x,seq(0,1,length.out=n),na.rm=TRUE,names=FALSE),Inf))
flag <- function(keep) !is.na(keep) & keep
train_models <- function(train, target) {
  truth <- if (target=="inference") train$regime=="inferentially_supported" else train$regime!="unsupported"
  fit <- function(candidates, classify) {
    error <- t(vapply(seq_len(nrow(candidates)),function(i) {
      r <- rates(classify(train,candidates[i,]),truth)
      c(balanced=sum(r),unsupported=r[1])
    },numeric(2)))
    candidates[order(error[,1],error[,2])[1L],,drop=FALSE]
  }
  single <- lapply(rule_types,function(v) {
    a <- fit(expand.grid(limit=grid(train[[v]]),sign=c(1,-1)),
             function(d,row) flag(row$sign*d[[v]]>=row$sign*row$limit))
    list(variable=v, limit=a$limit,sign=a$sign)
  })
  limits <- fit(expand.grid(variance=grid(train$response_variance,5L),
    geometry=grid(train$geometry,5L),rcond=grid(train$rcond,5L)),
    function(d,r) flag(d$response_variance>=r$variance & d$geometry>=r$geometry & d$rcond>=r$rcond))
  c(single,list(retention_breakpoint=list(variable="breakpoint"),combined=as.list(limits)))
}
apply_rule <- function(d,rule) {
  v <- rule$variable
  if (identical(v,"breakpoint")) return(flag(d$breakpoint))
  if (!is.null(v)) return(flag(rule$sign*d[[v]]>=rule$sign*rule$limit))
  with(d, flag(response_variance>=rule$variance & geometry>=rule$geometry & rcond>=rule$rcond))
}
generalise <- function(split,name) {
  ifelse(split=="principal","new_communities",ifelse(split=="validation","unseen_scenarios",
    ifelse(split=="test" & name %in% c("narrow","broad","dominant"),"continuous_familiar",
      ifelse(split=="test","continuous_unseen","training"))))
}
predictions <- do.call(rbind,lapply(c("inference","prediction"),function(target) {
  rules <- train_models(cases[cases$split=="train",],target)
  record <- do.call(rbind,lapply(names(rules),function(r) {
    a <- rules[[r]]
    data.frame(target, rule=r, parameters=paste(paste(names(a),unlist(a),sep="="),collapse=";"))
  }))
  write_result(record,file.path(output,paste0("rules_",target,".csv")))
  do.call(rbind,lapply(names(rules),function(r) {
    data.frame(cases,target,rule=r,predicted=apply_rule(cases,rules[[r]]),
      actual=if(target=="inference") cases$regime=="inferentially_supported" else cases$regime!="unsupported")
  }))
}))
predictions$evaluation <- generalise(predictions$split,predictions$name)
write_result(predictions,file.path(output,"selection_predictions.csv"))
peak <- do.call(rbind,lapply(split(predictions,predictions[c("name","n","response","evaluation","target","rule")],
                                  drop=TRUE),function(a) {
  supported <- sort(a$order[a$actual])
  predicted <- sort(a$order[a$predicted])
  top <- function(x) if(length(x)) max(x) else 1L
  # Nonconsecutive orders only mean holes among orders actually evaluated.
  holes <- function(x) as.integer(length(x)>1L && any(!a$order[a$order>=min(x) & a$order<=max(x)] %in% x))
  data.frame(a[1,c("name","n","response","evaluation","target","rule")], max_grid=max(a$order),
    highest_actual=top(supported),highest_predicted=top(predicted),error=top(predicted)-top(supported),
    supported_orders=paste(supported,collapse=","), predicted_orders=paste(predicted,collapse=","),
    eligible_grid=paste(sort(a$order),collapse=","),
    actual_has_gap=holes(supported),predicted_has_gap=holes(predicted))
}))
write_result(peak,file.path(output,"selection_maxima.csv"))
summary <- do.call(rbind,lapply(split(predictions,predictions[c("evaluation","target","rule")],drop=TRUE),function(a) {
  chosen <- a[a$predicted,]
  rejected <- a[a$actual,]
  errors <- peak$error[peak$rule==a$rule[1] & peak$target==a$target[1] & peak$evaluation==a$evaluation[1]]
  data.frame(a[1,c("evaluation","target","rule")], cells=nrow(a),eligible=nrow(rejected),
    recommended=nrow(chosen),unsupported_recommendations=sum(chosen$regime=="unsupported"),
    fraction_unsupported=if(nrow(chosen)) mean(chosen$regime=="unsupported") else NA_real_,
    fraction_overclaimed=if(nrow(chosen)) mean(!chosen$actual) else NA_real_,
    fraction_excluded=if(nrow(rejected)) mean(!rejected$predicted) else NA_real_,
    mean_max_error=mean(errors),max_order_MAE=mean(abs(errors)))
}))
write_result(summary,file.path(output,"selection_accuracy.csv"))

# Apply frozen screening rules to individual held-out communities, then compare
# with their scenario-level reference labels. Communities in a cell share a label.
per_community <- do.call(rbind,lapply(c("inference","prediction"),function(target) {
  rules <- train_models(cases[cases$split=="train",],target)
  d <- data[data$split!="train",]
  d$breakpoint <- d$break_rule
  do.call(rbind,lapply(names(rules),function(r) {
    data.frame(d[key], target,rule=r,regime=d$regime,predicted=apply_rule(d,rules[[r]]),
      actual=if(target=="inference") d$regime=="inferentially_supported" else d$regime!="unsupported")
  }))
}))
per_community$evaluation <- generalise(per_community$split,per_community$name)
one <- function(a) {
  recommended <- a[a$predicted,]
  data.frame(a[1,c("evaluation","name","n","response","target","rule")],cells=nrow(a),
    inferred=sum(a$actual),recommended=nrow(recommended),
    unsupported_if_recommended=if(nrow(recommended)) mean(recommended$regime=="unsupported") else NA_real_,
    overclaimed_if_recommended=if(nrow(recommended)) mean(!recommended$actual) else NA_real_,
    excluded_if_supported=if(any(a$actual)) mean(!a$predicted[a$actual]) else NA_real_)
}
group <- c("evaluation","name","n","response","target","rule")
write_result(do.call(rbind,lapply(split(per_community,per_community[group],drop=TRUE),one)),
             file.path(output,"selection_community_accuracy.csv"))
