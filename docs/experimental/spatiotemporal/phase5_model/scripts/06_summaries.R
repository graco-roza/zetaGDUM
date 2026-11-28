source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
s <- read_result("empirical_summaries")
a <- subset(s,task=="change" & metric=="incidence")
make <- function(z) {
  z <- z[order(z$initial),]
  cbind(z,cumulative_net=cumsum(z$net),cumulative_negative=cumsum(z$negative),
           cumulative_positive=cumsum(z$positive),
           support_estimand=if(z$dataset[1]=="moths_balanced") "fixed seven-trap panel" else "changing matched panels")
}
cumulative <- bind_rows(lapply(split(a,interaction(a$dataset,a$process,drop=TRUE)),make))
write_result(cumulative,"empirical_cumulative")
P <- subset(read_result("empirical_predictions"),task%in%c("change","local"))
summ <- lapply(split(P,interaction(P$dataset,P$task,drop=TRUE)),function(rows) {
  # Weighted conditional means are calculated within periods before averaging.
  ids <- unique(rows[c("initial","final")])
  ret <- lapply(c("obs","pred"),function(col) {
    pieces <- lapply(split(rows,interaction(rows$initial,rows$final,drop=TRUE)),function(r)
      aggregate(cbind(positive=pmax(r[[col]],0),negative=pmax(-r[[col]],0))~process,r,mean))
    x <- aggregate(cbind(positive,negative)~process,bind_rows(pieces),mean)
    x$net <- x$positive-x$negative
    x$gross <- x$positive+x$negative
    x$positive_share <- 100*.st_ratio(x$positive,sum(x$positive))
    x$negative_share <- 100*.st_ratio(x$negative,sum(x$negative))
    cbind(dataset=rows$dataset[1],task=rows$task[1],series=col,contrasts=nrow(ids),x)
  })
  bind_rows(ret)
})
write_result(bind_rows(summ),"empirical_observed_predicted_summaries")

perf <- read_result("simulation_performance")
I <- read_result("simulation_importance")
ref <- subset(I,scenario=="reference" & fold=="pooled")
test <- subset(I,scenario%in%c("multiple","null") & fold=="pooled")
truth <- aggregate(importance~response+feature,ref,mean)
truth <- bind_rows(lapply(split(truth,truth$response),function(x)
   transform(x,reference_rank=rank(-importance,ties.method="average"))))
rank <- bind_rows(lapply(split(test,interaction(test$scenario,test$response,test$replicate,drop=TRUE)),function(x)
  transform(x,fit_rank=rank(-importance,ties.method="average",na.last="keep"))))
rank <- merge(rank,truth[c("response","feature","reference_rank")],by=c("response","feature"),all.x=TRUE)
write_result(rank,"simulation_rank_comparison")
