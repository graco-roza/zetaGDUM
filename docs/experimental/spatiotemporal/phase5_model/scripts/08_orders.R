source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
source(file.path(phase,"scripts/archives.R"))
lake <- st5_lake_archive(file.path(base,"phase4_swedish_lakes"))
includes <- c(list(moths=st5_moth_archive(file.path(base,"phase4b_annual_metacommunities"))),
              lapply(lake$inputs,function(z) st5_input(z$p$records,incidence=z$p$x)))
collect <- function(p,dataset) {
  value <- list()
  for(k in 2:6) for(m in 1:5) {
    static <- st5_support(p,if(m==1L) "spatial" else "rectangle",k,m,budget=64L)
    sp <- st5_responses(p,static)
    sh <- st5_support(p,"shared",k,m,budget=40L)
    share <- st5_responses(p,sh,TRUE,draws=4)
    matched <- st5_support(p,"change",k,m,budget=40L)
    change <- st5_responses(p,matched,TRUE,draws=4)
    incidence <- sp$before[,"incidence"]
    zeta <- sp$before[,"zeta_presence"]
    mean_defined <- function(x) if(length(x)) mean(is.finite(x)) else NA_real_
    phi <- function(x,metric) if(is.null(x$phi)) NA_real_ else mean(is.finite(x$phi[,metric,1]))
    value[[paste(k,m)]] <- data.frame(dataset,k,m,static_rows=nrow(sp$rows),
      mean_shared=mean(zeta),zero_sharing=mean(zeta==0),static_incidence=mean_defined(incidence),
      static_turnover=mean_defined(sp$before[,"turnover"]),static_sd=sd(incidence,na.rm=TRUE),
      mean_incidence=mean(incidence,na.rm=TRUE),matched_rows=nrow(change$rows),
      joint_change_support=phi(change,"incidence"),shared_reorganisation_support=phi(share,"temporal_presence"),
      eligible_observation_fraction=st5_eligibility(p,static)$eligible_observation_fraction)
  }
  bind_rows(value)
}
results <- lapply(names(includes),function(d) collect(includes[[d]],d))
write_result(bind_rows(results),"empirical_order_support")
