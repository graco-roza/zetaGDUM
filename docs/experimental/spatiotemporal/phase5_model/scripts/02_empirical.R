source("docs/experimental/spatiotemporal/phase5_model/scripts/common.R")
moth_root <- file.path(base,"phase4b_annual_metacommunities")
lake_root <- file.path(base,"phase4_swedish_lakes")
source(file.path(phase,"scripts/archives.R"))
lake <- st5_lake_archive(lake_root)
lake_inputs <- lake$inputs
spatial_lake <- lake$blocks
write_result(bind_rows(lapply(names(lake_inputs),function(n) cbind(dataset=n,lake_inputs[[n]]$audit))),"lake_environment_audit")
p_moth <- st5_moth_archive(moth_root)
obs <- p_moth$records
x <- p_moth$x
moth_previous <- readRDS(file.path(moth_root,"data/processed/biotime-604.rds"))$panel
prev <- obs$site%in%dimnames(moth_previous$x)[[1]] & obs$time%in%moth_previous$times
lake_previous <- readRDS(file.path(lake_root,"data/processed/panels.rds"))
inputs <- c(list(moths=list(p=p_moth),
                 moths_balanced=list(p=st5_input(obs[prev,],incidence=x[prev,,drop=FALSE]))),lake_inputs)
run_dataset <- function(name) {
  coverage <- trajectory <- transitions <- summaries <- list()
  models <- importance <- predictions <- curves <- list()
  observations <- list()
  support_delta <- strict <- decomposed <- list()
  tables <- c("coverage","trajectory","transitions","summaries","models","importance","predictions",
              "curves","support_delta","strict","observations")
  finish <- function(env=parent.frame()) { force(env); setNames(lapply(tables,function(n) bind_rows(get(n,envir=env))),tables) }
  message(name)
  p <- inputs[[name]]$p
  exact <- grepl("^moths",name)
  observations[[name]] <- cbind(dataset=name,p$records[c("site","time")])
  settings <- list(
    spatial=st5_support(p,"spatial",budget=Inf),
    local=st5_support(p,"local",lag=1,budget=Inf),
    local_observed=st5_support(p,"local",budget=Inf),
    rectangle=st5_support(p,"rectangle",m=2,budget=160L),
    change=st5_support(p,"change",lag=1,budget=if(exact) Inf else 240L),
    change_joint=st5_support(p,"change",m=2,budget=96L),
    shared=st5_support(p,"shared",k=2,m=1,lag=1,budget=if(exact) Inf else 120L))
  for(task in names(settings)) {
    s <- settings[[task]]
    coverage[[paste(name,task)]] <- cbind(dataset=name,task,observed=nrow(p$x),st5_eligibility(p,s))
  }
  if(name=="moths") stopifnot(sum(vapply(settings$spatial,function(g)nrow(g$sets),0L))==1551L,
     sum(vapply(settings$local,function(g)nrow(g$sets),0L))==253L,
     sum(vapply(settings$change,function(g)nrow(g$sets),0L))==1293L)
  spatial <- st5_responses(p,settings$spatial)
  area <- if(!exact) {
    member <- lake_previous$membership
    preferred <- member$id[member$assemblage==name & member$design=="B" & member$included]
    years <- lake_previous$years[["B"]]
    if(length(preferred)>1L && length(years)>0L) {
      sub <- p$records$site%in%preferred & p$records$time%in%years
      q <- st5_input(p$records[sub,,drop=FALSE],abundance=p$a[sub,,drop=FALSE])
      st5_responses(q,st5_support(q,"spatial",budget=Inf))
    }
  }
  stat <- function(z) bind_rows(lapply(colnames(z$before),function(metric) {
    val <- tapply(z$before[,metric],z$rows$initial,mean,na.rm=TRUE)
    data.frame(time=as.numeric(names(val)),metric,value=as.numeric(val))
  }))
  trajectory[[name]] <- cbind(dataset=name,panel="response specific",stat(spatial))
  if(!is.null(area)) trajectory[[paste0(name,"_balanced")]] <- cbind(dataset=name,panel="Phase IV balanced",stat(area))
  for(task in setdiff(names(settings),"spatial")) {
    paired <- task %in% c("local","change","change_joint","shared")
    z <- st5_responses(p,settings[[task]],paired,draws=4L)
    if(!nrow(z$rows)) next
    groups <- split(seq_len(nrow(z$rows)),z$rows$period)
    for(metric in colnames(z$before)) {
      a <- st5_extract(z,metric)
      key <- paste(name,task,metric)
      transitions[[key]] <- bind_rows(lapply(groups,function(i) {
        i1 <- i[is.finite(a$y[i]) & is.finite(a$baseline[i])]
        data.frame(dataset=name,task,metric,initial=a$rows$initial[i[1]],final=a$rows$final[i[1]],
            lag=a$rows$lag[i[1]],eligible=a$rows$eligible[i[1]],sampled=length(i),
            defined=length(i1),before=mean(a$baseline[i1]),endpoint=mean(a$y[i1]),
            delta=mean(a$y[i1]-a$baseline[i1]))
      }))
      if(!is.null(a$phi)) summaries[[key]] <- bind_rows(lapply(groups,function(i) {
        ax <- st5_process_summary(a$phi[i,,drop=FALSE])
        ax$defined <- sum(complete.cases(a$phi[i,,drop=FALSE]))
        cbind(dataset=name,task,metric,initial=a$rows$initial[i[1]],final=a$rows$final[i[1]],ax)
      }))
    }
    if(task=="change") {
      tt <- p$times
      dates <- cbind(tt[-length(tt)],tt[-1])
      support_delta[[name]] <- cbind(dataset=name,st5_support_change(spatial,z,"incidence",dates))
    }
    if(!exact && task%in%c("local","change","change_joint","shared")) {
      metric <- switch(task,local="site_temporal_presence",shared="temporal_presence","incidence")
      decomposed[[task]] <- st5_extract(z,metric)
      if(task=="change") decomposed$crossed <- st5_crossed_data(z)
    }
  }
  if(exact) return(finish())
  # Static fitting and blocked validation use a predeclared uniform per-year budget.
  decomposed$spatial <- st5_extract(st5_responses(p,st5_support(p,"spatial",budget=240)),"incidence")
  decomposed$rectangle <- st5_extract(st5_responses(p,settings$rectangle),"incidence")
  for(task in names(decomposed)) {
    z <- decomposed[[task]]
    main <- if(task=="crossed") "crossed" else if(is.null(z$phi)) "fractional" else "coherent"
    key <- paste(name,task)
    methods <- if(task%in%c("change","local","spatial")) c("spatial","temporal","forward") else "temporal"
    for(method in methods) {
      message("  ",task," / ",method)
      # Four blocks allow forward folds to retain at least two inner time periods.
      result <- st5_validate(p,z,method,spatial_lake,mode=main,blocks=4L)
      if(!nrow(result$scores)) next
      models[[paste(key,method,main)]] <- cbind(dataset=name,task,result$scores)
      if(main=="coherent" && task=="change") {
        r <- st5_validate(p,z,method,spatial_lake,mode="ridge",blocks=4L)
        if(nrow(r$scores)) models[[paste(key,method,"ridge")]] <- cbind(dataset=name,task,r$scores)
      }
      if(method=="temporal" || (method=="spatial" && task=="change"))
        importance[[paste(key,method)]] <- cbind(dataset=name,task,method,bind_rows(list(
          st5_importance(result,"predictor"),st5_importance(result,"group"))))
      if(method=="temporal") {
        if(task=="spatial") {
          mono <- st5_validate(p,z,method,spatial_lake,mode=main,blocks=4L,monotonic=TRUE)
          if(nrow(mono$scores)) models[[paste(key,method,"monotonic")]] <- cbind(dataset=name,task="spatial_monotonic",mono$scores)
        }
        if(task=="local" || task=="change") {
          plot_pred <- lapply(result$folds,function(d) {
            fit <- .st5_refit(d)
            pred <- st5_predict(fit,d$te,d$b)
            i <- d$j
            rows <- d$te$rows[i,]
            do.call(rbind,lapply(seq_len(ncol(pred$phi)),function(j)
              data.frame(dataset=name,task,fold=d$label,initial=rows$initial,final=rows$final,
                process=colnames(pred$phi)[j],obs=d$te$phi[i,j],pred=pred$phi[i,j])))
          })
          predictions[[key]] <- bind_rows(plot_pred)
        }
        if(task=="change") {
          eq <- st5_validate(p,z,method,spatial_lake,mode=main,blocks=4L,equal_rows=TRUE)
          models[[paste(key,method,"equal_rows")]] <- cbind(dataset=name,task="change_equal_rows",eq$scores)
          for(d in result$folds) {
            fit <- .st5_refit(d)
            for(v in d$ref$vars) {
              target <- paste0(v,":dP")
              if(!target %in% colnames(d$a$x)) next
              median_row <- apply(d$a$x[d$i,,drop=FALSE],2,median)
              xx <- seq(quantile(d$a$x[d$i,target],.05),quantile(d$a$x[d$i,target],.95),length.out=31)
              design <- d$a
              design$x <- matrix(rep(median_row,each=length(xx)),nrow=length(xx),dimnames=list(NULL,names(median_row)))
              design$x[,target] <- xx
              design$missing <- design$extrapolation <- rep(FALSE,length(xx))
              pred <- st5_predict(fit,list(baseline=design$x[,"baseline"]),design)
              curves[[paste(key,d$label,v)]] <- bind_rows(lapply(seq_len(ncol(pred$phi)),function(j)
                data.frame(dataset=name,fold=d$label,predictor=v,process=colnames(pred$phi)[j],
                  change_sd_per_site=xx/unique(z$rows$k),contribution=pred$phi[,j])))
            }
          }
        }
      }
    }
  }
  # Reclassify trajectories inside the held-out region to compare the two estimands.
  for(block in sort(unique(spatial_lake))) {
    old <- decomposed$change
    old$rows$key <- vapply(old$rows$sites,function(s) paste(s,collapse="/"),character(1))
    keep <- p$records$site%in%names(spatial_lake)[spatial_lake==block]
    q <- st5_input(p$records[keep,],incidence=p$x[keep,,drop=FALSE])
    intra <- st5_extract(st5_responses(q,st5_support(q,"change",lag=1,budget=Inf),TRUE,draws=4L),"incidence")
    if(!nrow(intra$rows)) next
    ik <- vapply(intra$rows$sites,function(s) paste(s,collapse="/"),character(1))
    ii <- match(paste(old$rows$initial,old$rows$final,old$rows$key),paste(intra$rows$initial,intra$rows$final,ik))
    jj <- which(!is.na(ii))
    if(!length(jj)) next
    di <- intra$phi[ii[jj],,drop=FALSE]-old$phi[jj,1:3,drop=FALSE]
    strict[[paste(name,block)]] <- data.frame(dataset=name,block,rows=nrow(di),process=colnames(di),
      mean_reclassified_minus_regional=colMeans(di,na.rm=TRUE),mean_absolute_difference=colMeans(abs(di),na.rm=TRUE))
  }
  finish()
}
workers <- as.integer(Sys.getenv("ST5_EMPIRICAL_WORKERS","3"))
res <- if(workers>1L) parallel::mclapply(names(inputs),run_dataset,mc.cores=workers,mc.preschedule=FALSE) else
  lapply(names(inputs),run_dataset)
if(any(vapply(res,inherits,logical(1),"try-error"))) stop(res[[which(vapply(res,inherits,logical(1),"try-error"))[1]]])
for(n in names(res[[1]])) write_result(bind_rows(lapply(res,`[[`,n)),paste0("empirical_",n))
