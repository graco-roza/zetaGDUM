source("docs/final_gdum/design.R")
a <- commandArgs(TRUE)
scenario <- a[1L]
k <- as.integer(a[2L])
B <- as.integer(a[3L])
R <- as.integer(a[4L])
mode <- match.arg(a[5L], c("baseline","free","ridge","matched"))
assessment <- if (length(a) >= 6L) match.arg(a[6L], c("surface","observed")) else "surface"
n <- 70L
M <- 850L
ridge <- if(mode=="ridge") .01 else .001
monotone <- mode!="free"
workers <- as.integer(Sys.getenv("MS_GDUM_WORKERS","1"))
id <- scenario_id(scenario)
target <- read_target(scenario,k)
effect <- function(train,idx) {
  d <- ugdm_component_design(train$e,idx,"additive")
  y <- ugdm_component_response(train$x,idx)
  z <- ugdm_component_design(target$e,target$idx,reference=d)
  f <- ugdm_fit_components(d,y,"coherent",ridge=ridge,monotone_total=monotone)
  g <- importance(f,d,y,z,if (assessment=="observed") target$observed else target$expected)
  p <- g$point
  first <- p[p$predictor=="phosphorus",]
  second <- p[p$predictor=="secondary",]
  i <- match(responses,first$response)
  j <- match(responses,second$response)
  c(first$absolute[i],second$absolute[j],first$absolute[i]-second$absolute[j],
    active=sum(f$fits$total$coef[which(f$keep %in% f$direction)]<=1e-8), failed=g$failures)
}
labels <- paste(rep(c("phosphorus","secondary","contrast"),each=3L),responses,sep="_")
large <- lapply(1:2,function(j) {
  t <- generate(2100L,scenario,6817L+10L*id+j)
  effect(t,sampling(nrow(t$e),k,16000L,9029L+j+19L*k))[seq_along(labels)]
})
truth <- colMeans(do.call(rbind,large))
# The population game has no predictive gain under the known independence
# mechanism. Large finite communities still have chance alignments.
if (scenario=="null" && assessment=="observed") truth[] <- 0
one <- function(i) {
  seed <- 61883L+id*19831L+i*39L
  g <- generate(n,scenario,seed)
  point <- effect(g,sampling(n,k,M,seed+312L+k))
  pair <- effect(g,sampling(n,k,M,seed+2187L+k))
  finite <- effect(g,sampling(n,k,if(k==2L) choose(n,2L) else 16000L,seed+7353L+k))
  draws <- t(vapply(seq_len(R),function(j) {
    set.seed(seed+1017L+k+103L*j)
    idx <- sample_sets_weighted(site_weights(n,1L),k,M)
    effect(g,idx)[seq_along(labels)]
  },numeric(length(labels))))
  rows <- lapply(seq_along(labels),function(j) {
    x <- draws[is.finite(draws[,j]),j]
    valid <- length(x)>=.8*R && is.finite(point[j])
    q <- if(valid) quantile(x,c(.025,.10,.90,.975),names=FALSE) else rep(NA_real_,4)
    data.frame(scenario,order=k,n,replicate=i,seed,mode,M,draws=R,estimand=labels[j],assessment,
      truth=truth[j],point=point[j],repeat_point=pair[j],finite_target=finite[j],
      bootstrap_mean=if(valid) mean(x) else NA_real_,
      corrected=if(valid) 2*point[j]-mean(x) else NA_real_,
      lo95=q[1L],hi95=q[4L],lo80=q[2L],hi80=q[3L],
      basic_lo95=2*point[j]-q[4L],basic_hi95=2*point[j]-q[1L],
      basic_lo80=2*point[j]-q[3L],basic_hi80=2*point[j]-q[2L],
      active=point[["active"]],failed=point[["failed"]],valid_draws=length(x))
  })
  if (i %% 10L == 0L) message(scenario," ",k," ",mode,": ",i,"/",B)
  do.call(rbind,rows)
}
run <- if (workers>1L && .Platform$OS.type!="windows") {
  parallel::mclapply(seq_len(B),one,mc.cores=workers,mc.set.seed=FALSE)
} else lapply(seq_len(B),one)
stopifnot(all(vapply(run,is.data.frame,logical(1))))
atomic_csv(do.call(rbind,run),paste0("interval_",scenario,"_",k,"_",mode,if(assessment=="observed") "_observed", ".csv.gz"))
