source("docs/phase5/common.R")
args <- commandArgs(TRUE)
meta <- read_manifest()[as.integer(args[2]),]
d <- load_dd2(meta,args[1])
x <- d$x
env <- d$env
geo <- d$dist
n <- nrow(x)
orders <- unlist(meta[c("low","mid","high")],use.names=FALSE)
ref <- suppressWarnings(ugdm_design(env,cbind(seq_len(n-1L),2:n),geo,dist_range=range(as.dist(geo))))
focus <- c("spread","position","relative_position")
eco_limit <- c(.10,.20,.15)
abs_tau <- setNames(c(.012,.025,.055,.012,.025,.055),c(focus,paste0(focus,"_change")))
build <- function(M, sampler, k = orders[1:2]) {
  ss <- sample_sets_coupled(x,k,M,coupled=sampler!="independent", exact_limit=50000)
  ugdm_path_data(x,env,ss,geo,ref)
}
set.seed(412 + n)
a <- build(3500,"coupled",orders)
b <- build(3500,"coupled")
rows <- pilot <- budgets <- list()
for (component in c("incidence","relative")) {
  rr <- function(path, deletion) {
    lapply(path, function(z) {
      b <- model_data(z,component,"B")
      f <- fit5(b,"B")
      get <- function(a, fit) setNames(effect5(a,fit,env),effect_names)[focus]
      fit_summary <- setNames(effect5(b,f,env), effect_names)
      point <- fit_summary[focus]
      de <- if (deletion) t(vapply(seq_len(n), function(i) {
        b <- model_data(z,component,"B",keep=rowSums(z$design$idx==i)==0L)
        if (!length(b$y) || !f$converged) return(rep(NA_real_,length(focus)))
        get(b,fit5(b,"B",start=f))
      }, numeric(length(focus))))
      list(point=point,deletion=de, condition=unname(fit_summary[8]))
    })
  }
  d1 <- rr(a[1:2],TRUE)
  d2 <- rr(b,TRUE)
  reliable <- function(v,w) {
    se <- pmin(internal(".jackknife_sd")(v),internal(".jackknife_sd")(w))
    noise <- internal(".jackknife_sd")(v-w)/sqrt(2)
    list(se=se, noise=noise/se, stable=is.finite(se) & se>0 & is.finite(noise) & noise<=.50*se)
  }
  zz <- list(reliable(d1[[1]]$deletion,d2[[1]]$deletion),
    reliable(d1[[2]]$deletion,d2[[2]]$deletion),
    reliable(d1[[2]]$deletion-d1[[1]]$deletion,d2[[2]]$deletion-d2[[1]]$deletion))
  ec <- setNames(c(zz[[2]]$se,zz[[3]]$se),names(abs_tau))
  ok <- setNames(c(zz[[2]]$stable,zz[[3]]$stable),names(abs_tau))
  fit_once <- function(M, sampler) {
    z <- rr(build(M,sampler),FALSE)
    c(z[[2]]$point,setNames(z[[2]]$point-z[[1]]$point,paste0(focus,"_change")))
  }
  buy <- ugdm_contrast_budget(fit_once,c("coupled","independent"),c(600L,1500L,3500L),abs_tau,
    ecological_se=ec, calibrated=ok, fraction=.2, nrep=6,max_budget=15000L,seed=177+n)
  if (!is.null(buy$contrasts)) budgets[[length(budgets)+1L]] <- cbind(component,sampler=buy$sampler,buy$contrasts)
  pilot[[length(pilot)+1L]] <- cbind(component,buy$pilots)
  m <- if (component=="relative") .002 else .005
  condition <- lapply(a,function(s) sd(s[[component]],na.rm=TRUE))
  vals <- lapply(1:4,function(j) {
    pos <- if (j==4L) 3L else min(j,2L)
    earlier <- if (j==3L) orders[1] else NA_integer_
    terms <- if (j==3L) paste0(focus,"_change") else focus
    ec <- if (j==4L) rep(NA,length(focus)) else zz[[j]]$stable & zz[[j]]$se<=eco_limit
    if (j==1L) numerical <- rep(a[[1L]]$exact,length(focus)) else if (!is.null(buy$contrasts) && j<4L) {
      numerical <- buy$contrasts$passed[match(terms,buy$contrasts$effect)]
    } else numerical <- rep(NA,length(focus))
    st <- ugdm_order_status(rep(orders[pos],length(focus)),terms,
      rep(if (j==3L) min(condition[[1]],condition[[2]])>=m else condition[[pos]]>=m,length(focus)),ec,numerical)
    se <- if (j==4L) rep(NA_real_,length(focus)) else unname(zz[[j]]$se)
    noise <- if (j==4L) rep(NA_real_,length(focus)) else unname(zz[[j]]$noise)
    stable <- if (j==4L) rep(NA,length(focus)) else unname(zz[[j]]$stable)
    cond <- if (j==4L) NA_real_ else if (j==3L) max(d1[[1L]]$condition,d1[[2L]]$condition) else d1[[j]]$condition
    cbind(component,previous=earlier,response_sd=condition[[pos]],response_sd_threshold=m,
          ecological_se=se,ecological_se_limit=eco_limit,deletion_noise_ratio=noise,
          stable_deletions=stable,scaled_information_condition=cond,st$profile)
  })
  rows[[length(rows)+1L]] <- do.call(rbind,vals)
}
stamp <- function(z) cbind(community=meta$community,n=n,z)
write5(stamp(do.call(rbind,rows)),paste0("order_status_",meta$id))
if (length(budgets)) write5(stamp(do.call(rbind,budgets)),paste0("budget_",meta$id))
write5(stamp(do.call(rbind,pilot)),paste0("budget_pilot_",meta$id))
