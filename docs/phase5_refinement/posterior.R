source("docs/phase5_refinement/common.R")
args <- commandArgs(TRUE)
mode <- match.arg(args[1],c("calibration","reference","mixture","long_reference"))
signal <- args[2]
start <- if (length(args)>2) as.integer(args[3]) else 1L
end <- if (length(args)>3) as.integer(args[4]) else if (mode=="calibration") 16L else 4L
grid <- if (length(args)>4L) as.integer(strsplit(args[5],",",fixed=TRUE)[[1L]]) else
  if (mode!="calibration") c(900L,5000L) else c(900L,5000L,20000L,50000L)
calc <- diagnostics <- hull <- list()
for (r in start:end) for (M in grid) {
  n <- 70L
  q <- j_layout(n,r,signal)
  a <- q$pool(M,"uniform")
  b <- q$pool(M,"uniform")
  init <- unlist(lapply(a,function(z) suppressWarnings(glm.fit(z$X,z$y,family=binomial()))$coefficients),
                 use.names=FALSE)
  id <- data.frame(mode,signal,r,M,order1=q$orders[1],order2=q$orders[2])
  message(mode, " ", signal, " draw ", r, " M=", M)
  guarded <- function(fn) tryCatch(fn(),error=function(e) list(error=conditionMessage(e)))
  run <- list()
  if (mode=="mixture") {
    run$chain_mixture <- guarded(function() ugdm_jel_refined(a,draws=1400L,pilot=400L,
      pilot_kernel="chains",certify=TRUE,seed=521+r))
  } else if (mode!="long_reference") {
    run$adaptive <- guarded(function() ugdm_jel_refined(a,draws=if (mode=="reference") 1400L else 550L,
      pilot=200L, certify=TRUE,seed=7103+r))
  }
  if (mode=="reference") {
    run$legacy <- guarded(function() ugdm_jel(a,draws=1200L,seed=8931+r))
  }
  if (mode %in% c("reference","long_reference")) {
    long <- mode=="long_reference"
    run$mcmc <- guarded(function() ugdm_jel_refined(a,integration="mcmc",
      draws=if (long) 32000L else 3200L, warmup=if (long) 5500L else 1500L,
      adapt_covariance=TRUE, seed=if (long) 19421+r else 4128+r))
  } else if (mode=="calibration" && M %in% c(900L,5000L,50000L)) {
    pools <- Map(function(x,y) list(x,y),a,b)
    run$conditional_two <- guarded(function() ugdm_jel_refined(pools,score="conditional",
      draws=550L,pilot=200L,certify=TRUE,seed=8175+r))
  }
  if (mode=="reference" && M==900L) {
    candidate <- init+c(2,0,0,2,0,0)
    z <- do.call(cbind,lapply(1:2,function(j)
      ugdm_pool_scores(list(a[[j]]),candidate[(j-1L)*3L+1:3],"centered")))
    plain <- ugdm_empirical_likelihood(z,certify=TRUE)
    if (identical(plain$reason,"zero strictly outside the moment hull")) {
      adj <- ugdm_empirical_likelihood(z,adjusted=TRUE)
      hull[[length(hull)+1L]] <- data.frame(id,ordinary=plain$reason,adjusted_ok=adj$ok,
        adjusted_loglik=adj$loglik,artificial_weight=adj$weights[n+1L])
    }
  }
  posterior <- function(x, weights) c(mean=sum(x*weights), sd=sqrt(sum(weights*(x-sum(x*weights))^2)),
    p=sum(weights[x>0]), practical=sum(weights[x>.10]),
    qweight(x,weights,c(.025,.05,.10,.25,.75,.90,.95,.975)))
  for (method in names(run)) {
    out <- run[[method]]
    if (!is.null(out$error)) {
      diagnostics[[length(diagnostics)+1L]] <- data.frame(id,method,ess=NA,rhat=NA,
        hull=NA,unresolved=NA,rank=NA,p_gap=NA,stable=FALSE,seconds=NA,
        mode_convergence=NA,clipped=NA,accept=NA,error=out$error)
      next
    }
    theta <- out$draws%*%t(map6)
    all <- t(vapply(seq_len(ncol(theta)),function(j) posterior(theta[,j],out$weights),numeric(12)))
    p_gap <- rep(NA_real_,nrow(all))
    mean_gap <- rep(NA_real_,nrow(all))
    checked <- rep(FALSE,nrow(all))
    ess <- if (is.null(out$ess)) NA_real_ else min(out$ess)
    rhat <- if (is.null(out$rhat)) NA_real_ else max(out$rhat)
    count <- function(k) if (is.null(out$likelihood_status)) NA_real_ else sum(out$likelihood_status[names(out$likelihood_status)==k])
    if (method == "mcmc") {
      V <- internal(".jel_chain_diagnostics")(split(seq_len(nrow(theta)),out$group) |>
        lapply(function(i) theta[i,,drop=FALSE]))
      checked <- V$rhat<=1.05 & V$ess>=200
      ess <- min(c(out$ess,V$ess))
      rhat <- max(c(out$rhat,V$rhat))
      checked <- checked & ess>=200 & rhat<=1.05
    } else if (!is.null(out$group)) {
      chunks <- split(seq_len(nrow(theta)),out$group)
      tmp <- lapply(chunks,function(i) t(vapply(seq_len(ncol(theta)),function(j)
        posterior(theta[i,j],out$weights[i]/sum(out$weights[i])),numeric(12))))
      p_gap <- apply(abs(tmp[[1]][,3:4]-tmp[[2]][,3:4]),1,max)
      mean_gap <- abs(tmp[[1]][,1]-tmp[[2]][,1])/all[,2]
      # Keep per-contrast agreement after applying common posterior diagnostics.
      checked <- p_gap<=.10 & mean_gap<=.20 & ess>=120 & min(out$batch_ess)>=60 &
        out$mode_convergence==0 & out$clipped==0
    }
    diagnostics[[length(diagnostics)+1L]] <- data.frame(id,method,ess,rhat,
      hull=count("zero strictly outside the moment hull"),unresolved=count("no certified interior solution"),
      rank=count("rank-deficient moments"),p_gap=if (any(is.finite(p_gap))) max(p_gap) else NA,
      stable=all(checked),seconds=out$seconds,
      mode_convergence=if (is.null(out$mode_convergence)) NA else out$mode_convergence,
      clipped=if (is.null(out$clipped)) NA else out$clipped,
      accept=if (is.null(out$acceptance)) NA else mean(out$acceptance),error="")
    calc[[length(calc)+1L]] <- data.frame(id,method,effect=rownames(map6),
      truth=drop(map6%*%q$truth),stable=checked,p_gap,mean_gap,
      mean=all[,1],sd=all[,2],p=all[,3],practical=all[,4],
      l95=all[,5],l90=all[,6],l80=all[,7],l50=all[,8],h50=all[,9],
      h80=all[,10],h90=all[,11],h95=all[,12])
  }
}
key <- paste(mode,signal,start,end,paste(grid,collapse="-"),sep="_")
write_refinement(do.call(rbind,calc),paste0("posterior_",key))
write_refinement(do.call(rbind,diagnostics),paste0("postdiag_",key))
if(length(hull)) write_refinement(do.call(rbind,hull),paste0("hull_",key))
