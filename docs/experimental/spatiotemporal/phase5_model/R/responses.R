.st5_binary_metrics <- c("zeta_presence","incidence","turnover","nestedness",
                         "temporal_presence","site_temporal_presence")
.st5_abundance_metrics <- c("zeta_presence","zeta_abundance","zeta_dominance","incidence",
  "turnover","nestedness","abundance","occupancy","numerical","temporal_presence",
  "temporal_abundance","site_temporal_presence","site_temporal_abundance",
  "site_temporal_dominance","temporal_dominance","dominance")
.st5_abundance_only <- setdiff(.st5_abundance_metrics,.st5_binary_metrics)

st5_responses <- function(p, support, processes = FALSE, draws = 8L,
                          exact_limit = 32L, seed = 20L) {
  type <- attr(support,"type")
  d <- if (is.null(p$a)) 3L else 4L
  if (processes && type %in% c("spatial","rectangle")) stop("Supply paired support for attribution.")
  metrics <- if (d == 3L) .st5_binary_metrics else .st5_abundance_metrics
  metrics <- switch(type,
    local = intersect(c("site_temporal_presence","site_temporal_abundance","site_temporal_dominance"),metrics),
    shared = intersect(c("temporal_presence","temporal_abundance","temporal_dominance"),metrics),
    setdiff(metrics,grep("temporal",metrics,value=TRUE)))
  blocks <- lapply(seq_along(support),function(h) {
    g <- support[[h]]
    a <- st5_observation_ids(p,g$sets,g$a)
    b <- st5_observation_ids(p,g$sets,g$b)
    m <- length(g$a)
    n <- length(g$region)
    k <- ncol(g$sets)
    lag <- if (length(g$b)) mean(p$times[g$b]-p$times[g$a]) else NA_real_
    span <- diff(range(p$times[g$a]))
    span1 <- if (length(g$b)) diff(range(p$times[g$b])) else NA_real_
    flatten <- function(z) as.vector(p$index[g$region,z,drop=FALSE])
    row0 <- flatten(g$a)
    row1 <- if (length(g$b)) flatten(g$b) else row0
    sets <- matrix(match(g$sets,g$region),ncol=k)
    region_a <- if (d==3L) p$x else p$a
    A <- region_a[row0,,drop=FALSE]
    B <- region_a[row1,,drop=FALSE]
    active <- colSums((A>0) | (B>0)) > 0
    if (any(active)) { A <- A[,active,drop=FALSE]; B <- B[,active,drop=FALSE] }
    if (d==4L) {
      # Coalition edits use closed input shares. Empty rows stay zero until evaluation.
      close <- function(z) z / pmax(rowSums(z),.Machine$double.xmin)
      A <- close(A); B <- close(B)
    }
    alloc <- list(masks=list(matrix(1L,nrow(A),ncol(A))),exact=TRUE)
    if (processes) {
      panel <- st_panel(array(rbind(A,B),c(n,2L*m,ncol(A)),
                 dimnames=list(p$sites[g$region],NULL,colnames(A))), times=p$times[c(g$a,g$b)])
      # rbind() places all A observations before all B observations within each taxon.
      plan <- st3_plan(panel,seq_len(m),m+seq_len(m))
      alloc <- st3_allocations(plan,exact_limit,draws,seed+h)
    }
    if (d==3L) {
      storage.mode(A) <- storage.mode(B) <- "integer"
      r <- annual_incidence(A,B,sets,n,if (processes) alloc$masks else list())
      raw_metrics <- .st5_binary_metrics
      phi <- r$draws
      error <- r$error
    } else {
      idx <- .st3_site_indices(sets,n,m)
      indices <- lapply(seq_len(nrow(idx)),function(i) as.integer(idx[i,]))
      r <- st4_evaluate_cpp(A,B,indices,alloc$masks,
                           if(processes) .st3_coalitions() else matrix(0,16,0))
      raw_metrics <- .st5_abundance_metrics
      phi <- r$phi
      error <- r$efficiency_error
      if (processes) {
        aa <- A>0; bb <- B>0
        storage.mode(aa) <- storage.mode(bb) <- "integer"
        binary <- annual_incidence(aa,bb,sets,n,alloc$masks)
        keep <- match(.st5_binary_metrics,raw_metrics)
        if (max(c(0,abs(phi[,keep,1:3,,drop=FALSE]-binary$draws)),na.rm=TRUE)>1e-10)
          stop("Three-block and four-block incidence results disagree.")
        phi[,keep,1:3,] <- binary$draws
        error <- max(error,binary$error)
      }
    }
    jj <- match(metrics,raw_metrics)
    base <- r$before[,jj,drop=FALSE]
    final <- r$after[,jj,drop=FALSE]
    colnames(base) <- colnames(final) <- metrics
    normalized <- !startsWith(metrics,"zeta_")
    clip <- function(x) {
      if (any(x< -1e-10 | x>1+1e-10,na.rm=TRUE)) stop("Invalid normalized observation.")
      pmin(pmax(x,0),1)
    }
    base[,normalized] <- clip(base[,normalized])
    final[,normalized] <- clip(final[,normalized])
    gain <- loss <- shared <- rep(NA_real_,nrow(sets))
    if (type=="local") {
      x0 <- p$x[a[,1],,drop=FALSE]
      x1 <- p$x[b[,1],,drop=FALSE]
      gain <- rowSums(!x0 & x1); loss <- rowSums(x0 & !x1); shared <- rowSums(x0 & x1)
    }
    z <- data.frame(period=h,k=k,m=if(type=="local") 2L else m,
       initial=min(p$times[g$a]),final=if(length(g$b)) max(p$times[g$b]) else max(p$times[g$a]),
       lag=lag,span=span,final_span=span1,gain,loss,shared,
       gain_component=.st_ratio(gain,2*shared+gain+loss),
       loss_component=.st_ratio(loss,2*shared+gain+loss), regional_sites=n,
       eligible=g$eligible, allocations=length(alloc$masks), exact_allocation=alloc$exact)
    z$sites <- I(lapply(seq_len(nrow(sets)),function(i) p$sites[g$sets[i,]]))
    z$dates <- I(rep(list(p$times[c(g$a,g$b)]),nrow(sets)))
    z$a <- I(lapply(seq_len(nrow(a)),function(i) a[i,]))
    z$b <- I(lapply(seq_len(nrow(b)),function(i) b[i,]))
    contributions <- allocation_se <- NULL
    if (processes) {
      if (error>1e-10) stop("Shapley additivity failed.")
      phi <- phi[,jj,,,drop=FALSE]
      contributions <- apply(phi,c(1L,2L,3L),mean)
      allocation_se <- apply(phi,c(1L,2L,3L),function(x) {
        if (anyNA(x)) NA_real_ else if(alloc$exact) 0 else sd(x)/sqrt(length(x))
      })
      dimnames(contributions) <- dimnames(allocation_se) <- list(NULL,metrics,st3_processes[seq_len(d)])
    }
    list(rows=z,before=base,after=final,phi=contributions,allocation_se=allocation_se,
         error=if(processes) error else 0)
  })
  bind <- function(name) do.call(rbind,lapply(blocks,`[[`,name))
  rows <- bind("rows")
  reassemble <- function(name) {
    a <- lapply(blocks,`[[`,name)
    if (!length(a) || is.null(a[[1]])) return(NULL)
    out <- array(NA_real_,c(nrow(rows),length(metrics),d),dimnames=dimnames(a[[1]]))
    for (i in seq_along(a)) out[rows$period==i,,] <- a[[i]]
    out
  }
  if (!length(blocks)) return(list(rows=data.frame(),before=matrix(NA_real_,0,length(metrics),dimnames=list(NULL,metrics)),
     after=matrix(NA_real_,0,length(metrics),dimnames=list(NULL,metrics)),type=type,
     phi=NULL,allocation_se=NULL,efficiency_error=NA_real_))
  list(rows=rows,before=bind("before"),after=bind("after"),type=type,
       phi=reassemble("phi"),allocation_se=reassemble("allocation_se"),
       efficiency_error=max(vapply(blocks,`[[`,numeric(1),"error")))
}

st5_extract <- function(x, metric, processes = !is.null(x$phi)) {
  if (!metric %in% colnames(x$before)) stop("Response unavailable.")
  baseline <- if (x$type=="change") x$before[,metric] else rep(0,nrow(x$rows))
  endpoint <- if (x$type %in% c("change","shared","local")) x$after[,metric] else x$before[,metric]
  phi <- if (processes) {
    keys <- dimnames(x$phi)[[3]]
    if(metric %in% .st5_binary_metrics) keys <- intersect(keys,st3_processes[1:3])
    matrix(x$phi[,metric,keys],nrow(x$rows),dimnames=list(NULL,keys))
  }
  list(rows=x$rows, baseline=baseline, y=endpoint, phi=phi,type=x$type,metric=metric,
       weight_group=if(x$type=="local") vapply(x$rows$sites,`[[`,character(1),1L) else x$rows$period)
}

st5_subset <- function(z,i) {
  for (name in c("rows","baseline","y","phi","weight_group","composition","crossed","importance_weight"))
    if (!is.null(z[[name]])) z[[name]] <- if (length(dim(z[[name]]))) {
      if(length(dim(z[[name]])) == 3L) z[[name]][i,,,drop=FALSE] else z[[name]][i,,drop=FALSE]
    } else z[[name]][i]
  z
}

st5_crossed_data <- function(x) {
  z <- st5_extract(x,"abundance")
  z$composition <- cbind(before_occ=x$before[,"occupancy"],before_num=x$before[,"numerical"],
                         after_occ=x$after[,"occupancy"],after_num=x$after[,"numerical"])
  z$crossed <- aperm(x$phi[,c("occupancy","numerical"),,drop=FALSE],c(1,3,2))
  z
}

st5_process_summary <- function(phi,w=rep(1,nrow(phi))) {
  good <- complete.cases(phi) & is.finite(w) & w>0
  mean_w <- function(a) if (any(good)) colSums(a[good,,drop=FALSE]*w[good])/sum(w[good]) else rep(NA,ncol(a))
  negative <- mean_w(pmax(-phi,0))
  positive <- mean_w(pmax(phi,0))
  data.frame(process=colnames(phi),negative,positive,net=positive-negative,gross=positive+negative,
             negative_share=100*.st_ratio(negative,sum(negative)),
             positive_share=100*.st_ratio(positive,sum(positive)))
}

st5_support_change <- function(spatial,change,metric,times) {
  sample_means <- tapply(spatial$before[,metric],spatial$rows$initial,mean,na.rm=TRUE)
  do.call(rbind,lapply(seq_len(nrow(times)),function(j) {
    a <- times[j,1]; b <- times[j,2]
    all <- unname(sample_means[as.character(b)]-sample_means[as.character(a)])
    i <- which(change$rows$initial==a & change$rows$final==b)
    d <- change$after[i,metric]-change$before[i,metric]
    if (!is.null(change$phi)) d[rowSums(matrix(!is.finite(change$phi[i,metric,]),nrow=length(i),ncol=dim(change$phi)[3]))>0] <- NA_real_
    matched <- if (any(is.finite(d))) mean(d,na.rm=TRUE) else NA_real_
    data.frame(initial=a,final=b,all=all,matched=matched,discrepancy=all-matched,
               matched_valid=sum(is.finite(d)))
  }))
}
