source("docs/phase5_refinement/common.R")
Rcpp::sourceCpp("docs/manuscript_validation/calibration_helpers.cpp",
  cacheDir = Sys.getenv("CALIBRATION_CPP_CACHE", tempdir()))
args <- commandArgs(TRUE)
n <- as.integer(args[1])
signal <- args[2]
nr <- if (length(args)>2L) as.integer(args[3]) else 24L
grid <- c(900L,5000L,20000L,50000L)
out <- cover <- refs <- list()
point <- function(x, keep=rep(TRUE,length(x$y)), init=numeric(3))
  calibration_fit(x$X,x$y,unname(init),keep,rep(1,length(x$y)))
parameter <- function(pool) {
  a <- unlist(lapply(pool,point))
  deleted <- do.call(cbind,lapply(seq_along(pool),function(j) {
    x <- pool[[j]]
    idx <- x$design$idx
    sites <- split(rep(seq_len(nrow(idx)),ncol(idx)),as.vector(idx))
    t(vapply(seq_len(n),function(i) {
      keep <- rep(TRUE,nrow(idx))
      keep[sites[[as.character(i)]]] <- FALSE
      point(x,keep,a[(j-1L)*3L+1:3])
    },numeric(3)))
  }))
  list(point=a, pseudo = sweep(-(n-1)*deleted,2L,n*a,"+"))
}
score <- function(pool1,pool2,theta) {
  cache <- lapply(seq_along(pool1),function(j)
    internal(".prepare_score_pools")(list(pool1[[j]],pool2[[j]])))
  lapply(c("deletion","centered","inclusion","conditional"),function(m) {
    value <- function(which) {
      do.call(cbind,lapply(seq_along(cache),function(j) {
        p <- cache[[j]]
        b <- theta[(j-1L)*3L+1:3]
        mid <- stats::plogis(drop(p[[1L]]$X%*%b))
        hi <- stats::plogis(drop(p[[2L]]$X%*%b))
        H <- (crossprod(p[[1L]]$X*sqrt(mid*(1-mid))) / nrow(p[[1L]]$X) +
          crossprod(p[[2L]]$X*sqrt(hi*(1-hi))) / nrow(p[[2L]]$X)) / 2
        v <- internal(".pool_scores")(p[which],b,m)
        centers(v)%*%solve(H)
      }))
    }
    list(first=value(1),second=value(2),pooled=value(1:2))
  })
}
intervals <- function(point, V, weights, truth, method, id) {
  v <- centers(V)
  est <- drop(map6%*%point)
  J <- sqrt(diag(crossprod(v)/(n*(n-1))))
  dr <- sweep(crossprod(weights,v),2L,est,"+")
  result <- lapply(c(.5,.8,.9,.95),function(level) {
    if (grepl("^jack",method)) {
      lo <- est-qt((1+level)/2,n-1)*J
      hi <- est+qt((1+level)/2,n-1)*J
    } else {
      lim <- apply(dr,2L,quantile,probs=c((1-level)/2,(1+level)/2),names=FALSE)
      lo <- lim[1,]
      hi <- lim[2,]
    }
    data.frame(id,method,effect=rownames(map6),level,truth=unname(truth),
      estimate=unname(est),width=unname(hi-lo),covered=unname(lo<=truth & hi>=truth))
  })
  do.call(rbind,result)
}
for (r in seq_len(nr)) {
  ec <- j_layout(n,r,signal)
  tt <- drop(map6%*%ec$truth)
  W <- site_weights(n,599L)
  W <- sweep(W,2L,colSums(W),"/")
  full <- if (choose(n,max(ec$orders))<=60000) parameter(ec$pool(0,exact=TRUE))
  for (s in c("uniform","site")) for (M in grid) {
    t0 <- sum(proc.time()[1:2])
    a <- ec$pool(M,s)
    b <- ec$pool(M,s)
    f <- parameter(a)
    g <- parameter(b)
    pw1 <- f$pseudo %*% t(map6)
    pw2 <- g$pseudo %*% t(map6)
    id <- data.frame(n,signal,r,order1=ec$orders[1],order2=ec$orders[2],sampler=s,M)
    summarise_noise <- function(x,y,method,elapsed) {
      q <- ugdm_pseudovalue_noise(x,y)
      z <- q$numerical_fraction
      data.frame(id,method,effect=rownames(map6),mc_fraction=unname(z),
        MC_var=unname(diag(q$numerical_covariance)),JK_var=unname(diag(q$average_jackknife)),
        pooled_var=unname(diag(q$pooled_jackknife)),cpu=elapsed)
    }
    out[[length(out)+1L]] <- summarise_noise(pw1,pw2,"parameter",sum(proc.time()[1:2])-t0)
    avg <- (f$point+g$point)/2
    cover[[length(cover)+1L]] <- intervals(f$point,pw1,W,tt,"jackknife",id)
    cover[[length(cover)+1L]] <- intervals(f$point,pw1,W,tt,"parameter_BB",id)
    cover[[length(cover)+1L]] <- intervals(avg,(pw1+pw2)/2,W,tt,"parameter_BB_two",id)
    zz <- score(a,b,avg)
    names(zz) <- c("deletion","centered","inclusion","conditional")
    for (m in names(zz)) {
      l <- zz[[m]]
      out[[length(out)+1L]] <- summarise_noise(l$first%*%t(map6),l$second%*%t(map6),m,NA_real_)
      cover[[length(cover)+1L]] <- intervals(avg,l$pooled%*%t(map6),W,tt,paste0("score_",m),id)
    }
    if (!is.null(full)) {
      E <- centers(full$pseudo%*%t(map6))
      refs[[length(refs)+1L]] <- data.frame(id,effect=rownames(map6),
        exact=drop(map6%*%full$point),exact_jack=unname(internal(".jackknife_sd")(-E/(n-1))),
        point_error=drop(map6%*%(f$point-full$point)),
        pseudo_RMSE=sqrt(colMeans((centers(pw1)-E)^2)),
        pooled_pseudo_RMSE=sqrt(colMeans((centers((pw1+pw2)/2)-E)^2)))
    }
  }
  if (!is.null(full)) {
    id <- data.frame(n,signal,r,order1=ec$orders[1],order2=ec$orders[2],sampler="complete",M=0)
    cover[[length(cover)+1L]] <- intervals(full$point,full$pseudo%*%t(map6),W,tt,"jackknife",id)
    cover[[length(cover)+1L]] <- intervals(full$point,full$pseudo%*%t(map6),W,tt,"parameter_BB",id)
  }
}
key <- paste(n,signal,sep="_")
write_refinement(do.call(rbind,out),paste0("noise_",key))
write_refinement(do.call(rbind,cover),paste0("coverage_",key))
if (length(refs)) write_refinement(do.call(rbind,refs),paste0("exact_",key))
