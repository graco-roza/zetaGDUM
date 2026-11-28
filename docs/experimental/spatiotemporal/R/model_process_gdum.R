st3_basis <- function(x, env, training_times = NULL) {
  p <- .st_panel(x)
  e <- .st_environment(p, env)
  if (is.null(training_times)) training_times <- seq_along(p$times)
  train <- .st_orders(training_times, length(p$times), "training_times")
  flat <- .st_flat_env(e[, train, , drop = FALSE])
  d <- zetaGDUM::ugdm_design(flat, matrix(c(1L,2L), 1L))
  list(native = d$environment_basis, range = apply(flat, 2L, range))
}

st3_environment_design <- function(x, env, sets, windows, basis = NULL, limit = 1e5) {
  p <- .st_panel(x)
  e <- .st_environment(p, env)
  n <- dim(e)[1L]
  if (!is.matrix(sets) || ncol(sets) < 1L) stop("Supply one matrix of fixed-order spatial combinations.")
  k <- ncol(sets)
  sets <- .st_sets(n, k, sets, limit)
  if (!is.list(windows) || !length(windows) || length(windows) * as.double(nrow(sets)) > limit) {
    stop("Supply paired windows within the design row budget.")
  }
  m <- length(windows[[1L]]$before)
  before <- after <- matrix(NA_integer_, length(windows), m)
  for (i in seq_along(windows)) {
    w <- windows[[i]]
    a <- .st_orders(w$before, length(p$times), "before")
    b <- .st_orders(w$after, length(p$times), "after")
    if (length(a) != m || length(b) != m || any(diff(a) <= 0) || any(diff(b) <= 0) || max(a) >= min(b)) {
      stop("Compare chronological non-overlapping windows of the same temporal order.")
    }
    before[i, ] <- a
    after[i, ] <- b
  }
  if (is.null(basis)) basis <- st3_basis(p, e, sort(unique(c(before,after))))
  vars <- dimnames(e)[[3L]]
  if (!identical(vars, names(basis$native$center))) stop("Align predictors to the training basis.")
  means <- function(w) {
    z <- lapply(seq_len(nrow(w)), function(j) {
      matrix(colMeans(aperm(e[,w[j,],,drop=FALSE],c(2,1,3)), dims=1), n, length(vars))
    })
    out <- do.call(rbind,z)
    colnames(out) <- vars
    out
  }
  flat <- rbind(means(before),means(after))
  rows <- expand.grid(set=seq_len(nrow(sets)),contrast=seq_len(nrow(before)))
  idx <- sets[rows$set,,drop=FALSE]
  a <- idx + n*(rows$contrast-1L)
  b <- a + n*nrow(before)
  lag <- rowMeans(matrix(p$times[after]-p$times[before],nrow(before)))
  start <- rowMeans(matrix(p$times[before],nrow(before)))
  end <- rowMeans(matrix(p$times[after],nrow(after)))
  rename <- function(z,label) {
    colnames(z) <- paste0(label,".",colnames(z))
    z
  }
  if (k == 1L) {
    Z <- scale(flat,basis$native$center,basis$native$scale)
    P <- Z[a,,drop=FALSE]
    delta <- Z[b,,drop=FALSE]-P
    spread <- NULL
  } else {
    d <- zetaGDUM::ugdm_design(flat,rbind(a,b),ref=list(environment_basis=basis$native))
    i <- seq_len(nrow(a))
    j <- i+nrow(a)
    P <- d$Xpos[i,,drop=FALSE]
    delta <- d$Xpos[j,,drop=FALSE]-P
    S <- d$Xdir[i,,drop=FALSE]
    spread <- cbind(rename(S,"spread"),rename(d$Xdir[j,,drop=FALSE]-S,"spread_change"))
  }
  X <- cbind(rename(P,"position"),rename(delta,"position_change"),spread,
             log_lag=log(lag[rows$contrast]),rename(P*delta,"position_by_change"))
  outside <- sweep(flat,2,basis$range[1L,],"<") | sweep(flat,2,basis$range[2L,],">")
  extreme <- which(rowSums(outside)>0L)
  structure(list(X=X,sets=idx,rows=rows,basis=basis,k=k,m=m,
                 sites=dimnames(p$x)[[1L]],start=start,end=end,lag=lag,
                 environmental_extrapolation=rowSums(matrix(a %in% extreme | b %in% extreme,nrow(a)))>0),
            class="st3_design")
}

.st3_ridge <- function(X,Y,w,penalty) {
  weight <- sqrt(w)
  A <- crossprod(X*weight)
  G <- A+diag(penalty,ncol(X))
  coef <- solve(G,crossprod(X*weight,Y*weight))
  rownames(coef) <- colnames(X)
  colnames(coef) <- colnames(Y)
  list(coefficients=coef,gram=A,system=G)
}

st3_fit_process_gdum <- function(design, y, weights = rep(1,nrow(design$X)), ridge = 0.5,
                                 interaction = FALSE) {
  if (!inherits(design,"st3_design") || !is.matrix(y) || !is.numeric(y) ||
      nrow(y)!=nrow(design$X) || !ncol(y) || is.null(colnames(y)) ||
      anyDuplicated(colnames(y))) stop("Use a named numeric response matrix aligned to a Phase III design.")
  if (length(ridge)!=1L || !is.finite(ridge) || ridge<=0 ||
      length(weights)!=nrow(y) || any(!is.finite(weights) | weights<0)) stop("Invalid ridge or common weights.")
  keep <- rowSums(!is.finite(y))==0L & weights>0
  if (!any(keep)) stop("No common finite response support.")
  candidate <- colnames(design$X)
  if (!interaction) candidate <- candidate[!startsWith(candidate,"position_by_change.")]
  A <- design$X[keep,candidate,drop=FALSE]
  changing <- apply(A,2L,function(z) diff(range(z)) > 1e-12)
  if (interaction) {
    for (term in candidate[startsWith(candidate,"position_by_change.")]) {
      v <- sub("^position_by_change\\.","",term)
      if (!all(changing[paste0(c("position.","position_change."),v)])) changing[term] <- FALSE
    }
  }
  columns <- candidate[changing]
  X <- cbind(intercept=1,A[,changing,drop=FALSE])
  Y <- y[keep,,drop=FALSE]
  w <- weights[keep]
  penalty <- c(0,rep(ridge,ncol(X)-1L))
  fit <- .st3_ridge(X,Y,w,penalty)
  Z <- X*sqrt(w)
  decomposition <- qr(Z)
  profiled <- function(prefix) {
    j <- startsWith(colnames(Z),prefix)
    if (!any(j)) return(NA_real_)
    residual <- qr.resid(qr(Z[,!j,drop=FALSE]),Z[,j,drop=FALSE])
    sum(residual^2)/(sum(w)*sum(j))
  }
  structure(c(fit,list(design=design,columns=columns,dropped=candidate[!changing],keep=keep,
                       X=X,Y=Y,weights=w,penalty=penalty,ridge=ridge,
                       rank=decomposition$rank,singular_values=svd(Z,nu=0,nv=0)$d,
                       effective_df=sum(diag(solve(fit$system,fit$gram))),
                       effective_rows=sum(w)^2/sum(w^2),
                       profile_position=profiled("position."),profile_spread=profiled("spread."))),
            class="st3_fit")
}

st3_predict_process_gdum <- function(object, design = object$design, initial = NULL, tolerance = 1e-8) {
  if (!identical(design$basis,object$design$basis) || design$k!=object$design$k ||
      design$m!=object$design$m) stop("Prediction requires the training basis and matching orders.")
  X <- cbind(intercept=1,design$X[,object$columns,drop=FALSE])
  estimates <- X %*% object$coefficients
  limits <- apply(object$X,2L,range)
  extrapolation <- rowSums(sweep(X,2,limits[1L,],"<") | sweep(X,2,limits[2L,],">"))>0 |
    design$environmental_extrapolation
  total <- rowSums(estimates)
  bound <- if (is.null(initial)) rep(NA,nrow(X)) else {
    if (!is.numeric(initial) || length(initial)!=nrow(X) ||
        any(is.finite(initial) & (initial < -tolerance | initial > 1+tolerance))) {
      stop("Initial normalized dissimilarity must be in [0,1] with one value per row.")
    }
    initial+total < -tolerance | initial+total > 1+tolerance |
      rowSums(abs(estimates)>1+tolerance)>0
  }
  list(contributions=estimates,total=total,extrapolation=extrapolation,incompatible=bound)
}

# A common-support eight-column target for the crossed decomposition.
st3_model_response <- function(result, response = "abundance", crossed = FALSE) {
  if (!crossed) return(matrix(result$contributions[,response,], ncol=4L,
                               dimnames=list(NULL,st3_processes)))
  x <- st3_crossed(result)$contributions
  shape <- dim(x)
  out <- matrix(x,shape[1L],shape[2L]*shape[3L])
  colnames(out) <- as.vector(outer(dimnames(x)[[2L]],dimnames(x)[[3L]],paste,sep="."))
  out
}

.st3_moving <- function(n,width) {
  first <- sample.int(n-width+1L,ceiling(n/width),replace=TRUE)
  ids <- as.vector(outer(seq_len(width)-1L,first,"+"))
  tabulate(head(ids,n),n)
}

st3_contrast_uncertainty <- function(object, contrast, blocks, draws = 32L, width = 2L) {
  design <- object$design
  if (is.null(names(contrast)) || !identical(names(contrast),colnames(object$X)) ||
      length(blocks)!=length(design$sites) || anyNA(blocks) || length(unique(blocks))<3L) {
    stop("Supply an aligned design contrast and three or more named-site blocks.")
  }
  if (width<1L || width>length(design$start) || draws<2L) stop("Invalid temporal block width or draw count.")
  grid <- design$sets[object$keep,,drop=FALSE]
  membership <- matrix(match(blocks[grid],unique(blocks)),nrow(grid))
  ids <- seq_along(unique(blocks))
  chrono <- order(design$start,design$end)
  windows <- design$rows$contrast[object$keep]
  evaluate <- function(multiplier) {
    if (!any(multiplier>0)) return(rep(NA_real_,ncol(object$Y)))
    f <- .st3_ridge(object$X,object$Y,object$weights*multiplier,object$penalty)
    drop(contrast %*% f$coefficients)
  }
  deleted <- matrix(vapply(ids,function(j) evaluate(rowSums(membership==j)==0L),
                           numeric(ncol(object$Y))),nrow=ncol(object$Y))
  centre <- rowMeans(deleted)
  spatial <- sqrt((length(ids)-1L)*rowMeans((deleted-centre)^2))
  moving <- crossed <- matrix(NA_real_,ncol(object$Y),draws)
  unique_member <- lapply(seq_len(ncol(membership)),function(j) {
    if (j==1L) return(rep(TRUE,nrow(membership)))
    rowSums(membership[,seq_len(j-1L),drop=FALSE]==membership[,j])==0
  })
  for (j in seq_len(draws)) {
    counts <- integer(length(chrono))
    counts[chrono] <- .st3_moving(length(chrono),width)
    t <- counts[windows]
    moving[,j] <- evaluate(t)
    samples <- tabulate(sample.int(length(ids),length(ids),replace=TRUE),length(ids))
    space <- rep(1,nrow(membership))
    for (i in seq_len(ncol(membership))) {
      take <- unique_member[[i]]
      space[take] <- space[take]*samples[membership[take,i]]
    }
    crossed[,j] <- evaluate(t*space)
  }
  out <- cbind(spatial_block=spatial,temporal_block=apply(moving,1,sd,na.rm=TRUE),
               crossed=apply(crossed,1,sd,na.rm=TRUE))
  rownames(out) <- colnames(object$Y)
  list(estimate=drop(contrast %*% object$coefficients),se=out,
       finite_temporal=rowSums(is.finite(moving)),finite_crossed=rowSums(is.finite(crossed)))
}
