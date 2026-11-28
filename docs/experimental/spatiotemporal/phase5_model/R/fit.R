st5_weights <- function(group, importance=rep(1,length(group)), equal_rows=FALSE) {
  if (!length(group) || anyNA(group) || length(importance)!=length(group) ||
      any(!is.finite(importance) | importance<=0)) stop("Invalid weights or empty fit.")
  w <- if(equal_rows) rep(1,length(group)) else importance/ave(importance,group,FUN=sum)
  w/sum(w)
}

.st5_weight_response <- function(z,i=seq_along(z$y),equal_rows=FALSE) {
  importance <- if(is.null(z$importance_weight)) rep(1,length(i)) else z$importance_weight[i]
  st5_weights(z$weight_group[i],importance,equal_rows)
}

.st5_softplus <- function(z) pmax(z,0) + log1p(exp(-abs(z)))
.st5_ce <- function(y,p) {
  p <- pmin(pmax(p,1e-12),1-1e-12)
  -y*log(p)-(1-y)*log1p(-p)
}

.st5_prepare <- function(x,w) {
  center <- colSums(x*w)
  centered <- sweep(x,2,center)
  scale <- sqrt(colSums(centered^2*w))
  keep <- scale > 1e-10
  X <- cbind(intercept=1,sweep(centered[,keep,drop=FALSE],2,scale[keep],"/"))
  list(X=X,center=center[keep],scale=scale[keep],constant=colnames(x)[!keep],rank=qr(X)$rank)
}

.st5_fractional <- function(X,y,w,lambda,monotonic=FALSE) {
  if(all(y==0) || all(y==1))
    return(list(coef=rep(0,ncol(X)),constant=y[1],convergence=0))
  penalty <- c(0,rep(lambda,ncol(X)-1L))
  calc <- function(b,gradient=FALSE) {
    eta <- drop(X%*%b)
    if(gradient) return(drop(crossprod(X,w*(plogis(eta)-y)))+2*penalty*b)
    sum(w*(.st5_softplus(eta)-y*eta)) + sum(penalty*b^2)
  }
  b <- c(qlogis(sum(w*y)),rep(0,ncol(X)-1L))
  if(any(monotonic)) {
    fit <- optim(b,calc,function(b) calc(b,TRUE),method="L-BFGS-B",
                 lower=ifelse(monotonic,0,-Inf),control=list(factr=1e5,maxit=400))
  } else fit <- optim(b,calc,function(b) calc(b,TRUE),method="BFGS",
                      control=list(reltol=1e-10,maxit=400))
  list(coef=fit$par,convergence=fit$convergence)
}

.st5_ridge <- function(X,Y,w,lambda) {
  # The baseline is a slope, so it receives the same penalty as all other slopes.
  solve(crossprod(X,X*w)+diag(c(0,rep(lambda,ncol(X)-1L)),ncol(X)),crossprod(X,Y*w))
}

.st5_softmax <- function(eta) {
  exp_eta <- exp(eta-pmax(eta[,1],eta[,2],eta[,3]))
  exp_eta/rowSums(exp_eta)
}

.st5_multinomial <- function(X,y,w,lambda) {
  penalty <- c(0,rep(lambda,ncol(X)-1L))
  calc <- function(theta,gradient=FALSE) {
    b <- matrix(theta,ncol(X),2)
    eta <- cbind(0,X%*%b)
    prob <- .st5_softmax(eta)
    if(gradient) return(as.vector(crossprod(X,(prob-y)[,-1,drop=FALSE]*w)+2*penalty*b))
    z <- pmax(eta[,1],eta[,2],eta[,3])
    sum(w*(z+log(rowSums(exp(eta-z)))-rowSums(y*eta)))+sum(penalty*b^2)
  }
  avg <- colSums(y*w)
  active <- avg>0
  if(!all(active)) {
    # A category absent throughout training is a structural zero in this working fit.
    if(sum(active)==1) return(list(active=active,fit=NULL,convergence=0))
    fit <- .st5_fractional(X,y[,which(active)[2]],w,lambda)
    return(list(active=active,fit=fit,convergence=fit$convergence))
  }
  b <- matrix(0,ncol(X),2)
  b[1,] <- log(avg[-1]/avg[1])
  fit <- optim(as.vector(b),calc,function(b) calc(b,TRUE),method="BFGS",
               control=list(reltol=1e-10,maxit=500))
  list(coef=matrix(fit$par,ncol(X),2),active=active,convergence=fit$convergence)
}

.st5_valid <- function(z,design,mode) {
  good <- !design$missing & is.finite(z$baseline) & is.finite(z$y)
  if(mode %in% c("coherent","ridge","crossed")) good <- good & complete.cases(z$phi)
  if(mode=="crossed") good <- good & complete.cases(z$composition) &
    apply(is.finite(z$crossed),1,all)
  which(good)
}

st5_fit <- function(z,design,lambda=.01,lambda_contrast=lambda,
                    mode=if(is.null(z$phi)) "fractional" else "coherent",
                    columns=colnames(design$x),equal_rows=FALSE,monotonic=FALSE) {
  mode <- match.arg(mode,c("fractional","coherent","ridge","crossed"))
  if (any(!is.finite(c(lambda,lambda_contrast)) | c(lambda,lambda_contrast)<=0))
    stop("Use positive and finite ridge penalties.")
  if (monotonic && (z$type!="spatial" || mode!="fractional"))
    stop("The monotonic sensitivity applies to static spatial dissimilarity.")
  valid <- .st5_valid(z,design,mode)
  z <- st5_subset(z,valid)
  w <- .st5_weight_response(z,equal_rows=equal_rows)
  if(any(z$y<0 | z$y>1 | z$baseline<0 | z$baseline>1))
    stop("Use normalized responses.")
  prep <- .st5_prepare(design$x[valid,columns,drop=FALSE],w)
  X <- prep$X
  result <- list(mode=mode,prep=prep,lambda=lambda,lambda_contrast=lambda_contrast,
                 process=colnames(z$phi),n=length(valid),valid=valid,
                 rank=prep$rank,columns=colnames(X),unsupported=design$unsupported)
  delta <- z$y-z$baseline
  if (mode!="fractional" && max(abs(rowSums(z$phi)-delta))>1e-10) stop("Invalid observed attribution.")
  center <- function(phi,change) {
    b <- .st5_ridge(X,phi-change/ncol(phi),w,lambda_contrast)
    if(max(abs(rowSums(b)))>1e-10) stop("The fitted contrast coefficients must sum to zero.")
    b
  }
  if(mode %in% c("fractional","coherent")) {
    bounded <- if(monotonic) colnames(X)%in%design$columns$name[design$columns$group=="spatial_spread"] else FALSE
    result$total <- .st5_fractional(X,z$y,w,lambda,bounded)
  }
  if(mode=="coherent") result$contrast <- center(z$phi,delta)
  if(mode=="ridge") result$contrast <- .st5_ridge(X,z$phi,w,lambda_contrast)
  if(mode=="crossed") {
    q0 <- z$composition[,1:2,drop=FALSE]
    q1 <- z$composition[,3:4,drop=FALSE]
    if(max(abs(apply(z$crossed,c(1,3),sum)-(q1-q0)),
           abs(apply(z$crossed,c(1,2),sum)-z$phi))>1e-10)
      stop("Invalid observed crossed attribution.")
    y <- cbind(1-rowSums(q1),q1)
    y <- pmax(y,0)
    result$total <- .st5_multinomial(X,y/rowSums(y),w,lambda)
    result$crossed <- lapply(1:2,function(j) center(z$crossed[,,j],q1[,j]-q0[,j]))
  }
  result
}

st5_predict <- function(fit,z,design) {
  a <- design$x[,names(fit$prep$center),drop=FALSE]
  X <- cbind(intercept=1,sweep(sweep(a,2,fit$prep$center),2,fit$prep$scale,"/"))
  total <- fit$total
  baseline <- z$baseline
  missing <- design$missing | !is.finite(baseline)
  if(fit$mode=="crossed") missing <- missing | !complete.cases(z$composition[,1:2,drop=FALSE])
  X[missing,] <- NA_real_
  out <- list(extrapolation=design$extrapolation,missing=missing)
  if(fit$mode %in% c("fractional","coherent")) {
    endpoint <- if(is.null(total$constant)) plogis(drop(X%*%total$coef)) else
      ifelse(missing,NA_real_,total$constant)
  }
  if(fit$mode=="coherent") out$phi <- X%*%fit$contrast+(endpoint-baseline)/length(fit$process)
  if(fit$mode=="ridge") {
    out$phi <- X%*%fit$contrast
    endpoint <- baseline+rowSums(out$phi)
  }
  if(fit$mode=="crossed") {
    mu <- matrix(0,nrow(X),3)
    active <- total$active
    if(all(active)) mu <- .st5_softmax(cbind(0,X%*%total$coef))
    else if(sum(active)==1L) mu[,active] <- 1
    else {
      prob <- plogis(drop(X%*%total$fit$coef))
      mu[,active] <- cbind(1-prob,prob)
    }
    mu[missing,] <- NA_real_
    endpoint <- rowSums(mu[,-1,drop=FALSE])
    q0 <- z$composition[,1:2,drop=FALSE]
    change <- mu[,-1,drop=FALSE]-q0
    crossed <- lapply(1:2,function(j) X%*%fit$crossed[[j]]+change[,j]/length(fit$process))
    out$phi <- crossed[[1]]+crossed[[2]]
    out$crossed <- array(unlist(crossed),c(dim(out$phi),2L),
                         dimnames=list(NULL,fit$process,c("occupancy","numerical")))
    out$composition <- mu
  }
  if(!is.null(out$phi)) {
    colnames(out$phi) <- fit$process
    out$incompatible_process <- rowSums(abs(out$phi)>1+1e-10)>0
    if(!is.null(out$crossed)) out$incompatible_process <-
      out$incompatible_process | apply(abs(out$crossed)>1+1e-10,1,any)
  }
  out$endpoint <- endpoint
  out$change <- endpoint-baseline
  out$incompatible_endpoint <- !is.na(endpoint) & (endpoint<0 | endpoint>1)
  out
}
