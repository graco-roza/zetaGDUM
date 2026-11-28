st5_folds <- function(z, method=c("spatial","temporal","forward"), blocks=3L,
                      site_blocks=NULL, smallest=3L) {
  method <- match.arg(method)
  r <- z$rows
  if(!nrow(r)) return(list())
  if (method=="spatial") {
    sites <- sort(unique(unlist(r$sites)))
    if(is.null(site_blocks)) stop("Supply named spatial block memberships.")
    if(anyNA(site_blocks[sites])) stop("All observed sites need a spatial block.")
    labels <- split(sites,site_blocks[sites])
    if(length(labels)<2L) return(list())
    row <- function(group) {
      n <- vapply(r$sites,function(s) sum(s%in%group),integer(1))
      list(train=which(n==0L),test=which(n==r$k))
    }
  } else {
    dates <- sort(unique(unlist(r$dates)))
    if(length(dates)<2L*smallest) return(list())
    blocks <- min(blocks,floor(length(dates)/smallest))
    labels <- split(dates,cut(seq_along(dates),blocks,labels=FALSE))
    if(method=="forward") labels <- labels[-1L]
    row <- function(group) {
      n <- vapply(r$dates,function(s) sum(s%in%group),integer(1))
      train <- if(method=="forward") which(r$final<min(group)) else which(n==0L)
      list(train=train,test=which(n==lengths(r$dates)))
    }
  }
  folds <- lapply(labels,row)
  Filter(function(f) length(f$train)>=4L && length(f$test)>=1L,folds)
}

.st5_fold_data <- function(p,z,f,predictors,mode) {
  tr <- st5_subset(z,f$train)
  te <- st5_subset(z,f$test)
  ref <- st5_basis(p,tr,predictors)
  a <- st5_design(p,tr,ref)
  b <- st5_design(p,te,ref)
  i <- .st5_valid(tr,a,mode)
  j <- .st5_valid(te,b,mode)
  if(length(i)<4L || !length(j)) return(NULL)
  list(tr=tr,te=te,a=a,b=b,ref=ref,i=i,j=j)
}

.st5_loss <- function(z,pred,mode,i=seq_along(z$y)) {
  w <- .st5_weight_response(z,i)
  err <- (z$y[i]-pred$endpoint[i])^2
  out <- c(total=sum(w*.st5_ce(z$y[i],pred$endpoint[i])), total_change=sum(w*err))
  if(mode=="crossed") {
    q <- z$composition[i,3:4,drop=FALSE]
    out["composition"] <- -sum(rowSums(cbind(1-rowSums(q),q)*log(pmax(pred$composition[i,],1e-12)))*w)
  }
  if(!is.null(pred$phi)) {
    p <- colSums((z$phi[i,,drop=FALSE]-pred$phi[i,,drop=FALSE])^2*w)
    out <- c(out,p,joint=sum(p))
  }
  out
}

st5_tune <- function(p,z,folds,predictors,mode,grid=10^seq(-6,2),
                     equal_rows=FALSE,monotonic=FALSE) {
  grid <- sort(unique(grid))
  data <- Filter(Negate(is.null),lapply(folds,function(f) .st5_fold_data(p,z,f,predictors,mode)))
  if(!length(data)) return(list(total=.01,contrast=.01,inner_folds=0L))
  tune <- function(target,mode) {
    if(length(grid)==1L) return(grid)
    score <- vapply(grid,function(lambda) mean(vapply(data,function(d) {
      fit <- st5_fit(d$tr,d$a,lambda,lambda,mode,equal_rows=equal_rows,monotonic=monotonic)
      .st5_loss(d$te,st5_predict(fit,d$te,d$b),mode,d$j)[[target]]
    },numeric(1))),numeric(1))
    grid[which.min(score)]
  }
  total <- if(mode=="ridge") .01 else tune(if(mode=="crossed") "composition" else "total",
                                          if(mode=="crossed") "crossed" else "fractional")
  # Contrast tuning uses centered contrasts directly, independent of endpoint error.
  contrast <- if(mode=="fractional") total else {
    score <- vapply(grid,function(lambda) mean(vapply(data,function(d) {
      w <- .st5_weight_response(d$tr,d$i,equal_rows)
      prep <- .st5_prepare(d$a$x[d$i,,drop=FALSE],w)
      s <- d$b$x[d$j,names(prep$center),drop=FALSE]
      X <- cbind(1,sweep(sweep(s,2,prep$center),2,prep$scale,"/"))
      ww <- .st5_weight_response(d$te,d$j)
      score_one <- function(tr,te) sum(ww*rowSums((te-X%*%.st5_ridge(prep$X,tr,w,lambda))^2))
      centered <- function(z) z-rowMeans(z)
      if(mode=="crossed") return(sum(vapply(1:2,function(j)
        score_one(centered(d$tr$crossed[d$i,,j]),centered(d$te$crossed[d$j,,j])),numeric(1))))
      tr <- d$tr$phi[d$i,,drop=FALSE]; te <- d$te$phi[d$j,,drop=FALSE]
      if(mode=="coherent") score_one(centered(tr),centered(te)) else score_one(tr,te)
    },numeric(1))),numeric(1))
    grid[which.min(score)]
  }
  list(total=total,contrast=contrast,inner_folds=length(data))
}

st5_validate <- function(p,z,method=c("spatial","temporal","forward"),site_blocks=NULL,
                         mode=if(is.null(z$phi)) "fractional" else "coherent",
                         predictors=colnames(p$env),grid=10^seq(-6,2),blocks=3L,
                         equal_rows=FALSE,monotonic=FALSE,inner_blocks=NULL) {
  method <- match.arg(method)
  folds <- st5_folds(z,method,blocks,site_blocks)
  out <- lapply(names(folds),function(label) {
    d <- .st5_fold_data(p,z,folds[[label]],predictors,mode)
    if(is.null(d)) return(NULL)
    ib <- if(is.null(inner_blocks)) site_blocks else inner_blocks
    inner <- st5_folds(d$tr,method,blocks,ib)
    d$penalty <- st5_tune(p,d$tr,inner,predictors,mode,grid,equal_rows,monotonic)
    d$mode <- mode
    d$label <- label
    d$equal_rows <- equal_rows
    d$monotonic <- monotonic
    d
  })
  out <- Filter(Negate(is.null),out)
  if(!length(out)) return(list(folds=out,scores=data.frame(),method=method))
  baseline <- c("baseline","log_lag","initial_duration","change_duration")
  scores <- lapply(out,function(d) {
    options <- list(full=colnames(d$a$x),nuisance=intersect(baseline,colnames(d$a$x)),intercept=character())
    do.call(rbind,lapply(names(options),function(model) {
      fit <- .st5_refit(d,options[[model]])
      pred <- st5_predict(fit,d$te,d$b)
      loss <- .st5_loss(d$te,pred,mode,d$j)
      data.frame(method,fold=d$label,model,mode,train=fit$n,test=length(d$j),
        candidate_test=nrow(d$te$rows),rank=fit$rank,columns=length(fit$columns),
        unsupported=paste(fit$unsupported,collapse=";"),constant=paste(fit$prep$constant,collapse=";"),
        extrapolation=mean(pred$extrapolation[d$j]),lambda=d$penalty$total,
        lambda_contrast=d$penalty$contrast,inner_folds=d$penalty$inner_folds,
        invalid_endpoint=mean(pred$incompatible_endpoint[d$j]),
        invalid_process=if(is.null(pred$phi)) NA_real_ else mean(pred$incompatible_process[d$j]),
        convergence=if(is.null(fit$total)) NA_integer_ else fit$total$convergence,
        as.list(loss),rmse=sqrt(loss[["total_change"]]),check.names=FALSE)
    }))
  })
  list(folds=out,scores=do.call(rbind,scores),method=method)
}

.st5_refit <- function(d,columns=colnames(d$a$x))
  st5_fit(d$tr,d$a,d$penalty$total,d$penalty$contrast,d$mode,columns,
          d$equal_rows,d$monotonic)

st5_importance <- function(validation,by=c("predictor","group"),exact_limit=6L,
                           permutations=128L,seed=21L) {
  by <- match.arg(by)
  folds <- validation$folds
  if(!length(folds)) return(data.frame())
  terms <- function(d) d$a$columns
  vars <- sort(setdiff(unique(unlist(lapply(folds,function(d) terms(d)[[by]]))),
                       c("","nuisance")))
  if(!length(vars)) return(data.frame())
  p <- length(vars)
  old <- get0(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(is.null(old)) rm(".Random.seed",envir=.GlobalEnv)
          else assign(".Random.seed",old,envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  cached <- new.env(parent=emptyenv())
  evaluate <- function(set) {
    key <- paste0("S",paste(sort(set),collapse=","))
    if(!exists(key,cached,inherits=FALSE)) {
      f <- lapply(folds,function(d) {
        cols <- terms(d)
        keep <- cols$name[cols$group=="nuisance" | cols[[by]]%in%vars[set]]
        fit <- .st5_refit(d,keep)
        .st5_loss(d$te,st5_predict(fit,d$te,d$b),d$mode,d$j)
      })
      assign(key,do.call(rbind,f),cached)
    }
    get(key,cached,inherits=FALSE)
  }
  empty <- evaluate(integer())
  full <- evaluate(seq_len(p))
  out <- lapply(seq_len(nrow(empty)+1L),function(f) {
    test_loss <- function(x) if(f<=nrow(empty)) x[f,] else colMeans(x)
    div <- test_loss(empty)
    value <- function(z) 1-test_loss(evaluate(z))/ifelse(div>1e-14,div,NA_real_)
    phi <- matrix(0,p,ncol(empty))
    mc <- phi
    if(p<=exact_limit) {
      for(j in seq_len(p)) {
        other <- setdiff(seq_len(p),j)
        for(k in 0:length(other)) {
          sets <- if(k==0L) list(integer()) else combn(other,k,simplify=FALSE)
          for(s in sets) phi[j,] <- phi[j,]+(value(c(s,j))-value(s))/(p*choose(p-1,k))
        }
      }
    } else {
      draws <- array(NA_real_,c(p,ncol(empty),permutations))
      for(i in seq_len(permutations)) {
        order <- sample.int(p)
        old <- value(integer())
        for(j in seq_len(p)) {
          new <- value(order[seq_len(j)])
          draws[order[j],,i] <- new-old
          old <- new
        }
      }
      phi <- apply(draws,c(1,2),mean)
      mc <- apply(draws,c(1,2),sd)/sqrt(permutations)
    }
    improvement <- 1-test_loss(full)/ifelse(div>1e-14,div,NA_real_)
    if(max(c(0,abs(colSums(phi)-improvement)),na.rm=TRUE)>1e-10)
      stop("Predictive importance additivity failed.")
    data.frame(fold=if(f<=nrow(empty)) folds[[f]]$label else "pooled",feature_type=by,
               feature=rep(vars,ncol(empty)),response=rep(colnames(empty),each=p),
               importance=as.vector(phi),mc_se=as.vector(mc),
               total_improvement=rep(improvement,each=p))
  })
  do.call(rbind,out)
}

# Resample intact held-out blocks. This measures sensitivity of validation summaries
# to the chosen blocks and does not supply independent-row inferential intervals.
st5_block_stability <- function(importance,draws=499L,seed=23L) {
  u <- importance[importance$fold!="pooled",,drop=FALSE]
  if(!nrow(u)) return(data.frame())
  old <- get0(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(is.null(old)) rm(".Random.seed",envir=.GlobalEnv)
          else assign(".Random.seed",old,envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  by <- split(u,interaction(u$feature_type,u$response,drop=TRUE))
  do.call(rbind,lapply(by,function(z) {
    wide <- xtabs(importance~fold+feature,z,na.action=na.pass)
    if(nrow(wide)<2L || any(!is.finite(z$importance))) return(NULL)
    sampled <- replicate(draws,colMeans(wide[sample.int(nrow(wide),replace=TRUE),,drop=FALSE]))
    top <- colnames(wide)[apply(sampled,2,which.max)]
    q <- t(apply(sampled,1,quantile,c(.025,.975)))
    data.frame(response=z$response[1],feature_type=z$feature_type[1],feature=colnames(wide),
       blocks=nrow(wide),mean_fold_importance=colMeans(wide),lower=q[,1],upper=q[,2],
       top_frequency=as.numeric(table(factor(top,levels=colnames(wide))))/draws)
  }))
}
