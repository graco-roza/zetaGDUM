source("docs/phase5_refinement/common.R")
arg <- commandArgs(TRUE)
meta <- read_manifest()[as.integer(arg[2]),]
community <- load_dd2(meta,arg[1])
x <- community$x
env <- community$env
geo <- community$dist
n <- nrow(x)
orders <- unlist(meta[c("mid","high")],use.names=FALSE)
component <- "relative"
basis <- suppressWarnings(ugdm_design(env,cbind(seq_len(n-1),seq.int(2,n)),geo,dist_range=range(as.dist(geo))))
types <- c("spread","position","relative_position")
names9 <- c(paste0("lower_",types),paste0("higher_",types),paste0("paired_",types))
tol <- setNames(rep(c(.015,.030,.050),3),names9)
q <- (quantile(env[[1]],c(.25,.75),names=FALSE)-basis$environment_basis$center[1]) /
  basis$environment_basis$scale[1]
span <- .1*colSums(diff(predict(basis$environment_basis$splines[[1]],q)))
paths <- function(M,sampler) ugdm_path_data(x,env,
  sample_sets_coupled(x,orders,M,sampler=sampler,exact_limit=50000),geo,ref=basis)
model <- function(a) {
  b <- model_data(a,component,"A")
  d <- b$design
  D <- cbind(d$Xdir,d$Xdist)
  X <- cbind(1,D,d$Xpos/d$k)
  p <- ncol(D)
  spread <- c(1L,seq_len(p)+1L)
  pos <- c(1L,p+1L+seq_len(ncol(d$Xpos)))
  ind <- list(seq_len(ncol(X)),spread,pos)
  con <- list(spread[-1L],spread[-1L],integer())
  P <- lapply(ind,function(i) matrix(0,length(i),length(i)))
  y <- b$y
  w <- rep(1,length(y))
  solve_fit <- internal(".variant_solution")
  evaluate <- function(keep = seq_along(y), warm = NULL) {
    f <- lapply(seq_along(ind),function(j) solve_fit(X[keep,ind[[j]],drop=FALSE],
      y[keep],w[keep],P[[j]],con[[j]],if(!is.null(warm)) warm[[j]]))
    a <- f[[1L]]
    stat <- if (a$constant || a$null_deviance<=0 || !a$converged) rep(NA_real_,3) else {
      full <- 1-a$deviance/a$null_deviance
      e <- 1-f[[2]]$deviance/f[[2]]$null_deviance
      h <- 1-f[[3]]$deviance/f[[3]]$null_deviance
      c(sum(a$coef[2:4]*span),a$coef[pos[2L]]*diff(q)/d$k,
        if(all(vapply(f,function(g) g$converged,TRUE)) && full>1e-12) .5*(full+h-e)/full else NA_real_)
    }
    list(effect=stat, coefficients=lapply(f,function(g)g$coef))
  }
  idx <- d$idx
  list(evaluate=evaluate,idx=idx,X=X,y=y)
}
effects <- function(M,sampler) {
  d <- lapply(paths(M,sampler),model)
  a <- d[[1]]$evaluate()$effect
  b <- d[[2]]$evaluate()$effect
  setNames(c(a,b,b-a),names9)
}
high <- if (length(arg)>2L) as.integer(arg[3L]) else NA_integer_
if (is.na(high)) {
  message(meta$community," pilot")
  buy <- ugdm_contrast_budget(effects,c("uniform","site"),c(4000L,8000L,16000L),tol,
    nrep=6,max_budget=40000L,seed=912+n)
  if (is.null(buy$contrasts)) stop("No sampler gave a complete primary pilot.")
  tab <- buy$contrasts
  family <- buy$sampler
  M <- if (any(tab$passed)) max(tab$verified_at[tab$passed]) else NA_real_
  key <- meta$id
} else {
  # Keep the sampler selected in the original independent pilot.
  old_path <- file.path("docs/phase5_refinement/runs",paste0("emp_budget_",meta$id,".csv"))
  old <- if (file.exists(old_path)) read.csv(old_path) else
    subset(read.csv("docs/phase5_refinement/results/empirical_budget.csv.gz"),
      community==meta$community & tier=="point")
  family <- unique(old$sampler)
  stopifnot(length(family)==1L)
  tab <- old[,c("effect","tolerance","A","slope","predicted","candidate","verified_at",
                  "confirmed_sd","confirmed_upper","passed")]
  M <- high
  key <- paste0(meta$id,"_high")
}
id <- data.frame(community=meta$community,n,component,order1=orders[1],order2=orders[2],
  sampler=family,tier=if (is.na(high)) "point" else "larger")
paper <- cbind(id,tab,reference_M=M)
if (!is.finite(M)) {
  write_refinement(paper,paste0("emp_budget_",key))
} else {
  set.seed(1209+n+if(is.na(high)) 0L else 900000L)
  confirmation <- replicate(10L,effects(M,family))
  U <- apply(confirmation,1,sd)*sqrt(9/qchisq(.1,9))
  paper$fresh_upper <- unname(U[tab$effect])
  paper$fresh_pass <- as.logical(is.finite(U[tab$effect]) & U[tab$effect]<=tol[tab$effect])
  write_refinement(paper,paste0("emp_budget_",key))
  message(meta$community," parallel-pool matched deletion reference M=",M)
  if (any(paper$fresh_pass)) {
    jackknife <- function(pool) {
      d <- lapply(pool,model)
      fits <- lapply(d,function(a)a$evaluate())
      deletes <- lapply(seq_along(d),function(j) {
        a <- d[[j]]
        idx <- a$idx
        row <- seq_len(nrow(idx))
        sites <- split(rep(row,ncol(idx)),as.vector(idx))
        t(vapply(seq_len(n),function(i)
          a$evaluate(setdiff(row,sites[[as.character(i)]]),fits[[j]]$coefficients)$effect,numeric(3)))
      })
      p <- c(fits[[1]]$effect,fits[[2]]$effect,fits[[2]]$effect-fits[[1]]$effect)
      de <- cbind(deletes[[1]],deletes[[2]],deletes[[2]]-deletes[[1]])
      minimum <- vapply(d,function(a)nrow(a$idx)-max(tabulate(a$idx,n)),integer(1))
      list(point=p, deletions=de, min_exclusion=min(minimum))
    }
    first <- jackknife(paths(M,family))
    second <- jackknife(paths(M,family))
    keep <- colSums(is.finite(first$deletions)&is.finite(second$deletions))==n &
      is.finite(first$point)&is.finite(second$point)
    result <- lapply(which(keep),function(j) {
      v <- n*first$point[j]-(n-1)*first$deletions[,j,drop=FALSE]
      w <- n*second$point[j]-(n-1)*second$deletions[,j,drop=FALSE]
      cov <- ugdm_pseudovalue_noise(v,w)
      pool <- (v+w)/2
      average <- (first$point[j]+second$point[j])/2
      pseudo <- ugdm_pseudovalue_bootstrap(average,(n*average-pool)/(n-1),draws=999L)$draws[,1L]
      data.frame(id,M,effect=names9[j],estimate=average,
        first_SE=unname(sd(v[,1])/sqrt(n)),second_SE=unname(sd(w[,1])/sqrt(n)),
        pooled_SE=sqrt(cov$pooled_jackknife[1,1]),numerical_fraction=cov$numerical_fraction,
        pseudo_lo=quantile(pseudo,.025,names=FALSE),pseudo_hi=quantile(pseudo,.975,names=FALSE),
        verified=paper$fresh_pass[j],min_exclusion=min(first$min_exclusion,second$min_exclusion))
    })
    if(length(result)) write_refinement(do.call(rbind,result),paste0("emp_uncert_",key))
  }
}
