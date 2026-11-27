library(zetaGDUM)
outdir <- "docs/final_gdum/results"
mechanisms <- data.frame(name=c("occurrence","abundance","contrasting","null","weak","correlated","spatial"),
  occurrence=c(1,0,1,0,.23,.75,.7), abundance=c(0,1,1,0,.22,.75,.7),
  context=c(0,0,1,0,0,0,1), occupancy=c(-.45,1.3,.55,.35,.05,.90,.65),
  correlation=c(0,0,0,0,0,.82,.35), spatial=c(0,0,0,0,0,0,.60))
models <- c("spread_only", "additive", "interaction")
architectures <- c("coherent", "independent")
responses <- c("total", "occupancy", "heterogeneity")
scenario_id <- function(s) match(s, mechanisms$name)
environment <- function(n, scenario, seed) {
  set.seed(seed)
  p <- mechanisms[scenario_id(scenario),]
  z <- matrix(rnorm(n*2L), n)
  if (p$spatial>0) for (j in 2:n) {
    z[j,] <- p$spatial*z[j-1L,]+sqrt(1-p$spatial^2)*z[j,]
  }
  z[,2] <- p$correlation*z[,1]+sqrt(1-p$correlation^2)*z[,2]
  data.frame(phosphorus=10+75*pnorm(z[,1]), secondary=pnorm(z[,2]))
}
sites <- function(env, scenario, seed) {
  set.seed(seed)
  d <- mechanisms[scenario_id(scenario),]
  n <- nrow(env)
  s <- 48L
  x <- (env$phosphorus-10)/75
  opt <- seq(.03,.97,length.out=s)
  # Opposing changes in niche width govern occurrence and positive abundance.
  # Site normalisation and missing taxa then couple the observed components.
  w_occ <- .21*exp(2*d$context*(x-.5))
  w_ab <- .23*exp(-2*d$context*(x-.5))
  distance <- outer(x,opt,"-")
  niche_occ <- (distance/w_occ)^2
  niche_ab <- (distance/w_ab)^2
  marginal <- d$occupancy+.85*qnorm((seq_len(s)-.5)/s)
  p <- plogis(sweep(-d$occurrence*niche_occ,2,marginal,"+"))
  composition <- exp(.5*matrix(rnorm(n*s),n) - .7*d$abundance*niche_ab) *
    (matrix(runif(n*s),n) < p)
  cbind(1, composition)
}
generate <- function(n,scenario,seed) {
  e <- environment(n,scenario,seed)
  list(e=e,x=sites(e,scenario,seed+991L))
}
atomic_csv <- function(x,name) {
  target <- file.path(outdir,name)
  tmp <- tempfile(tmpdir=outdir)
  con <- if (endsWith(name,".gz")) gzfile(tmp,"wt") else file(tmp,"wt")
  write.csv(x,con,row.names=FALSE)
  close(con)
  stopifnot(file.rename(tmp,target))
}
reference_dir <- Sys.getenv("MS_GDUM_REFERENCE_DIR", file.path(outdir,"reference_cache"))
read_target <- function(s,k) readRDS(file.path(reference_dir,
                                             paste0(s,"_",k,".rds")))
sampling <- function(n,k,M,seed) sample_sets(n,k,M,seed)
truth_data <- function(s,k,sets=900L, communities=140L) {
  id <- scenario_id(s)
  e <- environment(480L,s,5101L+id)
  idx <- sampling(nrow(e),k,sets,8431L+k)
  expected <- Reduce("+", lapply(seq_len(communities),function(j) {
    unclass(ugdm_component_response(sites(e,s,36011L+id*1000L+j),idx))
  }))/communities
  class(expected) <- c("ugdm_component_response","matrix","array")
  attr(expected,"sets") <- idx
  observed <- ugdm_component_response(sites(e,s,62022L+id),idx)
  list(e=e,idx=idx,expected=expected,observed=observed)
}
importance <- function(f,d,y,z,t,resolution="predictor",permutations=48L) {
  a <- ugdm_component_importance(f,d,y,z,t,resolution=resolution,
      permutations=permutations,exact_max=if(length(f$vars)>3L) 0L else 3L,seed=if(length(f$vars)>3L) 10021L)
  x <- a$predictors[order(a$predictors$response,a$predictors$predictor),]
  list(point=x,hierarchy=a$hierarchy,failures=length(a$diagnostics$failed_coalitions))
}
accuracy <- function(p,t) vapply(responses,function(v) sqrt(mean((p[[v]]-t[,v])^2)),numeric(1))
