source("docs/phase5_refinement/common.R")
case <- as.integer(commandArgs(TRUE)[1L])
n <- 30L
s <- simulate_metacommunity(n,45,breadth=c(.28,.50)[case],abundance_mean=c(.7,1.3)[case],seed=628+case)
s$comm <- cbind(s$comm,diag(1L*(rowSums(s$comm)==0)))
x <- s$comm
env <- s$env
orders <- c(3L,4L)
z <- zeta_order_support(x,orders)
data <- function(M,sampler,exact_limit=0) ugdm_path_data(x,env,
  sample_sets_coupled(x,orders,M,sampler=sampler,exact_limit=exact_limit))
ref <- data(choose(n,4),"uniform",50000)
fitting <- function(a) {
  d <- ugdm_model_design(a$design,"B")
  fit <- ugdm_fit_variant(d,a$relative,model="B",sampling_weights=a$sampling_weights,free_sign=TRUE)
  theta <- c(fit$intercept,fit$beta,fit$lambda*d$k)
  list(fit=fit,theta=theta,design=d,y=a$relative,sampling_weights=a$sampling_weights)
}
reference <- lapply(ref,fitting)
scores <- function(a,theta) ugdm_pool_scores(list(a),theta,"conditional")
target <- Map(function(a,b) centers(scores(a,b)),reference,lapply(reference,function(s)s$theta))
theta <- lapply(reference,function(z)z$theta)
bench <- function(M,sampler) {
  start <- sum(proc.time()[1:2])
  d <- lapply(data(M,sampler),fitting)
  e <- Map(function(a,b,c) centers(scores(a,b))-c,d,theta,target)
  time <- sum(proc.time()[1:2])-start
  mods <- lapply(d,function(k) if(k$fit$converged) k$theta else rep(NA_real_,length(k$theta)))
  ref_s <- c(theta[[1]][2],theta[[2]][2],theta[[2]][2]-theta[[1]][2])
  ref_p <- c(theta[[1]][4],theta[[2]][4],theta[[2]][4]-theta[[1]][4])
  spread <- c(mods[[1]][2],mods[[2]][2],mods[[2]][2]-mods[[1]][2])
  position <- c(mods[[1]][4],mods[[2]][4],mods[[2]][4]-mods[[1]][4])
  err <- list(e[[1]],e[[2]],e[[2]]-e[[1]])
  true <- list(target[[1]],target[[2]],target[[2]]-target[[1]])
  data.frame(effect=c("lower","higher","paired"),spread_error=unname(spread-ref_s),
    position_error=unname(position-ref_p),
    score_error_ratio=vapply(seq_along(err),function(j) sqrt(sum(err[[j]]^2)/sum(true[[j]]^2)),0.0),
    CPU=time)
}
samplers <- c("uniform","site","species")
set.seed(8091+case)
base <- vapply(samplers,function(s) median(replicate(8,bench(2000,s)$CPU[1])),0.0)
budget <- setNames(pmax(500L,round(2000*base["uniform"]/base)),samplers)
out <- list()
for (s in samplers) for (j in 1:32) {
  a <- bench(budget[s],s)
  out[[length(out)+1L]] <- cbind(case,sampler=s,M=unname(budget[s]),r=j,
    target_CPU=unname(base["uniform"]),high_sharing_lower=z$sharing_lower[2],
    high_sharing_upper=z$sharing_upper[2],a)
}
write_refinement(do.call(rbind,out),paste0("equal_cpu_",case))
