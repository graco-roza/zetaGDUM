source("docs/order_information/design.R")
args <- commandArgs(TRUE)
name <- args[1L]
n <- as.integer(args[2L])
split <- match.arg(args[3L], c("train","test","principal","budget","direct","validation","assessment"))
response <- match.arg(args[4L], c("total","relative"))
B <- as.integer(args[5L])
draws <- as.integer(args[6L])
M <- as.integer(args[7L])
orders <- if (length(args) >= 8L) as.integer(strsplit(args[8L],",")[[1L]]) else order_grid(name)
workers <- as.integer(Sys.getenv("ZETA_ORDER_WORKERS", "1"))
id <- scenario_index(name)
stopifnot(!is.na(id), all(orders <= n), B >= 1, draws >= 10, M >= 2)
# Separate ecological training and test samples. The budget audit pairs principal seeds.
phase <- c(train=1L,test=2L,principal=3L,budget=3L,direct=3L,validation=4L,assessment=3L)[[split]]
heldout <- community(1800L,name, 70000L+id*100L)
assessment_sets <- if (split == "assessment") 10000L else 1500L
target <- lapply(orders,function(k) build(heldout,k,assessment_sets,response,4113L+k))
start <- proc.time()[["elapsed"]]
one <- function(i) {
  seed <- 811L+id*10000000L + phase*1000000L+n*2000L+i
  a <- community(n,name,seed)
  groups <- if (name=="spatial") a$blocks else seq_len(n)
  ans <- lapply(seq_along(orders),function(j) {
    k <- orders[j]
    timed <- proc.time()[["elapsed"]]
    kseed <- 883L+k*110L+17L*i
    train <- build(a,k,M,response,kseed)
    repeat_pool <- build(a,k,M,response,kseed+900123L)
    game <- fit_game(train)
    point <- evaluate_game(game,target[[j]])
    repeat_fit <- fit_game(repeat_pool)
    repeat_point <- evaluate_game(repeat_fit,target[[j]])
    set.seed(kseed+71011L)
    uncertainty <- tryCatch({
      if (split == "direct") {
        t(vapply(seq_len(draws), function(r) {
          w <- site_weights(n,1L)
          sample <- build(a,k,M,response,idx=sample_sets_weighted(w,k,M))
          refit <- fit_game(sample)
          if (!all(vapply(refit$fits,function(f) f$converged,logical(1)))) return(rep(NA_real_,4))
          drop(evaluate_game(refit,target[[j]]))
        }, numeric(4), USE.NAMES=FALSE)) |>
          structure(dimnames=list(NULL,colnames(point)))
      } else evaluate_game(game,target[[j]],score_draws(train,repeat_pool,game,groups,draws))
    }, error=identity)
    base <- data.frame(name,n,split,response,replicate=i,seed,order=k,M,draws,
      gradient=point[,"gradient"], secondary=point[,"secondary"], difference=point[,"difference"],
      heldout=point[,"D2"], converged=all(vapply(game$fits,function(f) f$converged,logical(1))),
      repeat_converged=all(vapply(repeat_fit$fits,function(f) f$converged,logical(1))),
      mc_gradient=(point[,"gradient"]-repeat_point[,"gradient"])/sqrt(2),
      mc_secondary=(point[,"secondary"]-repeat_point[,"secondary"])/sqrt(2),
      mc_difference=(point[,"difference"]-repeat_point[,"difference"])/sqrt(2))
    if (inherits(uncertainty,"error")) {
      base$error <- conditionMessage(uncertainty)
      uncertainty <- matrix(NA_real_, 0, 4, dimnames=list(NULL,colnames(point)))
    } else base$error <- ""
    good <- complete.cases(uncertainty)
    base$valid_draws <- sum(good)
    for (v in c("gradient","secondary","difference")) {
      q <- if (sum(good) >= ceiling(.8*draws)) quantile(uncertainty[good,v],c(.10,.90,.025,.975), names=FALSE) else rep(NA_real_,4)
      base[paste0(v,c("_lo80","_hi80","_lo95","_hi95"))] <- as.list(q)
    }
    diff <- uncertainty[good,"difference"]
    base$p_greater <- if (length(diff)) mean(diff>0) else NA_real_
    base$p_meaningful <- if (length(diff)) mean(diff>.02) else NA_real_
    base$p_reverse_meaningful <- if (length(diff)) mean(diff<-.02) else NA_real_
    base$seconds <- proc.time()[["elapsed"]]-timed
    base
  })
  if (i %% 10L == 0L) message(name," ",n," ",split,": replicate ",i,"/",B)
  do.call(rbind,ans)
}
run <- if (workers>1L && .Platform$OS.type!="windows") {
  parallel::mclapply(seq_len(B),one,mc.cores=workers,mc.set.seed=FALSE)
} else lapply(seq_len(B),one)
stopifnot(all(vapply(run,is.data.frame,logical(1))))
message("Elapsed including evaluation prep and fits: ", proc.time()[["elapsed"]]-start)
write_result(do.call(rbind,run),file.path(output,paste(name,n,split,response,"calibration.csv.gz",sep="_")))
