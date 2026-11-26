source("docs/order_information/design.R")
jobs <- expand.grid(replicate=1:3,n=c(70L,150L,250L),name=c("narrow","dominant"))
orders <- 2:30
run <- lapply(seq_len(nrow(jobs)),function(i) {
  j <- jobs[i,]
  a <- community(j$n,j$name,59001L+i)
  b <- community(1600L,j$name,89012L+i)
  start <- proc.time()[["elapsed"]]
  profile <- zeta_order_information(a$x,a$env,orders=orders,response="relative",
    method="sample",combinations=500,pools=2,seed=3921L+i)
  screening <- proc.time()[["elapsed"]]-start
  start <- proc.time()[["elapsed"]]
  z <- lapply(orders,function(k) {
    x <- build(a,k,500L,"relative",16082L+i+k)
    p <- build(a,k,500L,"relative",46082L+i+k)
    v <- build(b,k,1500L,"relative",76082L+i+k)
    g <- fit_game(x)
    check <- evaluate_game(g,v)
    res <- tryCatch(evaluate_game(g,v,score_draws(x,p,g,seq_len(j$n),39L)),error=identity)
    list(x=x,p=p,v=v,g=g,check=check,res=res)
  })
  inference <- proc.time()[["elapsed"]]-start
  p <- z[[1L]]
  f <- ugdm_fit_interactions(p$x$design,p$x$y,ridge=.001,distance=FALSE)
  public <- ugdm_within_order_importance(f,p$x$design,p$x$y,
    evaluation=p$v$design,evaluation_y=p$v$y)$predictors$absolute
  original <- p$check[1L,names(public)]
  stopifnot(max(abs(original-public)) < 1e-5)
  numerical <- rbind(original,
    evaluate_game(fit_game(p$p),p$v)[1L,names(public)],
    evaluate_game(fit_game(build(a,2L,500L,"relative",9917L+i)),p$v)[1L,names(public)])
  confirmation <- if(is.matrix(p$res)) ugdm_within_order_diagnostic(f,p$x$design,p$x$y,
    integration=numerical, draws=p$res[,names(public),drop=FALSE],heldout_D2=p$check[1L,"D2"],
    min_response_sd=0,min_heldout_D2=.025,max_mc_se=.04,max_width=.20,max_condition=1e8,
    calibrated=FALSE)$status else "score_failed"
  data.frame(j,orders=length(orders),sets_per_pool=500,score_draws=39,
    screening_seconds=screening,inference_seconds=inference,ratio=inference/screening,
    exact_game_error=max(abs(original-public)),pre_fit_candidates=sum(profile$candidate),
    confirmation)
})
write_result(do.call(rbind,run),file.path(output,"benchmark.csv"))
