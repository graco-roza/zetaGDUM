st3_scenarios <- c("expansion","contraction","redistribution","numerical",
                   "expansion_numerical","contraction_numerical","opposing",
                   "filtering","widespread","restricted","persistence","null","shared_shock")

st3_pair_windows <- function(times, m = 1L, stride = 2L*m, offset = 1L, gap = 0L) {
  m <- .st_one(m,length(times)%/%2L,"m")
  if (stride < 1L || offset < 1L || gap < 0L || offset + 2*m + gap - 1L > length(times)) {
    stop("No paired windows within the declared anchor and gap.")
  }
  lapply(seq.int(offset,length(times)-2L*m-gap+1L,by=stride),function(i) {
    list(before=i+seq_len(m)-1L,after=i+m+gap+seq_len(m)-1L)
  })
}

st3_simulate <- function(scenario = "null", n = 20L, surveys = 16L, seed = 1L,
                         env_seed = 314L, irregular = FALSE) {
  scenario <- match.arg(scenario,c(st3_scenarios,"mixed_allocation"))
  if (n<8L || surveys<8L || surveys%%2L!=0L) stop("Use at least eight sites and an even number of surveys >= 8.")
  S <- 12L
  sites <- paste0("s",seq_len(n))
  blocks <- rep(seq_len(5L),length.out=n)
  times <- if (irregular) cumsum(rep(c(1,2,1,3,1),length.out=surveys))-1 else seq_len(surveys)
  labels <- list(sites,as.character(times),paste0("sp",seq_len(S)))
  set.seed(env_seed)
  env <- array(rnorm(n*surveys*2L),c(n,surveys,2L),list(sites,labels[[2]],c("climate","soil")))
  for (j in seq(2L,surveys,2L)) {
    env[,j-1L,1L] <- env[,j-1L,1L]*.65 + (blocks-3)*.35
    env[,j,1L] <- env[,j-1L,1L] + .75*env[,j,1L] + rnorm(1)
    env[,j,2L] <- .7*env[,j-1L,2L]+sqrt(1-.7^2)*env[,j,2L]
  }
  set.seed(seed)
  out <- array(0,c(n,surveys,S),labels)
  site_memory <- rnorm(n,0,.45)
  for (j in seq(1L,surveys-1L,2L)) {
    base <- env[,j,1L]
    change <- env[,j+1L,1L]-base
    signal <- change + .25*env[,j,2L]
    if (scenario %in% c("null","shared_shock")) signal[] <- 0
    shock <- if (scenario=="shared_shock") rnorm(1,0,2.1)+rnorm(5,0,1.2)[blocks] else 0
    propensity <- plogis(-.65+1.35*signal+site_memory+shock)
    if (scenario=="opposing") {
      type <- base + .45*change + rnorm(1,0,.3) > 0
      magnitude <- runif(1,.33,.67)
      a <- matrix(0,n,S)
      b <- a
      z <- 2*magnitude/(1-magnitude)
      a[,1:2] <- 1
      a[type,3] <- z
      common <- z + runif(1,1.0,3.2)
      separation <- magnitude*(common+z)
      hi <- (common+separation)/2
      lo <- (common-separation)/2
      if (lo<=0) {
        common <- z*magnitude/(1-magnitude) + 2
        hi <- (common+magnitude*(common+z))/2
        lo <- common-hi
      }
      b[,1] <- ifelse(type,hi,lo)
      b[,2] <- ifelse(type,lo,hi)
      b[,3] <- z
    } else {
      occupancy <- outer(rep(1,n),c(1,1,.63,.38,.18,.78,.52,.31,.7,.36,.20,.55))
      if (scenario %in% c("filtering","widespread","restricted")) {
        occupancy[,3:S] <- plogis(qlogis(occupancy[,3:S])+
                                   outer(base,rep(c(-1,1),length.out=S-2L)))
      }
      x <- matrix(runif(n*S),n,S)<occupancy
      amounts <- exp(matrix(rnorm(n*S,1.15,.45),n,S))
      a <- x*amounts
      b <- a
      if (scenario=="contraction_numerical") {
        a[,1] <- exp(1.5+1.0*base)
        a[,2] <- exp(1.5-1.0*base)
        b <- a
      }
      gain_species <- switch(scenario,widespread=c(6L,9L),restricted=c(5L,11L),3:5)
      if (scenario %in% c("expansion","expansion_numerical","filtering","widespread","restricted","null","shared_shock","mixed_allocation")) {
        for (s in gain_species) {
          add <- !x[,s] & runif(n)<propensity
          b[add,s] <- amounts[add,s]
        }
      }
      if (scenario %in% c("contraction","contraction_numerical","filtering","null","shared_shock")) {
        for (s in 6:8) b[x[,s] & runif(n)<propensity,s] <- 0
      }
      if (scenario %in% c("redistribution","filtering","mixed_allocation")) {
        for (s in if (scenario=="mixed_allocation") 3:5 else 9:12) {
          from <- which(x[,s])
          to <- which(!x[,s] & b[,s]==0)
          moves <- min(length(from),length(to),2L)
          if (moves) {
            give <- from[sample.int(length(from),moves,prob=plogis(-signal[from]))]
            receive <- to[sample.int(length(to),moves,prob=plogis(signal[to]))]
            b[give,s] <- 0
            b[receive,s] <- amounts[receive,s]
          }
        }
      }
      if (scenario %in% c("numerical","expansion_numerical","filtering","widespread","restricted","null","shared_shock")) {
        magnitude <- if (scenario %in% c("null","shared_shock")) .35*site_memory+shock/3 else .7*signal
        b[,1] <- b[,1]*exp(magnitude)
        b[,2] <- b[,2]*exp(-magnitude)
      }
      if (scenario=="contraction_numerical") b[,1:2] <- exp(1.5)
    }
    out[,j,] <- a
    out[,j+1L,] <- b
  }
  if (scenario=="persistence") {
    index <- seq_len(surveys) > surveys/2L
    reference <- out[,which(!index)[1L],]
    for (i in which(index)) {
      enable <- env[,i,1L] > -.5
      out[enable,i,] <- reference[enable,]
    }
  }
  list(panel=st_panel(out,times),env=env,blocks=blocks,scenario=scenario)
}

st3_sample_sets <- function(n,k,budget) {
  if (choose(n,k)<=budget) return(.st_indices(n,k,budget))
  t(matrix(vapply(seq_len(budget),function(i) sort(sample.int(n,k)),integer(k)),k,budget))
}

st3_collect <- function(panel, sets, windows, response = "abundance", crossed = FALSE,
                        exact_limit = 12L, draws = 8L, seed = 1L, retain_draws = FALSE) {
  wanted <- if (crossed) c("abundance","occupancy","numerical") else response
  fields <- lapply(seq_along(windows),function(j) {
    w <- windows[[j]]
    p <- st3_plan(panel,w$before,w$after)
    x <- st3_shapley(p,sets,wanted,exact_limit,draws,seed+j,retain_draws)
    Y <- st3_model_response(x,response,crossed)
    se <- if (retain_draws && !x$exact) {
      a <- vapply(x$draws,function(p) {
        colMeans(matrix(p[,response,,drop=FALSE],ncol=4L),na.rm=TRUE)
      },numeric(4L))
      sd_row <- apply(a,1,sd)/sqrt(ncol(a))
      max(c(0,sd_row[is.finite(sd_row)]))
    } else max(c(0,x$allocation_mc_se[is.finite(x$allocation_mc_se)]))
    list(y=Y,before=x$before[,response],after=x$after[,response],
         valid=complete.cases(Y),error=x$max_efficiency_error,exact=x$exact,max_allocation_mc=se)
  })
  list(y=do.call(rbind,lapply(fields,`[[`,"y")),before=unlist(lapply(fields,`[[`,"before")),
       after=unlist(lapply(fields,`[[`,"after")),valid=unlist(lapply(fields,`[[`,"valid")),
       error=max(vapply(fields,`[[`,numeric(1),"error")),
       exact=all(vapply(fields,`[[`,logical(1),"exact")),
       allocation_mc=max(vapply(fields,`[[`,numeric(1),"max_allocation_mc")))
}

st3_subset <- function(panel, sites, times=seq_along(panel$times)) {
  st_panel(panel$x[sites,times,,drop=FALSE],panel$times[times])
}
