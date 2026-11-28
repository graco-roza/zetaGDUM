simulate_annual <- function(scenario, seed, n=60L, species=200L, years=30L) {
  set.seed(seed)
  process <- rep("fixed",species)
  if (scenario %in% c("expansion","opposing","reversal")) process[11:60] <- "expansion"
  if (scenario %in% c("contraction","opposing")) process[61:110] <- "contraction"
  if (scenario == "shocks") process[11:110] <- "redistribution"
  p <- runif(species,.08,.75)
  p[1:10] <- 1
  p[11:60] <- .17
  p[61:110] <- .8
  x <- array(FALSE,c(n,years,species))
  x[,1,] <- matrix(runif(n*species),n) < rep(p,each=n)
  block <- rep(seq_len(6),length.out=n)
  shock <- numeric(years)
  for (t in seq.int(2,years)) {
    previous <- matrix(x[,t-1,],n)
    x[,t,] <- previous
    shock[t] <- .8 * shock[t-1] + rnorm(1)
    if (scenario == "stable_redistribution") {
      # A common permutation preserves every full-combination spatial statistic.
      x[,t,] <- previous[sample.int(n),,drop=FALSE]
      next
    }
    if (scenario == "reversal" && t > ceiling(years/2)) {
      x[,t,] <- x[,years-t+1L,]
      next
    }
    for (s in which(process != "fixed")) {
      occupied <- which(previous[,s])
      absent <- which(!previous[,s])
      active <- (t+s) %% 2 == 0L
      if (process[s] == "expansion" && active && length(absent)) {
        x[absent[sample.int(length(absent),1)],t,s] <- TRUE
      }
      if (process[s] == "contraction" && active && length(occupied)) {
        x[occupied[sample.int(length(occupied),1)],t,s] <- FALSE
      }
      if (process[s] == "redistribution") {
        size <- min(length(occupied),length(absent),max(0,round(2+shock[t])))
        if (size) {
          x[occupied[sample.int(length(occupied),size)],t,s] <- FALSE
          x[absent[sample.int(length(absent),size)],t,s] <- TRUE
        }
      }
    }
  }
  observed <- x
  if (scenario == "detection") {
    # Same static latent panel; detectability fluctuates across whole years and blocks.
    for (t in seq_len(years)) {
      detection <- plogis(1.7 + .8*shock[t]+.5*cos(block))
      detected <- matrix(runif(n*species),n) < detection
      observed[,t,] <- x[,t,] & detected
    }
  }
  expected <- switch(scenario, expansion="positive raw zeta change", contraction="negative raw zeta change",
    opposing="opposing raw zeta contributions", stable_redistribution="invariant spatial diversity with local turnover",
    reversal="zero endpoint change with nonzero cumulative gross change", shocks="autocorrelated regional reshuffling",
    detection="zero latent but nonzero apparent local turnover", null="zero change")
  list(panel=st_panel(observed), latent=st_panel(x), blocks=block, expected=expected, shock=shock)
}

simulation_rep <- function(sim, scenario, replicate, sample_budget = 192L) {
  x <- sim$panel$x
  n <- dim(x)[1]
  T <- dim(x)[2]
  sets <- lapply(1:6,function(k) sets_for(n,k,sample_budget,replicate*11L+k))
  previous <- lapply(sets,function(s) matrix(0,nrow(s),3))
  gross <- previous
  observed_annual <- true_annual <- list()
  temporal <- list()
  for (t in seq_len(T-1L)) {
    plan <- st3_plan(sim$panel,t,t+1)
    alloc <- st3_allocations(plan,draws=8,seed=replicate+t*199L)
    ledger <- plan$ledger
    o0 <- colSums(x[,t,]); o1 <- colSums(x[,t+1,])
    if (scenario != "detection") {
      stopifnot(all(ledger$gains == 0 | ledger$losses == 0 | ledger$gains == ledger$losses))
    }
    true_annual[[t]] <- rbindlist(lapply(1:6,function(k) {
      diff <- .st_choose_prob(o1,k,n)-.st_choose_prob(o0,k,n)
      data.table(year=t+1L,k=k,process=processes, truth=c(sum(diff[o1>o0]),sum(diff[o0>o1]),0))
    }))
    for (k in 1:6) {
      s <- sets[[k]]
      v <- evaluate(plan,s,alloc)
      stopifnot(v$error < 1e-10)
      value <- v$phi[,1,,drop=FALSE]
      dim(value) <- c(nrow(s),3L)
      previous[[k]] <- previous[[k]] + value
      gross[[k]] <- gross[[k]] + abs(value)
      observed_annual[[paste(t,k)]] <- data.table(year=t+1L,k=k,process=processes,
        estimate=colMeans(value), negative_magnitude=colMeans(pmax(-value,0)), positive_magnitude=colMeans(pmax(value,0)))
      if (k==1L) temporal[[t]] <- data.table(year=t+1L, local_temporal=mean_valid(v$after[,6]),
        latent_local=mean_valid(evaluate(st3_plan(sim$latent,t,t+1),s)$after[,6]),
        occupancy_events=sum(ledger$gains+ledger$losses), shock=sim$shock[t+1L])
    }
  }
  plan <- st3_plan(sim$panel,1L,T)
  alloc <- st3_allocations(plan,draws=8,seed=replicate+987L)
  a <- merge(rbindlist(observed_annual),rbindlist(true_annual),by=c("year","k","process"))
  summary <- rbindlist(lapply(1:6,function(k) {
    z <- evaluate(plan,sets[[k]],alloc)
    v <- previous[[k]]
    stopifnot(max(abs(rowSums(v)-(z$after[,1]-z$before[,1]))) < 1e-10)
    original <- z$phi[,1,,drop=FALSE]
    dim(original) <- dim(v)
    rows_at_order <- which(a$k == k)
    true_contrast <- a[rows_at_order, .(truth=sum(truth)),by=process]
    uncertainty <- lapply(seq_along(processes),function(p) {
      y <- v[,p]
      s <- sets[[k]]
      combination <- if (choose(n,k)==nrow(s)) 0 else sd(y)/sqrt(length(y))
      blocked <- jackknife(vapply(unique(sim$blocks),function(b) {
        contains <- rowSums(matrix(sim$blocks[s]==b,nrow(s))) > 0
        mean_valid(y[!contains])
      },numeric(1)))
      data.table(combination_se=combination, spatial_block_sensitivity=blocked)
    })
    data.table(k=k,process=processes,accumulated=colMeans(v),direct=colMeans(original),
               cumulative_absolute=colMeans(gross[[k]])) |>
      cbind(rbindlist(uncertainty)) |> merge(true_contrast,by="process")
  }))
  tempo <- rbindlist(temporal)
  end <- evaluate(plan,sets[[1]],alloc)
  pair <- evaluate(plan,sets_for(n,2L,choose(n,2),1L))
  ecological <- data.table(scenario,replicate,n,years=T,species=dim(x)[3],
    mean_annual_local=mean(tempo$local_temporal),mean_annual_latent_local=mean(tempo$latent_local),
    endpoint_local=mean_valid(end$after[,6]), exact_pair_change=mean_valid(pair$after[,2]-pair$before[,2]),
    occupancy_event_lag1=cor(head(tempo$occupancy_events,-1),tail(tempo$occupancy_events,-1)))
  exact <- a[k==1]
  if (scenario != "detection") stopifnot(max(abs(exact$estimate-exact$truth)) < 1e-9)
  switch(scenario,
    expansion=stopifnot(sum(exact$truth[exact$process=="expansion"])>0),
    contraction=stopifnot(sum(exact$truth[exact$process=="contraction"])<0),
    opposing=stopifnot(any(exact$truth>0),any(exact$truth<0)),
    stable_redistribution=stopifnot(abs(ecological$exact_pair_change)<1e-12,ecological$mean_annual_local>.2),
    reversal=stopifnot(abs(ecological$exact_pair_change)<1e-12,sum(abs(exact$truth))>0),
    shocks=stopifnot(ecological$mean_annual_local>0),
    detection=stopifnot(ecological$mean_annual_local>.02,ecological$mean_annual_latent_local==0),
    null=stopifnot(ecological$mean_annual_local==0))
  list(years=a[, `:=`(scenario=scenario,replicate=replicate)],
       summary=summary[, `:=`(scenario=scenario,replicate=replicate,n=n,years=T)],
       temporal=tempo[, `:=`(scenario=scenario,replicate=replicate)], ecological=ecological)
}
