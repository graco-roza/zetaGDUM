library(data.table)

mean_valid <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
jackknife <- function(x) {
  if (length(x) < 3 || any(!is.finite(x))) return(NA_real_)
  sqrt((length(x)-1) * mean((x - mean(x))^2))
}
summarize_phi <- function(z, sets, blocks = NULL, exhaustive = FALSE, support = NULL) {
  summaries <- lapply(seq_along(responses), function(r) {
    valid <- apply(is.finite(z$phi[,r,,drop=FALSE]),1,all)
    if (!is.null(support)) valid <- valid & support[,r]
    do.call(rbind, lapply(seq_along(processes), function(p) {
      v <- z$phi[,r,p]
      net <- mean_valid(v[valid])
      h <- mean_valid(pmax(-v[valid],0))
      d <- mean_valid(pmax(v[valid],0))
      draw <- vapply(seq_len(dim(z$draws)[4]), function(i) mean_valid(z$draws[valid,r,p,i]), numeric(1))
      mc <- if (any(valid) && z$exact) 0 else sd(draw)/sqrt(length(draw))
      combination <- if (any(valid) && exhaustive) 0 else sd(v[valid])/sqrt(sum(valid))
      spatial <- if (is.null(blocks)) NA_real_ else jackknife(vapply(sort(unique(blocks)), function(group) {
        contains <- rowSums(matrix(blocks[sets] == group,nrow(sets))) > 0
        mean_valid(v[valid & !contains])
      }, numeric(1)))
      data.frame(response = responses[r], process = processes[p], net, homogenising = h,
                 differentiating = d, gross = h+d, supported = sum(valid), sets = nrow(sets),
                 allocation_se = mc, combination_se = combination, spatial_block_se = spatial,
                 endpoint_before = mean_valid(z$before[valid,r]), endpoint_after = mean_valid(z$after[valid,r]))
    }))
  })
  result <- rbindlist(summaries)
  result[, `:=`(homogenising_percent = 100*.st_ratio(homogenising,sum(homogenising)),
                 differentiating_percent = 100*.st_ratio(differentiating,sum(differentiating))), by=response]
  result
}

analyse_trajectory <- function(panel, blocks = NULL, budget = 512L, draws = 32L,
                               orders = 1:6, windows = TRUE, keep_sets = FALSE, exhaustive_limit = 100000) {
  n <- dim(panel$x)[1]
  years <- panel$times
  ny <- length(years)
  stopifnot(all(diff(years) == 1), ny >= 2, max(orders) <= n)
  sets <- lapply(orders, function(k) sets_for(n,k,budget,1504L+k))
  names(sets) <- orders
  exact <- vapply(seq_along(sets), function(i) choose(n,orders[i]) == nrow(sets[[i]]) &&
                    choose(n,orders[i]) <= budget, logical(1))
  run <- function(before, after, seed) {
    plan <- st3_plan(panel,before,after)
    alloc <- st3_allocations(plan,draws=draws,seed=seed)
    lapply(sets,evaluate,plan=plan,allocations=alloc)
  }
  annual <- lapply(seq_len(ny-1L), function(i) run(i,i+1L,5109+i))
  endpoint <- run(1L,ny,51100)
  summaries <- cumulative <- comparisons <- endpoints <- raw_sets <- spatial_annual <- list()
  errors <- c()
  for (i in seq_along(sets)) {
    k <- orders[i]
    visit <- lapply(annual, `[[`, i)
    direct <- endpoint[[i]]
    common <- Reduce(`&`, lapply(c(visit,list(direct)), function(z) apply(is.finite(z$phi),c(1,2),all)))
    for (r in 1:4) {
      summed <- Reduce(`+`,lapply(visit,function(z) z$after[,r]-z$before[,r]))
      errors <- c(errors, abs(summed[common[,r]] - (direct$after[,r]-direct$before[,r])[common[,r]]))
    }
    summaries[[i]] <- rbindlist(lapply(seq_len(ny-1L), function(j) {
      s <- summarize_phi(visit[[j]],sets[[i]],blocks,exact[i])
      s[, `:=`(k=k, initial=years[j], final=years[j+1L])]
      s
    }))
    csum <- rbindlist(lapply(seq_len(ny-1L), function(j)
      summarize_phi(visit[[j]],sets[[i]],blocks,exact[i],common)[,year:=years[j+1L]]))
    setorder(csum,year)
    csum[, c("cumulative_net","cumulative_homogenising","cumulative_differentiating") :=
          lapply(.SD,cumsum), by=.(response,process), .SDcols=c("net","homogenising","differentiating")]
    csum[,k:=k]
    cumulative[[i]] <- csum
    comparisons[[i]] <- merge(
      summarize_phi(direct,sets[[i]],blocks,exact[i],common)[, .(response,process,direct=net,
        direct_endpoint_reorganisation=endpoint_after,common_sets=supported,direct_allocation_se=allocation_se)],
      csum[year==max(year), .(response,process,accumulated=cumulative_net,
        cum_homogenising=cumulative_homogenising,cum_differentiating=cumulative_differentiating)], by=c("response","process"))
    comparisons[[i]][, `:=`(k=k,path_difference=accumulated-direct)]
    comparisons[[i]][response %in% responses[5:6], path_difference:=NA_real_]
    all_sets <- if (choose(n,k) <= exhaustive_limit) sets_for(n,k,exhaustive_limit,1) else NULL
    plan <- st3_plan(panel,1L,ny)
    all_endpoint <- if (!is.null(all_sets)) evaluate(plan,all_sets) else NULL
    endpoints[[i]] <- rbindlist(lapply(seq_along(responses), function(r) {
      baseline <- direct$before[,r]; end <- direct$after[,r]; use <- common[,r]
      direct_valid <- is.finite(baseline) & is.finite(end)
      observed <- if (r==1) c(sum(exact_zeta(panel$x[,1L,,drop=FALSE],k)),
                              sum(exact_zeta(panel$x[,ny,,drop=FALSE],k))) else if (!is.null(all_endpoint)) {
        col <- is.finite(all_endpoint$before[,r]) & is.finite(all_endpoint$after[,r])
        c(mean_valid(all_endpoint$before[col,r]),mean_valid(all_endpoint$after[col,r]))
      } else c(NA_real_,NA_real_)
      data.table(k=k,response=responses[r],exact_before=observed[1],exact_after=observed[2],
        exact_change=diff(observed),exact_endpoint_sets=if (!is.null(all_endpoint)) sum(is.finite(all_endpoint$before[,r]) &
          is.finite(all_endpoint$after[,r])) else if (r==1) choose(n,k) else NA_real_,
        sampled_endpoint_sets=sum(direct_valid), common_sets=sum(use), total_sets=choose(n,k),
        sampled_endpoint_change=mean_valid((end-baseline)[direct_valid]),
        common_support_change=mean_valid((end-baseline)[use]))
    }))
    if (keep_sets) {
      raw_sets[[i]] <- rbindlist(lapply(seq_along(visit), function(j) {
        q <- as.data.table(as.table(visit[[j]]$phi))
        setnames(q,c("set","response","process","contribution"))
        q[, `:=`(k=k, initial=years[j], final=years[j+1L],
                  sites=rep(apply(matrix(dimnames(panel$x)[[1]][sets[[i]]],ncol=k),1,paste,collapse=";"),18))]
        q
      }))
    }
    spatial_annual[[i]] <- rbindlist(lapply(seq_len(ny), function(j) {
      obs <- if (j==1) visit[[1]]$before else visit[[j-1]]$after
      vals <- if (!is.null(all_sets) && !exact[i] && j<ny) {
        evaluate(st3_plan(panel,j,j+1L),all_sets)$before
      } else if (!is.null(all_endpoint) && j==ny) all_endpoint$after else obs
      rbindlist(lapply(seq_len(4), function(r) data.table(year=years[j], k=k,
        response=responses[r], value=if (r==1) sum(exact_zeta(panel$x[,j,,drop=FALSE],k)) else mean_valid(vals[,r]),
        sets=if (r==1) choose(n,k) else nrow(vals), method=if (r==1) "analytical" else
          if (!is.null(all_sets)) "exhaustive" else "sampled")))
    }))
  }
  stopifnot(max(c(errors,0)) < 1e-11)
  diagnostic <- grid <- weights <- list()
  if (windows && ny >= 10) {
    for (m in 1:5) for (j in seq_len(ny-9L)) {
      before <- j+seq_len(m)-1L
      after <- before+5L
      z <- run(before,after,6104+j+50*m)
      for (i in seq_along(orders)) {
        k <- orders[i]
        key <- paste(k,m,j)
        v <- z[[i]]
        a <- exact_zeta(panel$x[,before,,drop=FALSE],k)
        b <- exact_zeta(panel$x[,after,,drop=FALSE],k)
        grid[[key]] <- summarize_phi(v,sets[[i]],blocks,exact[i])[, `:=`(k=k,m=m,start=years[j],lag=5L)]
        weights[[key]] <- data.table(k=k,m=m,start=years[j],taxon=dimnames(panel$x)[[3]],
                                      before=.st_ratio(a,sum(a)), after=.st_ratio(b,sum(b)))
        diagnostic[[key]] <- rbindlist(lapply(seq_along(responses),function(r) {
          finite <- apply(is.finite(v$phi[,r,,drop=FALSE]),1,all)
          data.table(k=k,m=m,start=years[j],lag=5L,duration=m-1L,response=responses[r],
            exact_zeta_before=sum(a),exact_zeta_after=sum(b),
            effective_taxa_before=.st_ratio(sum(a)^2,sum(a*a)),
            effective_taxa_after=.st_ratio(sum(b)^2,sum(b*b)),
            zero_sharing_before=mean(v$before[,1]==0),zero_sharing_after=mean(v$after[,1]==0),
            valid_fraction=mean(finite),fractions_at_one=mean_valid(v$before[finite,r]==1),
            nonzero_delta_fraction=mean_valid(abs(v$after[finite,r]-v$before[finite,r]) > 1e-12),
            allocations=dim(v$draws)[4],combination_sets=nrow(sets[[i]]),exact_allocation=v$exact,
            exact_combinations=exact[i])
        }))
      }
    }
  }
  dg <- rbindlist(diagnostic)
  if (nrow(dg)) {
    dg[, `:=`(spatial_retention_before=.st_ratio(exact_zeta_before,exact_zeta_before[match(k-1,k)]),
              spatial_retention_after=.st_ratio(exact_zeta_after,exact_zeta_after[match(k-1,k)])),
       by=.(m,start,response)]
    dg[, `:=`(temporal_retention_before=.st_ratio(exact_zeta_before,exact_zeta_before[match(m-1,m)]),
              temporal_retention_after=.st_ratio(exact_zeta_after,exact_zeta_after[match(m-1,m)])),
       by=.(k,start,response)]
  }
  list(annual=rbindlist(summaries), cumulative=rbindlist(cumulative), comparison=rbindlist(comparisons),
       endpoints=rbindlist(endpoints), spatial=rbindlist(spatial_annual), window=rbindlist(grid),
       diagnostics=dg, weights=rbindlist(weights), set_contributions=rbindlist(raw_sets),
       telescope_error=max(c(0,errors)),
       efficiency_error=max(sapply(c(unlist(annual,recursive=FALSE),endpoint), `[[`, "error")))
}
