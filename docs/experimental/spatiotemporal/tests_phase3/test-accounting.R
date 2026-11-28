test_that("positive-minimum partition agrees with direct relative Bray-Curtis", {
  a <- rbind(c(3,0,3,2,0), c(1,2,3,4,0), c(2,5,0,1,0))
  q <- a / rowSums(a)
  for (k in 2:3) {
    x <- a[seq_len(k), ]
    z <- st3_abundance_partition(x)
    expect_equal(unname(z$diversity["abundance"]), 1 - sum(apply(q[seq_len(k), ], 2, min)))
    expect_equal(sum(z$diversity[c("occupancy","numerical")]), unname(z$diversity["abundance"]))
    expect_true(all(z$species$occupancy >= 0 & z$species$numerical >= -1e-14))
    expect_equal(st3_abundance_partition(x * seq_len(k))$diversity, z$diversity)
  }
  b <- st3_abundance_partition(a[1:2, ])$species
  expect_equal(b$occupancy + b$numerical, abs(q[1, ] - q[2, ]) / 2)
  expect_equal(b$occupancy[b$shared], rep(0, sum(b$shared)))
  expect_equal(b$numerical[!b$shared], rep(0, sum(!b$shared)))
  expect_gt(st3_abundance_partition(a)$species$numerical[3], 0)
  expect_equal(st3_abundance_partition(rbind(c(1,2),c(2,1)))$diversity[["occupancy"]], 0)
  expect_equal(st3_abundance_partition(rbind(c(1,0,1),c(0,1,1)))$diversity[["numerical"]], 0)
  expect_equal(unname(st3_abundance_partition(rbind(c(1,2),c(1,2)))$diversity), rep(0,3))
})

test_that("native zeta and spatial indices are preserved, including Baselga at order two", {
  source("../phase2_load.R", local = TRUE)
  aux <- new.env()
  st2_load("..", aux)
  x <- array(c(3,2,0,1, 3,1,1,2, 1,0,2,1, 0,1,2,1, 2,3,2,1, 1,1,2,2), c(4,2,3))
  sets <- t(combn(4,2))
  a <- matrix(x[,1,],4)
  b <- matrix(x[,2,],4)
  plan <- st3_plan(x,1,2)
  z <- st3_shapley(plan,sets)
  pair <- zetaGDUM::pair_partition(a)
  expect_equal(unname(z$before[,"incidence"]), pair$sorensen)
  expect_equal(unname(z$before[,"turnover"]), pair$simpson)
  expect_equal(unname(z$before[,"nestedness"]), pair$nestedness)
  expect_equal(mean(z$after[,"zeta_abundance"]), zetaGDUM::zeta_curve(b,2,"abundance")$zeta[2])
  y <- array(rep(x,2),c(4,4,3))
  zz <- st3_shapley(st3_plan(y,1:2,3:4),sets)
  surf <- aux$st2_zeta(x,2,2)
  for (type in c("presence","abundance","dominance")) {
    expect_equal(mean(zz$before[,paste0("zeta_",type)]),
                 subset(surf$surface,response==type&spatial_order==2&temporal_order==2)$zeta)
  }
})

test_that("allocations average over admissible sites, preserving endpoint matrices and process budgets", {
  a <- cbind(c(1,0,0,0,2,0), c(2,2,2,2,2,2), c(0,2,2,2,0,0))
  b <- cbind(c(0,1,1,1,2,0), c(2,3,1,2,2,2), c(4,0,0,1,0,0))
  x <- array(cbind(a,b),c(6,3,2))
  x <- aperm(x,c(1,3,2))
  plan <- st3_plan(x,1,2)
  al <- st3_allocations(plan,exact_limit=10)
  expect_true(al$exact)
  expect_equal(length(al$masks),6)
  q <- plan$ledger
  expect_equal(q$gains-q$losses,q$expansion-q$contraction)
  expect_equal(q$gains+q$losses,q$expansion+q$contraction+2*q$redistribution)
  assignments <- lapply(al$masks, function(mask) {
    expect_equal(st3_apply(plan,mask,0),plan$before)
    expect_equal(st3_apply(plan,mask,15),plan$after)
    for (j in seq_along(q$species)) {
      expect_equal(sum(mask[,j]==1),q$expansion[j])
      expect_equal(sum(mask[,j]==2),q$contraction[j])
      expect_equal(sum(mask[,j]==4),2*q$redistribution[j])
    }
    which(mask[,1]==1)
  })
  expect_equal(as.integer(table(unlist(assignments))),rep(4L,3))
  combo <- t(combn(6,2))
  result <- st3_shapley(plan,combo,retain_draws=TRUE)
  expect_lt(result$max_efficiency_error,1e-12)
  for (one in result$draws) {
    expect_equal(unname(apply(one,c(1,2),sum)[result$valid]),unname(result$delta[result$valid]))
  }
  expect_equal(result$contributions[, "incidence", "numerical"], rep(0,nrow(combo)), tolerance=1e-12)
  cross <- st3_crossed(result)
  expect_equal(unname(cross$total_by_process),unname(result$contributions[,"abundance",]))
  expect_equal(unname(apply(cross$contributions,1,sum)),unname(cross$delta))
  sampled <- st3_shapley(plan,combo,exact_limit=1,draws=40,seed=6)
  expect_identical(sampled,st3_shapley(plan,combo,exact_limit=1,draws=40,seed=6))
  expect_true(any(sampled$allocation_mc_se > 0,na.rm=TRUE))
  expect_equal(mean(result$contributions[,"zeta_presence","redistribution"]),0,tolerance=1e-12)
})

test_that("redistribution changes temporal composition with fixed mean zeta", {
  a <- cbind(1,c(1,1,0,0),c(0,0,1,1))
  b <- a[c(3,4,1,2),]
  p <- st3_plan(aperm(array(cbind(a,b),c(4,3,2)),c(1,3,2)),1,2)
  for(k in 2:4) {
    d <- st3_shapley(p,t(combn(4,k)))
    expect_equal(mean(d$delta[,"zeta_presence"]),0)
    expect_equal(unname(mean(d$delta[,"zeta_presence"])),sum(.st_choose_prob(colSums(b>0),k,4)-.st_choose_prob(colSums(a>0),k,4)))
    expect_true(all(d$delta[,"site_temporal_presence"]>0))
  }
  expect_true(all(st3_site_reorganisation(p)$incidence>0))
})

test_that("undefined intermediate subsets invalidate the complete affected decomposition", {
  a <- rbind(c(1,0),c(1,1),c(1,1))
  b <- rbind(c(0,1),c(0,1),c(1,1))
  p <- st3_plan(aperm(array(cbind(a,b),c(3,2,2)),c(1,3,2)),1,2)
  d <- st3_shapley(p,rbind(c(1,2),c(2,3)))
  expect_true(is.finite(d$delta[1,"abundance"]))
  expect_false(d$valid[1,"abundance"])
  expect_true(all(is.na(d$contributions[1,"abundance",])))
  expect_true(d$valid[2,"abundance"])
  expect_true(d$valid[1,"zeta_presence"])
})

test_that("positive and negative attribution summaries do not cancel prematurely", {
  x <- matrix(c(-2,2,1,-1),2,byrow=TRUE)
  y <- st3_signed_summary(x)
  expect_equal(y$by_combination$net,c(0,0))
  expect_equal(y$by_combination$gross,c(4,2))
  expect_equal(y$mean_absolute,c(1.5,1.5))
  expect_equal(y$absolute_mean,c(.5,.5))
  expect_equal(unname(y$importance),matrix(.5,2,2))
})
