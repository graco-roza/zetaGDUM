source("docs/order_information/design.R")
out <- lapply(seq_len(12L),function(i) {
  a <- community(150L,"narrow",5898L+i)
  b <- community(1800L,"narrow",4851L+i)
  set.seed(9900L+i)
  shuffled_a <- a
  shuffled_b <- b
  shuffled_a$env <- a$env[sample.int(nrow(a$env)),,drop=FALSE]
  shuffled_b$env <- b$env[sample.int(nrow(b$env)),,drop=FALSE]
  orders <- c(2L,6L,12L,24L)
  x <- zeta_order_information(a$x,a$env,orders=orders,method="sample",seed=715)
  y <- zeta_order_information(a$x,shuffled_a$env,orders=orders,method="sample",seed=715)
  order <- lapply(orders,function(k) {
    measure <- function(c,d) {
      s <- build(c,k,2000L,"relative",5534L+k)
      t <- build(d,k,4000L,"relative",5554L+k)
      evaluate_game(fit_game(s),t)[1L,"D2"]
    }
    data.frame(replicate=i,order=k,aligned=measure(a,b),shuffled=measure(shuffled_a,shuffled_b))
  })
  cbind(do.call(rbind,order),variance=x$response_variance,
    difference_in_profile_variance=x$response_variance-y$response_variance,
    difference_in_spread=x$spread_variance-y$spread_variance,
    difference_in_position=x$position_variance-y$position_variance,
    sampled_condition=x$design_condition,shuffled_sampled_condition=y$design_condition)
})
write_result(do.call(rbind,out),file.path(output,"permutation_control.csv"))
