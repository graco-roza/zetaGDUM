# 02_orders.R: zeta_order_information() for the three communities (candidate informative orders)
source("lib6.R"); ORD <- c(2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30)
res <- list(); for (g in names(D)) { d <- D[[g]]; t0 <- Sys.time()
  oi <- zeta_order_information(d$x, get_env(d), get_xy(d), orders = ORD, response = "all", seed = 1, lonlat = TRUE, combinations = 1000L)
  cat(g, "elapsed", format(Sys.time() - t0), "\n"); saveRDS(oi, sprintf("out/orderinfo_%s.rds", g)); res[[g]] <- oi }
print(names(res$PP)); str(res$PP, max.level = 1, give.attr = FALSE); print(head(as.data.frame(res$PP), 15))
