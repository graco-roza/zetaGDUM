source("docs/order_information/design.R")
name <- commandArgs(TRUE)[1]
id <- scenario_index(name)
M <- if (length(commandArgs(TRUE)) >= 2) as.integer(commandArgs(TRUE)[2]) else 16000L
a <- lapply(1:3, function(j) community(2200L, name, 14122L+id*100L+j))
rows <- lapply(order_grid(name), function(k) {
  out <- lapply(c("total","relative"), function(response) {
    pools <- lapply(1:3,function(j) build(a[[j]],k,M,response,11801L+j*100L+k))
    do.call(rbind,lapply(1:3,function(j) {
      game <- fit_game(pools[[j]])
      est <- evaluate_game(game,pools[[j %% 3L+1L]])
      data.frame(name,order=k,response,reference=j,sets=M,
        converged=all(vapply(game$fits,function(z) z$converged,logical(1))),est)
    }))
  })
  do.call(rbind,out)
})
write_result(do.call(rbind,rows), file.path(output,paste0("reference_",name,".csv")))
