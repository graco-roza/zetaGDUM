source("docs/phase5/common.R")
file <- "docs/phase5/runs/coupling_2.csv"
bench <- if (file.exists(file)) read.csv(file) else
  subset(read.csv("docs/phase5/results/coupling.csv.gz"), case == 2L)
n <- 70L
s <- simulate_metacommunity(n,58,breadth=.36,spatial=.55,abundance_mean=1.4,seed=1503)
s$comm <- cbind(s$comm,diag(1L*(rowSums(s$comm)==0)))
orders <- c(4L,8L)
geo <- coords_dist(s$coords)
x <- subset(bench, component=="relative" & effect=="spread" & comparison=="change" & budget==2400)
cost <- aggregate(cpu ~ sampler + coupled, x, median)
targ <- cost$cpu[cost$sampler=="uniform" & !cost$coupled]
cost$M <- pmax(300L,round(2400*targ/cost$cpu))
fit <- function(a) {
  d <- model_data(a,"relative","B")
  f <- fit5(d,"B")
  u <- if (f$converged && is.finite(f$D2)) ugdm_local_distinctiveness(s$comm,f,d$design,
    "relative",d$sampling_weights)
  setNames(effect5(d,f,s$env,u),effect_names)
}
set.seed(614)
rows <- list()
truth <- subset(bench,component=="relative" & comparison=="change" & r==1L &
                budget==2400 & sampler=="uniform" & !coupled)
truth <- setNames(truth$reference,truth$effect)
for (i in seq_len(nrow(cost))) for (r in 1:24) {
  a <- cost[i,]
  t0 <- sum(proc.time()[1:2])
  s1 <- sample_sets_coupled(s$comm,orders,a$M,a$sampler,coupled=a$coupled,exact_limit=0)
  d <- ugdm_path_data(s$comm,s$env,s1,geo)
  e <- fit(d[[2]])-fit(d[[1]])
  elapsed <- sum(proc.time()[1:2])-t0
  rows[[length(rows)+1L]] <- data.frame(sampler=a$sampler,coupled=a$coupled,M=a$M,r,
    effect=names(truth),estimate=unname(e[names(truth)]),reference=unname(truth),cpu=elapsed,target_cpu=targ)
}
write5(do.call(rbind,rows),"equal_time")
