source("docs/final_gdum/design.R")
rows <- do.call(rbind,lapply(mechanisms$name,function(s) {
  do.call(rbind,lapply(1:12,function(j) {
    seed <- 380881L+scenario_id(s)*10000L+31L*j
    a <- generate(110L,s,seed)
    z <- zeta_order_information(a$x,a$e,orders=2:12,response="all",method="sample",
      combinations=320,pools=2,seed=41L+seed)
    data.frame(scenario=s,replicate=j,z)
  }))
}))
atomic_csv(rows,"simulation_profiles.csv.gz")
