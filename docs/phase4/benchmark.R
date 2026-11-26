source("docs/phase4/common.R")
args <- commandArgs(TRUE)
case <- as.integer(args[1L])
cfg <- data.frame(n = c(18, 70, 102), species = c(30, 65, 140),
  breadth = c(.90,.33,.49), abundance = c(1.8,1.1,1.4), low = c(3,4,8), high = c(4,8,15))
stopifnot(case %in% c(seq_len(nrow(cfg)), 4L))
if (case == 4L) {
  a <- read_manifest()[5L, ]
  d <- load_dd2(a, args[2L])
  c <- data.frame(n = a$sites, low = a$mid, high = a$high)
  s <- list(comm = d$x)
  env <- d$env[1L]
} else {
  c <- cfg[case, ]
  s <- simulate_metacommunity(c$n, c$species, breadth = c$breadth, abundance_mean = c$abundance,
                              spatial = .4, use_env = c(TRUE, FALSE), seed = 9212 + case)
  s$comm <- safe_community(s$comm)
  env <- s$env["env1"]
}
p <- phase4_input(s$comm, env)
spec <- internal(".limit_effect_spec")(env, p$basis, names(env)[1L])
samplers <- c("conventional", "species", "site", "stratified")
budgets <- c(500L, 2000L, 8000L)
k_grid <- c(c$low, c$high)
reps <- if (case == 4L) 12L else 18L
set.seed(3251 + case)
output <- list()
for (component in c("total", "relative")) {
  for (k in k_grid) {
    M_ref <- 60000L
    ref <- internal(".limit_draw")(p, k, M_ref, "conventional", .5, 60000L, 10000L)
    f_ref <- internal(".limit_fit")(p, ref, component, "constant", FALSE, spec, TRUE)
    truth <- f_ref$effects
    for (M in budgets) for (m in samplers) for (r in seq_len(reps)) {
      start <- sum(proc.time()[1:2])
      a <- internal(".limit_draw")(p, k, M, m, .5, 0L, max(4000L, M))
      f <- internal(".limit_fit")(p, a, component, "constant", FALSE, spec, TRUE)
      cpu <- sum(proc.time()[1:2]) - start
      d <- f$diagnostics
      output[[length(output) + 1L]] <- data.frame(case, n = c$n, k, component,
        M, sampler = m, r, effect = names(truth), point = unname(f$effects),
        reference = unname(truth), reference_exact = ref$exact,
        sharing = d[["sharing_estimate"]], shared_draws = d[["shared_draws"]],
        response_sd = d[["response_sd"]], exclusion_ess = d[["minimum_delete_ess"]],
        fallback = d[["sampler_fallback"]], regular = f$regular, cpu,
        mass_source = if (is.null(a$mass_source)) NA_character_ else a$mass_source,
        mass_draws = if (identical(a$mass_source, "independent_pilot")) max(4000L, M) else 0,
        proposals = if (is.null(a$proposals)) 0 else sum(a$proposals))
    }
  }
}
phase4_write(do.call(rbind, output), "comparison", if (case == 4L) "fish" else paste0("sim", case))
