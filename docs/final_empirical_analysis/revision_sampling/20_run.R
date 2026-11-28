# 20_run.R: run the sampling-requirement simulation. Usage (from final/, ZG_LIB6 = principal package library): SIM=meta|known NCORES=9 Rscript revision_sampling/20_run.R
# Series (unique units; a unit = site draw at one order, evaluated with BOTH combination samplers):
#   S  sites:          n in {100,250,500,1000,2500,5000}, 4,000 combinations, orders 2,3,5,10,20, 3 environments, 6 replicates
#   C  combinations:   n = 500, C in {250,1000,4000,16000}, orders 2,3,5,10,20, 3 environments, 6 replicates
#   F  factorial:      n x C (6 x 4) at order 5, 3 environments, 3 replicates (separates sites from combinations in two dimensions)
SIM <- Sys.getenv("SIM", "meta"); NC <- as.integer(Sys.getenv("NCORES", "9")); source(if (SIM == "meta") "revision_sampling/10_meta_rep.R" else "revision_sampling/11_known_rep.R")
NS <- c(100L, 250L, 500L, 1000L, 2500L, 5000L); CS <- c(250L, 1000L, 4000L, 16000L)
U <- rbind(cbind(expand.grid(n = NS, k = KS, env = ENVS, rep = 1:6, C = 4000L, stringsAsFactors = FALSE), series = "S"), cbind(expand.grid(n = 500L, k = KS, env = ENVS, rep = 1:6, C = CS, stringsAsFactors = FALSE), series = "C"),
  cbind(expand.grid(n = NS, k = 5L, env = ENVS, rep = 1:3, C = CS, stringsAsFactors = FALSE), series = "F"))
key <- function(u) paste(u$n, u$k, u$env, u$rep, u$C, sep = "_"); U$key <- key(U)
UU <- U[!duplicated(U$key), ]; UU$in_S <- UU$key %in% U$key[U$series == "S"]; UU$in_C <- UU$key %in% U$key[U$series == "C"]; UU$in_F <- UU$key %in% U$key[U$series == "F"]
UU$cost <- UU$n^.3 * UU$C^.7 * UU$k^.5; UU <- UU[order(-UU$cost), ]; D <- file.path("revision_sampling/logs", paste0("units_", SIM)); dir.create(D, FALSE, TRUE); cat(nrow(UU), "unique units\n")
one <- function(i) { u <- UU[i, ]; f <- file.path(D, paste0(u$key, ".rds")); if (file.exists(f)) return(invisible(NULL)); t0 <- Sys.time()
  r <- tryCatch((if (SIM == "meta") run_meta else run_known)(u$n, u$env, u$rep, u$k, u$C, paste0(SIM, "_", u$key)), error = function(e) { message("UNIT ERR ", u$key, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) { r$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); saveRDS(r, f) }; invisible(NULL) }
invisible(mclapply(seq_len(nrow(UU)), one, mc.cores = NC, mc.preschedule = FALSE))
fl <- list.files(D, full.names = TRUE); R <- lapply(fl, readRDS); b <- function(nm) bind_rows(lapply(R, `[[`, nm))
attach_series <- function(x) { x$key <- paste(x$n_sites, x$order, x$env, x$rep, x$n_combos_target, sep = "_"); x <- left_join(x, UU[, c("key", "in_S", "in_C", "in_F")], by = "key"); x$key <- NULL; x }
dir.create("revision_sampling/tables", FALSE, TRUE); for (nm in c("fit", "sup", "con")) write.csv(attach_series(b(nm)), sprintf("revision_sampling/tables/sim_%s_%s.csv", SIM, nm), row.names = FALSE)
cat("done", nrow(UU), "units;", length(fl), "completed\n")
