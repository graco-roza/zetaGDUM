# Collect validation-only independent chains, leaving the Phase 5 archives intact.
folder <- "docs/phase5_refinement/runs"
keys <- c("null_1_1_900", "null_2_2_900", "weak_2_2_900", "weak_3_3_900")
out <- "docs/interaction_refinement/results"
collect <- function(prefix) {
  files <- file.path(folder, paste0(prefix, "_long_reference_", keys, ".csv"))
  if (!all(file.exists(files))) stop("Run the four additional independent chains before collecting.")
  do.call(rbind, lapply(files, read.csv))
}
write.csv(collect("posterior"), file.path(out, "bayes_additional_jel.csv"), row.names = FALSE)
write.csv(collect("postdiag"), file.path(out, "bayes_additional_jel_diagnostics.csv"), row.names = FALSE)
