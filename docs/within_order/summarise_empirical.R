out <- "docs/within_order/results"
manifest <- read.csv("docs/manuscript_validation/communities.csv")
merge_outputs <- function(kind) do.call(rbind, lapply(manifest$id, function(id)
  read.csv(file.path(out, paste0("empirical_", kind, "_", id, ".csv")))))
d <- merge_outputs("diagnostics")
a <- subset(merge_outputs("importance"), group == "predictor")
lead <- function(labels, value) {
  if (all(is.finite(value))) labels[which.max(value)] else NA_character_
}
key <- c("id", "community", "response", "order")
groups <- split(a, a[key], drop = TRUE)
p <- do.call(rbind, lapply(groups, function(x) {
  fit <- lead(x$term, x$absolute)
  common <- lead(x$term, x$common_assessment_absolute)
  second <- lead(x$term, x$second_pool_absolute)
  data.frame(x[1L, key], leading = fit, leading_on_assessment = common,
    leading_on_second_pool = second, same_lead_when_reintegrated = common == second,
    maximum_permutation_mc_se = if (all(is.finite(x$allocation_mc_se)))
      max(x$allocation_mc_se) else NA_real_)
}))
report <- merge(d, p, by = key, all.x = TRUE, sort = FALSE)
report <- report[order(match(report$id, manifest$id), report$response, report$order), ]
write.csv(report, file.path(out, "empirical_overview.csv"), row.names = FALSE)
