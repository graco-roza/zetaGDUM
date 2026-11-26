.local_retention_weights <- function(n, k) {
  probability <- matrix(1, n, k)
  if (k > 1L) {
    for (j in seq_len(k - 1L)) {
      probability[, j + 1L] <- probability[, j] * pmax(seq_len(n) - j, 0) / (n - j)
    }
  }
  probability
}

.local_retention_components <- function(x, n, k) {
  probability <- .local_retention_weights(n, k)
  full_weight <- n / (n - seq_len(n) + 1L)
  truncated_weight <- if (k == n) full_weight else rowSums(probability)
  site <- matrix(0, n, k + 2L)
  weight <- cbind(probability, truncated_weight, full_weight)

  for (j in seq_len(ncol(x))) {
    occupied <- which(x[, j] > 0)
    if (!length(occupied)) next
    sorted <- sort(x[occupied, j])
    levels <- rle(sorted)
    before <- c(0L, cumsum(levels$lengths))[seq_along(levels$values)]
    increments <- diff(c(0, levels$values))
    by_threshold <- weight[length(occupied) - before, , drop = FALSE] * increments
    # Cumulative threshold mass up to each focal site's own abundance.
    if (length(increments) > 1L) {
      for (i in 2:length(increments)) {
        by_threshold[i, ] <- by_threshold[i, ] + by_threshold[i - 1L, ]
      }
    }
    index <- match(x[occupied, j], levels$values)
    site[occupied, ] <- site[occupied, , drop = FALSE] + by_threshold[index, , drop = FALSE]
  }
  site
}

#' Exact site-level retention and retention loss
#'
#' A site's conditional zeta profile is expected shared composition among
#' distinct sets containing that site. Occurrence divides expected shared
#' richness by focal-site richness. Relative abundance uses within-site
#' proportions and species-wise minima. Both are survival profiles starting
#' at one for sites with defined composition.
#'
#' At each occupancy or abundance threshold with \eqn{m} eligible sites
#' including the focal site, the expected retained order count is
#' \eqn{n/(n-m+1)}. The count is averaged over the focal site's species or
#' its abundance-threshold increments. `T` is the sum of conditional zeta
#' through `max_order`. `T_full` uses the closed form through all orders.
#' \eqn{T_+=1+(T-1)/\zeta_{2|i}} is the characteristic retained order count
#' conditional on survival through order two. It is undefined when sharing
#' at order two is zero or order two is unavailable.
#'
#' @param x Numeric site-by-species abundances. Rows are stations and may
#'   carry unique names. Missing abundances are invalid. Relative abundance
#'   requires positive row totals.
#'   Occurrence return values are missing at sites with zero richness.
#' @param type Focal-richness-normalised `"occurrence"` or `"relative"`
#'   abundance.
#' @param max_order Optional maximum retained order. Defaults to all sites.
#' @return A list with `summary` (site, richness, site count, maximum order,
#'   `T`, `T_plus` and `T_full`), `profile` (site by order conditional
#'   zeta), and `loss` (site by loss order). Loss at order \eqn{k+1} is
#'   \eqn{\zeta_k-\zeta_{k+1}}. The last column is right-censored mass
#'   surviving through `max_order`.
#' @export
zeta_local_retention <- function(x, type = c("occurrence", "relative"), max_order = NULL) {
  type <- match.arg(type)
  if (anyNA(x)) stop("Omit unobserved sites and supply complete abundance rows.")
  x <- .conditional_community(x, type == "relative")
  n <- nrow(x)
  k <- if (is.null(max_order)) n else .zeta_order_integer(max_order, n)
  site <- if (is.null(rownames(x))) seq_len(n) else rownames(x)
  if (anyNA(site) || anyDuplicated(site) || any(!nzchar(site))) {
    stop("Site names must be unique and nonempty.")
  }
  richness <- rowSums(x > 0)
  if (type == "occurrence") x <- 1 * (x > 0)
  all <- .local_retention_components(x, n, k)
  if (type == "occurrence") all <- all / ifelse(richness > 0, richness, NA_real_)
  profile <- all[, seq_len(k), drop = FALSE]
  # Remove tiny floating-point departures at order one.
  profile[richness > 0, 1L] <- 1
  dimnames(profile) <- list(site, seq_len(k))
  loss <- if (k == 1L) profile else cbind(
    pmax(profile[, -k, drop = FALSE] - profile[, -1L, drop = FALSE], 0),
    profile[, k])
  dimnames(loss) <- list(site, c(seq_len(k - 1L) + 1L, paste0(">", k)))
  T <- all[, k + 1L]
  second <- if (k > 1L) profile[, 2L] else rep(NA_real_, n)
  tail <- rowSums(profile[, -1L, drop = FALSE])
  T_plus <- 1 + tail / ifelse(second > 0, second, NA_real_)
  summary <- data.frame(site, n, max_order = k, richness, T, T_plus,
                        T_full = all[, k + 2L], row.names = NULL)
  list(summary = summary, profile = profile, loss = loss)
}

.retention_means <- function(x, keys, columns) {
  identifiers <- lapply(x[keys], function(z) match(z, unique(z)))
  groups <- split(seq_len(nrow(x)), interaction(identifiers, drop = TRUE))
  combined <- lapply(groups, function(i) {
    info <- x[i[1L], keys, drop = FALSE]
    info$occasions <- length(i)
    for (column in columns) {
      values <- x[[column]][i]
      values <- values[is.finite(values)]
      info[[paste0(column, "_n")]] <- length(values)
      info[[paste0(column, "_mean")]] <- if (length(values)) mean(values) else NA_real_
      info[[paste0(column, "_sd")]] <- if (length(values) > 1L) stats::sd(values) else NA_real_
    }
    info
  })
  out <- do.call(rbind, combined)
  rownames(out) <- NULL
  out
}

#' Compare local retention between sampling occasions
#'
#' Each occasion is evaluated separately. Use stations observed on every
#' occasion for directly comparable communities, or declare a common
#' maximum order and compare truncated profiles of all observed stations.
#' Communities and site counts can differ in the truncated comparison.
#' Seasonal means use each station's observed values. Standard deviations
#' require at least two observed values. Full-order `T_full` is retained
#' by occasion but excluded from seasonal averages.
#'
#' @param occasions Nonempty named list of site-by-species matrices with
#'   unique, nonempty station row names. Each matrix contains observed
#'   stations only.
#' @param season Optional vector of seasonal labels in occasion order.
#'   Defaults to a common group.
#' @param comparison Either `"common_sites"` or `"truncated"`.
#' @param type Focal-richness-normalised occurrence or relative abundance.
#' @param max_order Optional maximum order for common-site comparisons.
#'   Required for truncated comparisons and at most the smallest site count.
#' @return A list with `occasions` (retention fits), `seasonal` and
#'   `seasonal_profile`, `seasonal_loss` (station-level summaries),
#'   `comparison` and
#'   `max_order`.
#' @export
zeta_temporal_retention <- function(occasions, season = NULL,
                                    comparison = c("common_sites", "truncated"),
                                    type = c("occurrence", "relative"), max_order = NULL) {
  type <- match.arg(type)
  comparison <- match.arg(comparison)
  labels <- names(occasions)
  if (!is.list(occasions) || !length(occasions) || is.null(labels) ||
      anyNA(labels) || anyDuplicated(labels) || any(!nzchar(labels))) {
    stop("Supply a nonempty list of uniquely named occasions.")
  }
  stations <- lapply(occasions, rownames)
  valid <- vapply(stations, function(i) length(i) > 0L && !anyNA(i) &&
                    !anyDuplicated(i) && all(nzchar(i)), logical(1))
  if (!all(valid)) stop("Each occasion needs unique and nonempty station row names.")
  if (is.null(season)) season <- rep("all", length(occasions))
  if (!is.atomic(season) || length(season) != length(occasions) ||
      anyNA(season) || any(!nzchar(season))) {
    stop("Supply one nonempty seasonal label per occasion.")
  }
  if (comparison == "common_sites") {
    common <- Reduce(intersect, stations)
    if (!length(common)) stop("The occasions have no stations in common.")
    occasions <- lapply(occasions, function(x) x[common, , drop = FALSE])
  } else if (is.null(max_order)) {
    stop("Declare a common maximum order for truncated comparisons.")
  }
  smallest <- min(vapply(occasions, nrow, integer(1)))
  k <- if (is.null(max_order)) smallest else .zeta_order_integer(max_order, smallest)
  fits <- lapply(occasions, zeta_local_retention, type = type, max_order = k)
  totals <- do.call(rbind, Map(function(s, a) data.frame(a$summary, season = s),
                               as.character(season), fits))
  seasonal_curve <- function(field, order, value) {
    curves <- do.call(rbind, Map(function(s, a) {
      z <- a[[field]]
      data.frame(site = rep(a$summary$site, k), order = rep(order, each = nrow(z)),
                 measure = as.vector(z), season = s)
    }, as.character(season), fits))
    names(curves)[names(curves) == "measure"] <- value
    .retention_means(curves, c("site", "season", "order"), value)
  }
  list(occasions = fits, seasonal = .retention_means(totals, c("site", "season"), c("T", "T_plus")),
       seasonal_profile = seasonal_curve("profile", seq_len(k), "zeta"),
       seasonal_loss = seasonal_curve("loss", colnames(fits[[1L]]$loss), "mass"),
       comparison = comparison, max_order = k)
}
